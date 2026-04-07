#include <fairmq/Device.h>
#include <fairmq/runDevice.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <sstream>

// ROOT Headers
#include "TSystem.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "THttpServer.h"

// Project Headers from spadi/nestdaq-user-impl
#include "TimeFrameHeader.h"
#include "SubTimeFrameHeader.h"
#include "HeartbeatFrameHeader.h"
#include "FilterHeader.h"
#include "KTimer.cxx"

// Decoder Headers from nestdaq-user-impl (as used in decode_slice.cc)
#include "AmQStrTdcData.h"

namespace bpo = boost::program_options;

struct OnlineSlicerNode : fair::mq::Device
{
    struct OptionKey {
        static constexpr std::string_view InputChannelName {"in-chan-name"};
    };

    enum {
        kAmQData,
        kAmQHeartbeat,
        kAmQHeartbeat2nd
    };

    enum {
        kHRTDC = 5,
        kLRTDC = 6
    };

    OnlineSlicerNode() {}

    void InitTask() override
    {
        using opt = OptionKey;
        fInputChannelName = fConfig->GetValue<std::string>(opt::InputChannelName.data());
        LOG(info) << "OnlineSlicerNode: InitTask Started. Input: " << fInputChannelName;

        gROOT->SetBatch(kTRUE);
        if (!gApplication) {
            static int argc = 1;
            static char* argv[] = {(char*)"OnlineSlicerNode", nullptr};
            new TApplication("OnlineSlicerNode", &argc, argv);
        }

        if (!fServer) {
            fServer = new THttpServer("http:8889"); 
            fServer->SetReadOnly(kTRUE);
            LOG(info) << "THttpServer started on port 8889";
        }

        fH1TrigInterval = new TH1F("h1_trig_interval", "Trigger Interval;Time [us];Counts", 1000, 0, 100);
        fH1TrigTime     = new TH1F("h1_trig_time", "Trigger Time distribution;Time [s];Counts", 1000, 0, 1.0);
        
        if(fServer) {
            fServer->Register("/Summary", fH1TrigInterval);
            fServer->Register("/Summary", fH1TrigTime);
            fServer->SetItemField("/", "_monitoring", "1000");
        }

        fDrawTimer.SetDuration(100); 
    }

    int CheckAMANEQHeader(AmQStrTdc::Data::v1::Bits bits) {
        auto head = bits.head;
        if (head == AmQStrTdc::Data::v1::HeadTypes::Data) return kAmQData;
        if (head == AmQStrTdc::Data::v1::HeadTypes::Heartbeat) return kAmQHeartbeat;
        if (head == AmQStrTdc::Data::v1::HeadTypes::Heartbeat2nd) return kAmQHeartbeat2nd;
        return -1;
    }

    bool ConditionalRun() override
    {
        try {
            fair::mq::Parts parts;
            if (Receive(parts, fInputChannelName, 0, 100) <= 0) return true;

            size_t total_size = 0;
            for (const auto& msg : parts) total_size += msg->GetSize();
            if (total_size == 0) return true;

            size_t n_words = (total_size + 7) / 8;
            std::vector<uint64_t> merged_buffer(n_words, 0);
            uint8_t* dest_ptr = reinterpret_cast<uint8_t*>(merged_buffer.data());
            size_t offset = 0;
            for (const auto& msg : parts) {
                std::memcpy(dest_ptr + offset, msg->GetData(), msg->GetSize());
                offset += msg->GetSize();
            }

            uint64_t* data = merged_buffer.data();
            uint64_t* end_ptr = data + n_words;
            uint64_t* current_ptr = data;

            std::vector<uint64_t> trg_times;
            size_t current_slice_idx = -1; // Index 0 when first SLICE seen
            uint32_t fem_type = 0;
            uint32_t fem_id = 0;

            while (current_ptr < end_ptr) {
                uint64_t word = *current_ptr;
                uint64_t magic = word & 0x00FFFFFFFFFFFFFF; 

                if (magic == TimeFrame::v1::MAGIC || magic == TimeFrame::v0::MAGIC) {
                    auto* h_tf = reinterpret_cast<TimeFrame::v1::Header*>(current_ptr);
                    if (h_tf->type == TimeFrame::META) {
                        trg_times.clear();
                        current_slice_idx = -1; 
                    } else if (h_tf->type == TimeFrame::SLICE) {
                        current_slice_idx++;
                    }
                    current_ptr += 3;
                    continue;
                }

                if (magic == Filter::v1::MAGIC) {
                    auto* h_flt = reinterpret_cast<Filter::v1::Header*>(current_ptr);
                    uint16_t num_trigs = h_flt->numTrigs;
                    current_ptr += 7; 
                    
                    if (current_ptr < end_ptr && ((*current_ptr & 0x00FFFFFFFFFFFFFF) == Filter::TDC_MAGIC)) {
                        current_ptr += 2; 
                        for (int i = 0; i < num_trigs; ++i) {
                            if (current_ptr >= end_ptr) break;
                            Filter::v1::TrgTime tt;
                            std::memcpy(&tt, current_ptr, 8);
                            uint64_t t_ps = static_cast<uint64_t>(tt.time) * 4000;
                            trg_times.push_back(t_ps);
                            
                            if (fLastTrgTime > 0) fH1TrigInterval->Fill((t_ps - fLastTrgTime) * 1e-6);
                            fLastTrigTime = t_ps;
                            current_ptr++;
                        }
                    }
                    continue;
                }

                if (magic == SubTimeFrame::v1::MAGIC || magic == SubTimeFrame::v0::MAGIC) {
                    auto* h_stf = reinterpret_cast<SubTimeFrame::v1::Header*>(current_ptr);
                    fem_type = h_stf->femType;
                    fem_id = h_stf->femId;
                    current_ptr += 6;
                    continue;
                }

                if (magic == HeartbeatFrame::MAGIC) {
                    auto* h_hbf = reinterpret_cast<HeartbeatFrame::Header*>(current_ptr);
                    uint32_t length = h_hbf->length / 8 - 2;
                    current_ptr += 2;

                    uint64_t current_trig = (current_slice_idx < trg_times.size()) ? trg_times[current_slice_idx] : 0;

                    for (uint32_t i = 0; i < length; ++i) {
                        if (current_ptr >= end_ptr) break;
                        AmQStrTdc::Data::v1::Bits h_amq;
                        std::memcpy(&h_amq, current_ptr, 8);
                        
                        auto flag = CheckAMANEQHeader(h_amq);
                        if (flag == kAmQData) {
                            uint64_t hit_ps = 0;
                            uint64_t tot = 0;
                            uint64_t ch = 0;
                            bool is_hr = (fem_type == kHRTDC);

                            if (fem_type == kLRTDC) {
                                hit_ps = static_cast<uint64_t>(h_amq.tdc) * 1024;
                                tot = h_amq.tot;
                                ch = h_amq.ch;
                            } else if (fem_type == kHRTDC) {
                                hit_ps = h_amq.hrtdc;
                                tot = h_amq.hrtot;
                                ch = h_amq.hrch;
                            }
                            
                            if (fem_id != 0) {
                                CheckAndCreateHistograms(fem_id, is_hr);
                                fMapHitPattern[fem_id]->Fill(ch);
                                if (current_trig > 0) {
                                    long long diff = (long long)hit_ps - (long long)current_trig;
                                    if (ch < fMapTDC[fem_id].size()) {
                                        fMapTDC[fem_id][ch]->Fill(is_hr ? diff : diff / 1000.0);
                                        fMapTOT[fem_id][ch]->Fill(tot);
                                    }
                                }
                            }
                        }
                        current_ptr++;
                    }
                    continue;
                }
                current_ptr++;
            }
            if (fDrawTimer.Check()) UpdateDisplay();
        } catch (const std::exception& e) {
            LOG(warn) << "OnlineSlicerNode: " << e.what();
        }
        return true;
    }

    void CheckAndCreateHistograms(uint32_t fem_id, bool is_hr) {
        if (fMapHitPattern.count(fem_id) > 0) return;

        std::string ip_str = ToIPAddress(fem_id);
        int max_ch = is_hr ? 64 : 128;
        int max_pages = (max_ch + 15) / 16;
        LOG(info) << "Registering FEM ID: " << fem_id << " (" << ip_str << ")";

        fMapHitPattern[fem_id] = new TH1F(Form("h1_hitpat_%s", ip_str.c_str()), Form("FEM %s Hit Pattern", ip_str.c_str()), max_ch, 0, max_ch);
        fMapTDC[fem_id].resize(max_ch, nullptr);
        fMapTOT[fem_id].resize(max_ch, nullptr);

        for (int p = 0; p < max_pages; ++p) {
            auto* c_tdc = new TCanvas(Form("c_tdc_%s_p%d", ip_str.c_str(), p), Form("FEM %s TDC p%d", ip_str.c_str(), p), 1200, 800);
            c_tdc->Divide(4, 4);
            fCanvases[fem_id].push_back(c_tdc);
            fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tdc);
        }

        for (int i = 0; i < max_ch; ++i) {
            if (is_hr) {
                fMapTDC[fem_id][i] = new TH1F(Form("h1_tdc_%s_ch%d", ip_str.c_str(), i), "Rel Time [ps]", 2000, -100000, 100000);
                fMapTOT[fem_id][i] = new TH1F(Form("h1_tot_%s_ch%d", ip_str.c_str(), i), "TOT", 1000, 0, 200000);
            } else {
                fMapTDC[fem_id][i] = new TH1F(Form("h1_tdc_%s_ch%d", ip_str.c_str(), i), "Rel Time [ns]", 2000, -1000, 1000);
                fMapTOT[fem_id][i] = new TH1F(Form("h1_tot_%s_ch%d", ip_str.c_str(), i), "TOT", 1000, 0, 1000);
            }
            fServer->Register(Form("/FEM_%s/TDC", ip_str.c_str()), fMapTDC[fem_id][i]);
            fServer->Register(Form("/FEM_%s/TOT", ip_str.c_str()), fMapTOT[fem_id][i]);

            int page = i / 16;
            int pad = (i % 16) + 1;
            fCanvases[fem_id][page]->cd(pad);
            fMapTDC[fem_id][i]->Draw();
        }
        fServer->Register(Form("/FEM_%s", ip_str.c_str()), fMapHitPattern[fem_id]);
    }

    void UpdateDisplay() {
        for (auto const& [id, canv_vec] : fCanvases) {
            for (auto* c : canv_vec) {
                if (c) { 
                    for(int i=1; i<=16; ++i) if(c->GetPad(i)) c->GetPad(i)->Modified();
                    c->Update(); 
                }
            }
        }
        if (gSystem) gSystem->ProcessEvents();
    }

    std::string ToIPAddress(uint32_t fem_id) {
        uint8_t* ip = reinterpret_cast<uint8_t*>(&fem_id);
        return std::to_string(ip[3]) + "." + std::to_string(ip[2]) + "." + std::to_string(ip[1]) + "." + std::to_string(ip[0]);
    }

private:
    std::string fInputChannelName;
    THttpServer* fServer = nullptr;
    KTimer fDrawTimer;
    uint64_t fLastTrgTime = 0;

    TH1F *fH1TrigInterval, *fH1TrigTime;
    std::map<uint32_t, TH1F*> fMapHitPattern;
    std::map<uint32_t, std::vector<TH1F*>> fMapTDC, fMapTOT;
    std::map<uint32_t, std::vector<TCanvas*>> fCanvases;
};

void addCustomOptions(bpo::options_description& options)
{
    using opt = OnlineSlicerNode::OptionKey;
    options.add_options()
        (opt::InputChannelName.data(), bpo::value<std::string>()->default_value("in"), "Input channel name");
}

std::unique_ptr<fair::mq::Device> getDevice(fair::mq::ProgOptions& /*config*/)
{
    return std::make_unique<OnlineSlicerNode>();
}
