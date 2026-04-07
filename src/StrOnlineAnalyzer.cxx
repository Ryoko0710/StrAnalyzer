#include <fairmq/Device.h>
#include <fairmq/runDevice.h>

#include <iomanip>
#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <thread>

#include "TSystem.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "THttpServer.h"

// Headers from spadi/nestdaq-user-impl
#include "TimeFrameHeader.h"
#include "SubTimeFrameHeader.h"
#include "HeartbeatFrameHeader.h"
#include "FileSinkHeader.h"
#include "FileSinkTrailer.h"
#include "FilterHeader.h"
#include "KTimer.cxx"

// Headers from decoder_ichikawa
#include "AmQStrTdcData_RM.h"

namespace bpo = boost::program_options;

struct StrOnlineNode : fair::mq::Device
{
    struct OptionKey {
        static constexpr std::string_view InputChannelName {"in-chan-name"};
        static constexpr std::string_view SamplingMode     {"sampling-mode"};
        static constexpr std::string_view Prescale         {"prescale"};
    };

    enum {
        kAmQData,
        kAmQHeartbeat,
        kAmQHeartbeat2nd,
        kAmqStreamingRM
    };

    enum {
        kHRTDC = 5,
        kLRTDC = 6
    };

    StrOnlineNode() {}

    void InitTask() override
    {
        using opt = OptionKey;
        fInputChannelName = fConfig->GetValue<std::string>(opt::InputChannelName.data());
        fSamplingMode     = fConfig->GetValue<std::string>(opt::SamplingMode.data());
        fPrescale         = fConfig->GetValue<int>(opt::Prescale.data());
        LOG(info) << "StrOnlineNode: InitTask Started. Channel: " << fInputChannelName << " Mode: " << fSamplingMode << " Prescale: " << fPrescale;

        gROOT->SetBatch(kTRUE);

        if (!gApplication) {
            static int argc = 1;
            static char* argv[] = {(char*)"StrOnlineNode", nullptr};
            new TApplication("StrOnlineNode", &argc, argv);
        }

        if (!fServer) {
            fServer = new THttpServer("http:8888");
            fServer->SetReadOnly(kTRUE);
            LOG(info) << "THttpServer started on port 8888";
        }

        fDrawTimer.SetDuration(100); 

        fH2HitPattern = new TH2F("h2_hitpat", "Hit Pattern;Channel;FEM ID", 128, 0, 128, 10, 0, 10);
        if(fServer) {
            fServer->Register("/Summary", fH2HitPattern);
            fServer->SetItemField("/", "_monitoring", "1000");
            fServer->SetItemField("/Summary", "_monitoring", "1000");
        }
        LOG(info) << "StrOnlineNode: InitTask Finished.";
    }

    int CheckAMANEQHeader(AmQStrTdc::Data::v1::Bits bits) {
        auto head = bits.head;
        if (head == AmQStrTdc::Data::v1::HeadTypes::Data) return kAmQData;
        if (head == AmQStrTdc::Data::v1::HeadTypes::Heartbeat) return kAmQHeartbeat;
        if (head == AmQStrTdc::Data::v1::HeadTypes::Heartbeat2nd) return kAmQHeartbeat2nd;
        if (head == AmQStrTdc::Data::v1::HeadTypes::StreamingRM) return kAmqStreamingRM;
        return -1;
    }

    bool ConditionalRun() override
    {
        try {
            // Check if the channel is registered
            bool channel_exists = false;
            for (auto const& [name, channel] : GetChannels()) {
                if (name == fInputChannelName) {
                    channel_exists = true;
                    break;
                }
            }

            if (!channel_exists) {
                std::stringstream ss;
                ss << "Channel " << fInputChannelName << " is not registered! (Available: ";
                for (auto const& [name, channel] : GetChannels()) ss << name << " ";
                ss << ")";
                LOG(error) << ss.str();
                std::this_thread::sleep_for(std::chrono::seconds(1));
                return true; 
            }

            fair::mq::Parts parts;
            if (Receive(parts, fInputChannelName, 0, 100) > 0) {
                if (fSamplingMode == "latest") {
                    while(true) {
                        fair::mq::Parts newer;
                        if (Receive(newer, fInputChannelName, 0, 0) > 0) {
                            parts = std::move(newer);
                        } else {
                            break;
                        }
                    }
                }
                
                fProcessCount++;
                if (fProcessCount % fPrescale != 0) {
                    if (fDrawTimer.Check()) UpdateDisplay();
                    return true;
                }

                // Merge all parts
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

                uint32_t fem_type = 0;
                uint32_t fem_id = 0;

                while (current_ptr < end_ptr) {
                    uint64_t buf = *current_ptr;
                    uint64_t magic = buf & 0x00FFFFFFFFFFFFFF; 

                    if (magic == TimeFrame::v1::MAGIC || magic == TimeFrame::v0::MAGIC) {
                        current_ptr += 3;
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
                        uint32_t hbf_length_words = h_hbf->length / 8;
                        if (hbf_length_words < 2) {
                            current_ptr++;
                            continue;
                        }
                        uint32_t data_length = hbf_length_words - 2; 
                        
                        current_ptr += 2; 
                        for (uint32_t i = 0; i < data_length; ++i) {
                            if (current_ptr >= end_ptr) break;
                            uint64_t word = *current_ptr;
                            AmQStrTdc::Data::v1::Bits h_amq;
                            std::memcpy(&h_amq, &word, 8);
                            
                            auto amq_flag = CheckAMANEQHeader(h_amq);
                            if (amq_flag == kAmQData) {
                                uint64_t tdc = 0, tot = 0, ch = 0;
                                bool is_hr = (fem_type == kHRTDC);
                                
                                if (fem_type == kLRTDC) {
                                    tdc = h_amq.tdc;
                                    tot = h_amq.tot;
                                    ch  = h_amq.ch;
                                } else if (fem_type == kHRTDC) {
                                    tdc = h_amq.hrtdc;
                                    tot = h_amq.hrtot;
                                    ch  = h_amq.hrch;
                                } else {
                                    current_ptr++;
                                    continue;
                                }

                                CheckAndCreateHistograms(fem_id, is_hr);
                                if (fH2HitPattern) fH2HitPattern->Fill(ch, fem_id);
                                if (fMapHitPattern.count(fem_id)) fMapHitPattern[fem_id]->Fill(ch);
                                if (fMapTDC.count(fem_id) && ch < fMapTDC[fem_id].size()) {
                                    fMapTDC[fem_id][ch]->Fill(tdc);
                                    fMapTOT[fem_id][ch]->Fill(tot);
                                }
                            }
                            current_ptr++;
                        }
                        continue;
                    }

                    if (magic == FileSinkHeader::v1::MAGIC) {
                        current_ptr += 38; 
                        continue;
                    }
                    if (magic == FileSinkTrailer::v1::MAGIC) {
                        current_ptr += 38;
                        continue;
                    }
                    if (magic == Filter::v1::MAGIC) {
                         current_ptr += 7;
                         continue;
                    }
                    current_ptr++;
                }
            }
            
            if (fDrawTimer.Check()) UpdateDisplay();

        } catch (const std::exception& e) {
            LOG(warn) << "Exception in ConditionalRun: " << e.what();
            LOG(warn) << "Entering Keep-Alive mode to maintain HTTP Server.";
            while (GetCurrentState() == fair::mq::State::Running) {
                if (fDrawTimer.Check()) UpdateDisplay();
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
        }
        return true;
    }

    std::string ToIPAddress(uint32_t fem_id) {
        uint8_t* ip = reinterpret_cast<uint8_t*>(&fem_id);
        return std::to_string(ip[3]) + "." + std::to_string(ip[2]) + "." + std::to_string(ip[1]) + "." + std::to_string(ip[0]);
    }

    void CheckAndCreateHistograms(uint32_t fem_id, bool is_hr) {
        if (fMapTDC.count(fem_id) > 0) return;

        std::string ip_str = ToIPAddress(fem_id);
        int max_ch = is_hr ? 64 : 128;
        LOG(info) << "Registering New FEM ID: " << fem_id << " (" << ip_str << ") - " << (is_hr ? "HR" : "LR");

        fMapTDC[fem_id].resize(max_ch, nullptr);
        fMapTOT[fem_id].resize(max_ch, nullptr);
        
        fMapHitPattern[fem_id] = new TH1F(Form("h1_hitpat_%s", ip_str.c_str()), 
                                          Form("FEM %s Hit Pattern;Channel;Counts", ip_str.c_str()), 
                                          max_ch, 0, max_ch);

        for (int i = 0; i < max_ch; ++i) {
            double tdc_max = is_hr ? 1000000000.0 : 500000.0;
            double tot_max = 1000.0;

            fMapTDC[fem_id][i] = new TH1F(Form("h1_tdc_%s_ch%d", ip_str.c_str(), i), 
                                          Form("FEM %s TDC Ch%d;TDC;Counts", ip_str.c_str(), i), 
                                          1000, 0, tdc_max);
            fMapTOT[fem_id][i] = new TH1F(Form("h1_tot_%s_ch%d", ip_str.c_str(), i), 
                                          Form("FEM %s TOT Ch%d;TOT;Counts", ip_str.c_str(), i), 
                                          1000, 0, tot_max);
            
            if (fServer) {
                fServer->Register(Form("/FEM_%s/TDC", ip_str.c_str()), fMapTDC[fem_id][i]);
                fServer->Register(Form("/FEM_%s/TOT", ip_str.c_str()), fMapTOT[fem_id][i]);
            }
        }
        
        if (fServer) {
            fServer->Register(Form("/FEM_%s", ip_str.c_str()), fMapHitPattern[fem_id]);
            fServer->SetItemField(Form("/FEM_%s", ip_str.c_str()), "_monitoring", "1000");
        }

        // --- Create Canvases (Matching OnlineAnalysisTFBuilder.cxx layout) ---
        TCanvas* c_tdc_1 = new TCanvas(Form("c_%s_tdc_0_31", ip_str.c_str()), Form("FEM %s TDC 0-31", ip_str.c_str()), 1600, 1000);
        c_tdc_1->Divide(8, 4);
        TCanvas* c_tdc_2 = new TCanvas(Form("c_%s_tdc_32_63", ip_str.c_str()), Form("FEM %s TDC 32-63", ip_str.c_str()), 1600, 1000);
        c_tdc_2->Divide(8, 4);
        TCanvas* c_tdc_3 = nullptr; TCanvas* c_tdc_4 = nullptr;
        if (max_ch > 64) {
            c_tdc_3 = new TCanvas(Form("c_%s_tdc_64_95", ip_str.c_str()), Form("FEM %s TDC 64-95", ip_str.c_str()), 1600, 1000);
            c_tdc_3->Divide(8, 4);
            c_tdc_4 = new TCanvas(Form("c_%s_tdc_96_127", ip_str.c_str()), Form("FEM %s TDC 96-127", ip_str.c_str()), 1600, 1000);
            c_tdc_4->Divide(8, 4);
        }

        TCanvas* c_tot_1 = new TCanvas(Form("c_%s_tot_0_31", ip_str.c_str()), Form("FEM %s TOT 0-31", ip_str.c_str()), 1600, 1000);
        c_tot_1->Divide(8, 4);
        TCanvas* c_tot_2 = new TCanvas(Form("c_%s_tot_32_63", ip_str.c_str()), Form("FEM %s TOT 32-63", ip_str.c_str()), 1600, 1000);
        c_tot_2->Divide(8, 4);
        TCanvas* c_tot_3 = nullptr; TCanvas* c_tot_4 = nullptr;
        if (max_ch > 64) {
            c_tot_3 = new TCanvas(Form("c_%s_tot_64_95", ip_str.c_str()), Form("FEM %s TOT 64-95", ip_str.c_str()), 1600, 1000);
            c_tot_3->Divide(8, 4);
            c_tot_4 = new TCanvas(Form("c_%s_tot_96_127", ip_str.c_str()), Form("FEM %s TOT 96-127", ip_str.c_str()), 1600, 1000);
            c_tot_4->Divide(8, 4);
        }

        // Set logy for TOT canvases
        for (int i = 1; i <= 32; ++i) {
            if (auto* p = c_tot_1->GetPad(i)) p->SetLogy(1);
            if (auto* p = c_tot_2->GetPad(i)) p->SetLogy(1);
            if (c_tot_3) { if (auto* p = c_tot_3->GetPad(i)) p->SetLogy(1); }
            if (c_tot_4) { if (auto* p = c_tot_4->GetPad(i)) p->SetLogy(1); }
        }

        // Register Canvases
        if (fServer) {
            fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tdc_1);
            fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tdc_2);
            if(c_tdc_3) fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tdc_3);
            if(c_tdc_4) fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tdc_4);
            fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tot_1);
            fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tot_2);
            if(c_tot_3) fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tot_3);
            if(c_tot_4) fServer->Register(Form("/FEM_%s/Canvases", ip_str.c_str()), c_tot_4);
        }

        fCanvases.push_back(c_tdc_1); fCanvases.push_back(c_tdc_2);
        if(c_tdc_3) fCanvases.push_back(c_tdc_3); if(c_tdc_4) fCanvases.push_back(c_tdc_4);
        fCanvases.push_back(c_tot_1); fCanvases.push_back(c_tot_2);
        if(c_tot_3) fCanvases.push_back(c_tot_3); if(c_tot_4) fCanvases.push_back(c_tot_4);
        
        // Initial Draw
        for (int i = 0; i < 32; ++i) {
            c_tdc_1->cd(i+1); fMapTDC[fem_id][i]->Draw();
            c_tot_1->cd(i+1); fMapTOT[fem_id][i]->Draw();
            c_tdc_2->cd(i+1); fMapTDC[fem_id][i+32]->Draw();
            c_tot_2->cd(i+1); fMapTOT[fem_id][i+32]->Draw();
            if (max_ch > 64) {
                c_tdc_3->cd(i+1); fMapTDC[fem_id][i+64]->Draw();
                c_tot_3->cd(i+1); fMapTOT[fem_id][i+64]->Draw();
                c_tdc_4->cd(i+1); fMapTDC[fem_id][i+96]->Draw();
                c_tot_4->cd(i+1); fMapTOT[fem_id][i+96]->Draw();
            }
        }
    }

    void UpdateDisplay()
    {
        for (auto* c : fCanvases) {
            if (c) {
                for(int i=0; i<32; ++i) { 
                   if(auto* pad = c->GetPad(i+1)) pad->Modified(); 
                }
                c->Update();
            }
        }
        if(gSystem) gSystem->ProcessEvents(); 
    }

private:
   std::string fInputChannelName;
   std::string fSamplingMode; // "all" or "latest"
   int fPrescale = 1;
   int fProcessCount = 0;
   KTimer fDrawTimer;

   TH2F* fH2HitPattern = nullptr;
   
   std::map<uint32_t, TH1F*> fMapHitPattern;
   std::map<uint32_t, std::vector<TH1F*>> fMapTDC;
   std::map<uint32_t, std::vector<TH1F*>> fMapTOT;

   std::vector<TCanvas*> fCanvases;
   THttpServer* fServer = nullptr;
};

void addCustomOptions(bpo::options_description& options)
{
    using opt = StrOnlineNode::OptionKey;
    options.add_options()
        (opt::InputChannelName.data(), bpo::value<std::string>()->default_value("in"), "Name of the input channel")
        (opt::SamplingMode.data(), bpo::value<std::string>()->default_value("latest"), "Sampling mode: 'all' or 'latest'")
        (opt::Prescale.data(), bpo::value<int>()->default_value(1), "Process every N-th message");
}

std::unique_ptr<fair::mq::Device> getDevice(fair::mq::ProgOptions& /*config*/)
{
    return std::make_unique<StrOnlineNode>();
}
