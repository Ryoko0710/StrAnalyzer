cmake_minimum_required(VERSION 3.16)
project(StrAnalyzer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------------------------
# Find Packages
# ------------------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH "/gpfs/home/had/kino/repos/spadi")

find_package(FairMQ REQUIRED)
find_package(ROOT REQUIRED COMPONENTS RIO RHTTP Hist)
find_package(Boost REQUIRED COMPONENTS iostreams thread program_options filesystem date_time regex)

# Dependencies from FairMQ (This follows the pattern in NestdaqAnalyzer)
if(FairMQ_FOUND)
  foreach(dep IN LISTS FairMQ_PACKAGE_DEPENDENCIES)
    if(FairMQ_${dep}_COMPONENTS)
      if(${dep} MATCHES "Boost")
        find_package(${dep} ${FairMQ_${dep}_VERSION} COMPONENTS ${FairMQ_${dep}_COMPONENTS} iostreams thread filesystem date_time regex)
      else()
        find_package(${dep} ${FairMQ_${dep}_VERSION} COMPONENTS ${FairMQ_${dep}_COMPONENTS})
      endif()
    else()
      find_package(${dep} ${FairMQ_${dep}_VERSION})
    endif()
  endforeach()
endif()

if(NOT TARGET FairLogger::FairLogger)
    find_package(FairLogger PATHS "/gpfs/home/had/kino/repos/spadi")
endif()

# ------------------------------------------------------------------------------
# Include Directories
# ------------------------------------------------------------------------------
set(NESTDAQ_USER_IMPL_DIR "/gpfs/home/had/kino/repos/spadi/src/nestdaq-user-impl")
set(DECODER_ICHIKAWA_DIR "/gpfs/home/had/kino/repos/decoder_ichikawa")

include_directories(
    ${NESTDAQ_USER_IMPL_DIR}
    ${DECODER_ICHIKAWA_DIR}
    ${FairMQ_INCDIR}
    ${ROOT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# ------------------------------------------------------------------------------
# Executables
# ------------------------------------------------------------------------------
add_executable(StrOnlineAnalyzer src/StrOnlineAnalyzer.cxx)
target_link_libraries(StrOnlineAnalyzer
    FairMQ::FairMQ
    FairLogger::FairLogger
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
install(TARGETS StrOnlineAnalyzer RUNTIME DESTINATION bin)
