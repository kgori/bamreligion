cmake_minimum_required(VERSION 3.7)
project(bamreligion)

if(NOT "${CMAKE_SOURCE_DIR}/deps/bamtools/include/api/BamAlignment.h")
    execute_process(COMMAND ${CMAKE_SOURCE_DIR}/build_deps.sh)
endif()

set(CMAKE_CXX_STANDARD 14)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
link_libraries(Threads::Threads)

find_package(Boost 1.50.0 COMPONENTS filesystem system program_options)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
        link_directories(${Boost_LIBRARY_DIR_DEBUG})
        link_libraries(
                ${Boost_SYSTEM_LIBRARY_DEBUG}
                ${Boost_FILESYSTEM_LIBRARY_DEBUG}
                ${Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG})
    ELSE()
        link_directories(${Boost_LIBRARY_DIR_RELEASE})
        link_libraries(
                ${Boost_SYSTEM_LIBRARY_RELEASE}
                ${Boost_FILESYSTEM_LIBRARY_RELEASE}
                ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE})
    ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)

endif()

set(SOURCE_FILES
        src/main.cpp src/PileupUtils.cpp src/BamfileIO.cpp src/Utils.cpp
        deps/bamtools/src/utils/bamtools_pileup_engine.cpp)
include_directories(deps/bamtools/include deps/bamtools/src src)
link_directories(deps/bamtools/lib)
link_libraries(bamtools)
add_executable(bamreligion ${SOURCE_FILES})

add_subdirectory(tests)
