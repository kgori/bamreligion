include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
set(SOURCE_FILES
        test.cpp
        ../../src/PileupUtils.cpp
        ../../src/BamfileIO.cpp
        ../../src/Utils.cpp)

add_executable(runTests ${SOURCE_FILES})
link_directories(../../deps/bamtools/lib)
target_link_libraries(runTests gtest gtest_main)
target_link_libraries(runTests bamtools)