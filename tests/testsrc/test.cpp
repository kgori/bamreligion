//
// Created by Kevin Gori on 25/02/2017.
//

#include <BamfileIO.h>
#include "gtest/gtest.h"
#include "Utils.h"

namespace fs = boost::filesystem;

TEST(test, test_filter_bam) {
    int written = filter_bam(fs::path("../data/query.bam"),
                             fs::path("../data/subject.bam"),
                             fs::path("."),
                             fs::path("../data/result.bam"),
                             5);

    std::vector<std::string> names;
    ClosingBamReader checker(fs::path("../data/result.bam"));
    BamTools::BamAlignment read;
    while (checker.GetNextAlignment(read)) {
        names.push_back(read.Name);
    }
    fs::remove(fs::path("../data/result.bam"));

    ASSERT_EQ(written, 7);  // Expecting 7 reads written

    // This is the expected order the reads should come out
    ASSERT_STREQ(names[0].c_str(), "SOLEXA-1GA-2_2_FC20EMB:5:195:284:685");
    ASSERT_STREQ(names[1].c_str(), "SOLEXA-1GA-2_2_FC20EMB:5:35:583:827");
    ASSERT_STREQ(names[2].c_str(), "SOLEXA-1GA-2_2_FC20EMB:5:248:130:724");
    ASSERT_STREQ(names[3].c_str(), "SOLEXA-1GA-2_2_FC20EMB:5:236:644:107");
    ASSERT_STREQ(names[4].c_str(), "SOLEXA-1GA-2_2_FC20EMB:5:165:628:70");
    ASSERT_STREQ(names[5].c_str(), "SOLEXA-1GA-2_2_FC20EMB:5:108:485:455");
    ASSERT_STREQ(names[6].c_str(), "SOLEXA-1GA-2_2_FC20EMB:5:240:501:237");
}