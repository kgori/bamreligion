//
// Created by Kevin Gori on 24/02/2017.
//

#include <unordered_set>
#include <iomanip>
#include <future>
#include <boost/filesystem.hpp>
#include <iostream>
#include "PileupUtils.h"

std::string Region::toString() {
    std::stringstream s;
    s << "Region(" << ref_id << ", "
      << start << ", "
      << end << ")";
    return s.str();
}

bool Region::fullyLeftOf(const BamTools::BamAlignment &r) {
    return (ref_id < r.RefID || (ref_id == r.RefID && end < r.Position));
}

void CoverageVisitor::Visit(const BamTools::PileupPosition& pileupdata) {
    if (pileupdata.PileupAlignments.size() >= mincoverage) {
        if (!initialised) {
            Initialise(pileupdata);
        }
        end = pileupdata.Position;
    }
    else {
        if (initialised) {
            Finalise();
        }
        return;
    }
}

void CoverageVisitor::Initialise(const BamTools::PileupPosition& pileupdata) {
    initialised = true;
    refid = pileupdata.RefId;
    start = end = pileupdata.Position;
}

void CoverageVisitor::Finalise() {
    if (initialised) {
        regions.emplace_back(refid, start, end);
        initialised = false;
        refid = start = end = 0;
        // std::cout << regions.back().toString() << std::endl;
    }
}

bool Region::overlaps(const BamTools::BamAlignment &read) {
    bool use_mate_info = false;
    if (read.IsMapped()) use_mate_info = false;
    else if (read.IsMateMapped()) use_mate_info = true;
    else return false;
    auto read_start = use_mate_info ? read.MatePosition : read.Position;
    auto read_end = use_mate_info ? read_start + read.Length : read.GetEndPosition();
    auto read_ref_id = use_mate_info ? read.MateRefID : read.RefID;
    return (read_ref_id == this->ref_id &&
            read_start <= this->end &&
            read_end > this->start);
}