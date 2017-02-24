//
// Created by Kevin Gori on 24/02/2017.
//
#include <api/BamWriter.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <utils/bamtools_pileup_engine.h>

#ifndef SIDEKICK_PILEUPUTILS_H
#define SIDEKICK_PILEUPUTILS_H

struct Region {
    Region(unsigned long r, unsigned long s, unsigned long e): ref_id(r), start(s), end(e) {}
    std::string toString();
    bool fullyLeftOf(const BamTools::BamAlignment &r);
    bool overlaps(const BamTools::BamAlignment &r);
    unsigned long ref_id;
    unsigned long start;
    unsigned long end;
};

struct CoverageVisitor : public BamTools::PileupVisitor {
public:
    CoverageVisitor(int mincoverage) : PileupVisitor(), mincoverage(mincoverage) {}
    ~CoverageVisitor() {}
    void Visit(const BamTools::PileupPosition& pileupdata);
    void Initialise(const BamTools::PileupPosition& pileupdata);
    void Finalise();

    std::vector<Region> regions;
    int mincoverage;

private:
    bool initialised = false;
    unsigned long refid = 0;
    unsigned long start = 0;
    unsigned long end = 0;
    unsigned long nreads = 0;
};

bool overlaps(const BamTools::BamAlignment &read, const Region &region, bool use_mate_info);

#endif //SIDEKICK_PILEUPUTILS_H
