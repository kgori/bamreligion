#include <iostream>
#include <boost/program_options.hpp>
#include "BamfileIO.h"
#include "PileupUtils.h"
#include "Utils.h"
#include <future>
#include <iomanip>

using namespace BamTools;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

template <typename Type>
void log_warning(Type &variable, const std::string &variable_name, const Type minval) {
    if (variable < minval) {
        variable = minval;
        auto t = time(nullptr);
        auto tm = *localtime(&t);
        std::cerr << "[" << std::put_time(&tm, "%d-%m-%Y %H:%M:%S") << "] "
                  << "Warning - value of " << variable_name << " set to minimum of "
                  << minval << std::endl;
    }
}

unsigned long write_overlaps(const fs::path &infile, const fs::path &outfile,
                             const std::vector<Region> &regions) {
    if (regions.empty()) {
        return 0;
    }
    unsigned long nreads = 0;
    ClosingBamReader reader(infile);
    size_t i = 0;
    auto region = regions[i];

    ClosingBamWriter writer(outfile, reader.GetConstSamHeader(), reader.GetReferenceData());
//    std::cout << "[write_overlaps] writing to " << writer.GetFilename() << std::endl;
    BamAlignment read;
    while (reader.GetNextAlignmentCore(read)) {
        if (region.overlaps(read)) {
            writer.SaveAlignment(read);
            nreads++;
        }
        else if (region.fullyLeftOf(read)) {
            i++;
            if (i > regions.size()) break;
            region = regions[i];
        }
    }
    return nreads;
}

bool passes_initial_checks(const BamAlignment &r) {
    return (!r.IsDuplicate() &&
            r.IsPaired() &&
            !r.IsProperPair() &&
            !r.IsFailedQC() &&
            (r.AlignmentFlag & 0x0800) == 0 &&
            r.IsPrimaryAlignment() &&
            (!r.IsMapped() || !r.IsMateMapped()));
}

double avg_base_quality(const BamAlignment &r) {
    double totalqual = 0;
    for (auto basequal : r.Qualities) {
        totalqual += basequal - 33;
    }
    return totalqual / r.Qualities.length();
}

bool passes_quality_checks(const BamAlignment &r, int base_qual, int map_qual) {
    return r.IsMapped() ? r.MapQuality >= map_qual : avg_base_quality(r) >= base_qual;
}

int main(int argc, char** argv) {

    // Options
    int MAPQUAL = 30;
    int BASEQUAL = 10;
    int MINCOV = 1;
    bool delete_wdir = false;
    std::string _working_dir_ = ".bamreligion.tmp";
    std::string _input_file_;
    std::string _mapped_file_;
    std::string _unmapped_file_;
    std::string _all_file_;

    po::options_description desc("Allowed options");
    desc.add_options()
    ("mapqual,q", po::value<int>(&MAPQUAL)->default_value(MAPQUAL), "Minimum mapping quality")
    ("basequal,b", po::value<int>(&BASEQUAL)->default_value(BASEQUAL), "Minimum mapping quality")
    ("coverage,c", po::value<int>(&MINCOV)->default_value(MINCOV), "Minimum mapped read coverage")
    ("working-dir,w", po::value<std::string>(&_working_dir_)->default_value(_working_dir_), "Working dir")
    ("delete,d", po::value<bool>(&delete_wdir)->default_value(delete_wdir), "Delete working dir")
    ("input,i", po::value<std::string>(&_input_file_)->required(), "Path to input file")
    ("mapped,m", po::value<std::string>(&_mapped_file_)->required(), "Path to output half-mapped file")
    ("unmapped,u", po::value<std::string>(&_unmapped_file_)->required(), "Path to output half-unmapped file")
    ("all,a", po::value<std::string>(&_all_file_)->required(), "Path to output both-unmapped file")
    ("help,h", "Show help")
    ;
    po::variables_map vm;

    try {
        po::store(po::parse_command_line(argc, (const char *const *) argv, desc), vm); // can throw
        /**
        --help option
        */
        if ( vm.count("help") ) {
            std::cout << "Basic Command Line Parameter App"
                      << std::endl
                      << desc
                      << std::endl;
            return 0;
        }
        po::notify(vm); // throws on error, so do after help in case
        // there are any problems
    }
    catch(po::error& e) {
        std::cerr << "ERROR: "
                  << e.what()
                  << std::endl
                  << std::endl;
        std::cerr << desc << std::endl;
        return 1;
    }

    fs::path working_dir(_working_dir_);
    fs::create_directories(working_dir);

    // Initialise all filepaths
    fs::path filename(_input_file_);
    fs::path hm_out(_mapped_file_);
    fs::path hu_out(_unmapped_file_);
    fs::path au_out(_all_file_);
    assert(fs::exists(filename));
    assert(fs::exists(fs::system_complete(hm_out).parent_path()));
    assert(fs::exists(fs::system_complete(hu_out).parent_path()));
    assert(fs::exists(fs::system_complete(au_out).parent_path()));

    // Initialise tempfiles
    fs::path tmp_mapped = working_dir / fs::path("tmp_mapped.bam");
    fs::path tmp_unmapped = working_dir / fs::path("tmp_unmapped.bam");
    fs::path tmp_both_1 = working_dir / fs::path("tmp_both_1.bam");
    fs::path tmp_both_2 = working_dir / fs::path("tmp_both_2.bam");
    fs::path tmp_mapped_filtered = working_dir / fs::path("tmp_mapped_filtered.bam");
//    fs::path tmp_unmapped_filtered = working_dir / fs::path("tmp_unmapped_filtered.bam");
    fs::path tmp_both_1_filtered = working_dir / fs::path("tmp_both_1_filtered.bam");
    fs::path tmp_both_2_filtered = working_dir / fs::path("tmp_both_2_filtered.bam");

    // Set option limits
//    if (MAPQUAL < 1) {
//        auto t = std::time(nullptr);
//        auto tm = *std::localtime(&t);
//        std::cerr << "[" << std::put_time(&tm, "%d-%m-%Y %H:%M:%S") << "] "
//                  << "Warning - value of " << MAPQUAL << "set to minimum of 1"
//                  << std::endl;
//    }
//
//    if (MAPQUAL < 0) {
//        auto t = std::time(nullptr);
//        auto tm = *std::localtime(&t);
//        std::cerr << "[" << std::put_time(&tm, "%d-%m-%Y %H:%M:%S") << "] "
//                  << "Warning - value of " << MAPQUAL << "set to minimum of 1"
//                  << std::endl;
//    }

    log_warning(BASEQUAL, "BASEQUAL", 0);
    log_warning(MINCOV, "MINCOV", 1);
    log_warning(MAPQUAL, "MAPQUAL", 0);

    // Print all option values
    std::cout << "MAPQUAL " << MAPQUAL << std::endl;
    std::cout << "BASEQUAL " << BASEQUAL << std::endl;
    std::cout << "MINCOV " << MINCOV << std::endl;
    std::cout << "Input file " << fs::system_complete(filename).string() << std::endl;
    std::cout << "Mapped file " << fs::system_complete(hm_out).string() << std::endl;
    std::cout << "Unmapped file " << fs::system_complete(hu_out).string() << std::endl;
    std::cout << "All file " << fs::system_complete(au_out).string() << std::endl;


    ClosingBamReader reader(filename);
    const SamHeader header = reader.GetConstSamHeader();
    const RefVector references = reader.GetReferenceData();

    std::vector<unsigned long> counts {0, 0, 0, 0, 0, 0};

    {
        ClosingBamWriter tmp_mapped_writer(tmp_mapped, header, references);
        ClosingBamWriter tmp_unmapped_writer(tmp_unmapped, header, references);
        ClosingBamWriter tmp_both_1_writer(tmp_both_1, header, references);
        ClosingBamWriter tmp_both_2_writer(tmp_both_2, header, references);
        ClosingBamWriter debug(working_dir / fs::path("debug.bam"), header, references);


        std::cout << "Scanning " << reader.GetFilename() << " for unmapped reads\n" << std::endl;
        int nreads = 0;

        BamAlignment read;
        while (reader.GetNextAlignmentCore(read)) {
            nreads++;
            if (nreads % 1000000 == 0) {
                std::cout << "Read " << nreads << " reads\r";
                std::flush(std::cout);
            }

            if (passes_initial_checks(read)) {
                read.BuildCharData();
                debug.SaveAlignment(read);
                counts[0]++;
                if (passes_quality_checks(read, BASEQUAL, MAPQUAL)) { // At least one of pair is unmapped
                    counts[1]++;
                    if (!read.IsMapped()) { // Current read is unmapped
                        if (!read.IsMateMapped()) { // Mate is unmapped
                            if (read.IsFirstMate()) { // Both unmapped, first in pair
                                tmp_both_1_writer.SaveAlignment(read);
                                counts[2]++;
                            }
                            else { // Both unmapped, second in pair
                                tmp_both_2_writer.SaveAlignment(read);
                                counts[3]++;
                            }
                        } else { // Current read unmapped, mate is mapped
                            tmp_unmapped_writer.SaveAlignment(read);
                            counts[4]++;
                        }
                    } else { // Current read is mapped, mate is unmapped
                        assert (!read.IsMateMapped());
                        tmp_mapped_writer.SaveAlignment(read);
                        counts[5]++;
                    }
                }
            }
        }
//        std::cout << std::endl;
//        std::cout << counts[0] << '\t'
//                  << counts[1] << '\t'
//                  << counts[2] << '\t'
//                  << counts[3] << '\t'
//                  << counts[4] << '\t'
//                  << counts[5] << std::endl;
    }

    std::vector<std::future<int>> filter_results;
    filter_results.push_back(
            std::async(filter_bam, tmp_unmapped, tmp_mapped, working_dir, tmp_mapped_filtered, 1000000));
//    filter_results.push_back(
//            std::async(filter_bam, tmp_mapped, tmp_unmapped, working_dir, tmp_unmapped_filtered, 1000000));
    filter_results.push_back(
            std::async(filter_bam, tmp_both_1, tmp_both_2, working_dir, tmp_both_2_filtered, 1000000));
    filter_results.push_back(
            std::async(filter_bam, tmp_both_2, tmp_both_1, working_dir, tmp_both_1_filtered, 1000000));

    // Pileup
    std::cout << "Checking coverage of filtered reads" << std::endl;

    auto visitor = std::make_unique<CoverageVisitor>(MINCOV);
    {
        PileupEngine pileup(true);
        pileup.AddVisitor(visitor.get());
        filter_results[0].get();  // this result needs to be ready now
        ClosingBamReader pileup_reader(tmp_mapped_filtered);
        pileup_reader.CreateIndex();
        std::cout << "Piling up " << pileup_reader.GetFilename() << std::endl;
        BamAlignment read;
        while (pileup_reader.GetNextAlignmentCore(read)) {
            pileup.AddAlignment(read);
            //std::cout << "\t" << read.RefID << "\t" << read.Position << std::endl;
        }
        pileup.Flush();
    }
    visitor->Finalise();

    if (visitor->regions.empty()) {
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::cerr << "[" << std::put_time(&tm, "%d-%m-%Y %H:%M:%S") << "] "
                  << "Error - No qualifying reads were found." << std::endl;
        exit(1);
    }

    // Make sure files are ready
    for (auto it = std::next(filter_results.begin()); it != filter_results.end(); it++) {
        (*it).get();
    }
      // Replace this with a final filter_bam call
//    std::cout << "Writing " << write_overlaps(tmp_unmapped_filtered, hu_out, visitor->regions)
//              << " half-unmapped reads tied to high coverage areas" << std::endl;

    std::cout << "Wrote " << write_overlaps(tmp_mapped_filtered, hm_out, visitor->regions)
              <<" half-mapped reads tied to high coverage areas" << std::endl;

    // One last filter
    std::cout << "Wrote " << filter_bam(hm_out, tmp_unmapped, working_dir, hu_out, 1000000)
              << " half-unmapped reads tied to high coverage areas" << std::endl;

    std::vector<fs::path> paths{tmp_both_1_filtered, tmp_both_2_filtered};
    ClosingBamMultiReader consolidate_unmapped_reader(paths);
    ClosingBamWriter consolidate_unmapped_writer(au_out, consolidate_unmapped_reader.GetHeader(),
                                                 consolidate_unmapped_reader.GetReferenceData());
    BamAlignment read;
    unsigned int n = 0;
    while (consolidate_unmapped_reader.GetNextAlignment(read)) {
        consolidate_unmapped_writer.SaveAlignment(read);
        n++;
    }
    std::cout << "Wrote " << n << " both-unmapped reads to "
              << consolidate_unmapped_writer.GetFilename() << std::endl;

    if (delete_wdir) {
        fs::remove_all(working_dir);
    }
    return 0;
}
