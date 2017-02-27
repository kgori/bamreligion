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

struct FilePathException : public std::runtime_error {
    FilePathException(const std::string &what) : std::runtime_error(what) {}
};

struct NoResultsException : public std::runtime_error {
    NoResultsException(const std::string &what) : std::runtime_error(what) {}
};

struct FilePaths {
    FilePaths(const std::string &inputfile_,
              const std::string &working_dir_,
              const std::string &halfmapped_,
              const std::string &halfunmapped_,
              const std::string &allunmapped_,
              const std::string &filtered_,
              bool cleanup
              )
            : inputfile(fs::system_complete(inputfile_)),
              halfmapped(fs::system_complete(halfmapped_)),
              halfunmapped(fs::system_complete(halfunmapped_)),
              bothunmapped(fs::system_complete(allunmapped_)),
              working_dir(fs::path(working_dir_)),
              filtered(fs::system_complete(filtered_)),
              cleanup(cleanup)
    {
        // Create a temporary directory path
        if (working_dir.string().empty()) {
            working_dir = fs::temp_directory_path();
        }
        working_dir /= fs::unique_path();

        this->created = fs::create_directories(working_dir); // TODO: only create if all checks are OK
        if (created) std::cout << "Created the tmp path" << std::endl;
        else std::cout << "tmp path already existed" << std::endl;

        tmp_mapped = working_dir / fs::path("tmp_mapped.bam");
        tmp_mapped_filtered = working_dir / fs::path("tmp_mapped_filtered.bam");
        tmp_unmapped = working_dir / fs::path("tmp_unmapped.bam");
        tmp_both_1 = working_dir / fs::path("tmp_both_1.bam");
        tmp_both_2 = working_dir / fs::path("tmp_both_2.bam");
        tmp_both_1_filtered = working_dir / fs::path("tmp_both_1_filtered.bam");
        tmp_both_2_filtered = working_dir / fs::path("tmp_both_2_filtered.bam");

        // Checks
        if (!fs::exists(inputfile)) {
            throw FilePathException(std::string("File ") + inputfile.string() + " not found");
        }
        if (!fs::exists(halfmapped.parent_path())) {
            throw FilePathException(std::string("Path ") + halfmapped.parent_path().string()
                                     + " not found");
        }
        if (!fs::exists(halfunmapped.parent_path())) {
            throw FilePathException(std::string("Path ") + halfunmapped.parent_path().string()
                                     + " not found");
        }
        if (!fs::exists(bothunmapped.parent_path())) {
            throw FilePathException(std::string("Path ") + bothunmapped.parent_path().string()
                                     + " not found");
        }
        if (!filtered.empty() && !fs::exists(filtered.parent_path())) {
            throw FilePathException(std::string("Path ") + filtered.parent_path().string()
                                     + " not found");
        };
        if (!fs::exists(working_dir)) {
            throw FilePathException(std::string("Path ") + working_dir.string() + " not found");
        };
        std::vector<fs::path> outfiles{halfmapped, halfunmapped, bothunmapped, filtered};
        for (const auto &outfile : outfiles) {
            if (inputfile == outfile) {
                throw FilePathException(std::string("Input file " + inputfile.string()
                                                         + std::string(" will be overwritten by output file ")
                                                         + outfile.string() + " causing loss of data."));
            }
        }
    }
    ~FilePaths() {
        if (cleanup && created && fs::exists(working_dir)) {
            std::cout << "Cleaning up " << working_dir.string() << std::endl;
            fs::remove_all(working_dir);
        }
    }

    // Input file
    fs::path inputfile;

    // Output files
    fs::path halfmapped;
    fs::path halfunmapped;
    fs::path bothunmapped;
    fs::path filtered;

    // Temp files
    fs::path working_dir;
    fs::path tmp_mapped;
    fs::path tmp_mapped_filtered;
    fs::path tmp_unmapped;
    fs::path tmp_both_1;
    fs::path tmp_both_2;
    fs::path tmp_both_1_filtered;
    fs::path tmp_both_2_filtered;

private:
    bool cleanup = true;
    bool created = true;
};

std::string time_now() {
    auto t = time(nullptr);
    auto tm = *localtime(&t);
    std::stringstream ss;
    ss << std::put_time(&tm, "%d-%m-%Y %H:%M:%S");
    return ss.str();
}

template <typename Type>
void log_warning(Type &variable, const std::string &variable_name, const Type minval) {
    if (variable < minval) {
        variable = minval;
        auto t = time(nullptr);
        auto tm = *localtime(&t);
        std::cerr << "[" << time_now() << "] "
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

unsigned long initial_extraction(ClosingBamReader &reader, const FilePaths &paths, int base_qual,
                        int map_qual, bool write_filtered = false, int update_freq = 1000000)
{
    update_freq = update_freq < 1000 ? 1000 : update_freq;

    auto header = reader.GetConstSamHeader();
    auto references = reader.GetReferenceData();

    ClosingBamWriter tmp_mapped_writer(paths.tmp_mapped, header, references);
    ClosingBamWriter tmp_unmapped_writer(paths.tmp_unmapped, header, references);
    ClosingBamWriter tmp_both_1_writer(paths.tmp_both_1, header, references);
    ClosingBamWriter tmp_both_2_writer(paths.tmp_both_2, header, references);
    ClosingBamWriter filtered_writer;
    if (write_filtered) {
        ClosingBamWriter filtered_writer(paths.filtered, header, references);
    }


    std::cout << "[initial_extraction] [" << time_now()
              << "] Scanning " << reader.GetFilename()
              << " for unmapped reads" << std::endl;
    unsigned long nreads = 0;
    unsigned long nfiltered = 0;

    BamAlignment read;
    while (reader.GetNextAlignmentCore(read)) {
        nreads++;
        if (nreads % update_freq == 0) {
            std::cout << "Read " << nreads << " reads\r";
            std::flush(std::cout);
        }
        if (nreads % (update_freq*50) == 0) {
            std::cout << "\n[initial_extraction] [" << time_now() << "]" << std::endl;
        }

        if (passes_initial_checks(read)) {
            read.BuildCharData();
            if (write_filtered) {
                filtered_writer.SaveAlignment(read);
                nfiltered++;
            }

            if (passes_quality_checks(read, base_qual, map_qual)) { // At least one of pair is unmapped
                if (!read.IsMapped()) { // Current read is unmapped
                    if (!read.IsMateMapped()) { // Mate is unmapped
                        if (read.IsFirstMate()) { // Both unmapped, first in pair
                            tmp_both_1_writer.SaveAlignment(read);
                        }
                        else { // Both unmapped, second in pair
                            tmp_both_2_writer.SaveAlignment(read);
                        }
                    } else { // Current read unmapped, mate is mapped
                        tmp_unmapped_writer.SaveAlignment(read);
                    }
                } else { // Current read is mapped, mate is unmapped
                    assert (!read.IsMateMapped());
                    tmp_mapped_writer.SaveAlignment(read);
                }
            }
        }
    }
    return nfiltered;
}

int main(int argc, char** argv) {

    // Options
    int MAPQUAL = 30;
    int BASEQUAL = 10;
    int MINCOV = 1;
    bool delete_wdir = false;
    std::string _working_dir_;
    std::string _input_file_;     // input bam file
    std::string _mapped_file_;    // destination for mapped reads with unmapped mates, passing quality and coverage checks
    std::string _unmapped_file_;  // destination for unmapped reads with mapped mates, passing quality and coverage checks
    std::string _all_file_;       // destination for reads with both mates unmapped, passing quality checks
    std::string _filtered_file_;  // destination for all reads passing 'passes_initial_checks'

    po::options_description desc("Allowed options");
    desc.add_options()
    ("mapqual,q", po::value<int>(&MAPQUAL)->default_value(MAPQUAL), "Minimum mapping quality")
    ("basequal,b", po::value<int>(&BASEQUAL)->default_value(BASEQUAL), "Minimum mapping quality")
    ("coverage,c", po::value<int>(&MINCOV)->default_value(MINCOV), "Minimum mapped read coverage")
    ("working-dir,w", po::value<std::string>(&_working_dir_), "Working dir")
    ("delete,d", po::value<bool>(&delete_wdir)->default_value(delete_wdir), "Delete working dir")
    ("input,i", po::value<std::string>(&_input_file_)->required(), "Path to input file")
    ("mapped,m", po::value<std::string>(&_mapped_file_)->required(), "Path to output half-mapped file")
    ("unmapped,u", po::value<std::string>(&_unmapped_file_)->required(), "Path to output half-unmapped file")
    ("all,a", po::value<std::string>(&_all_file_)->required(), "Path to output both-unmapped file")
    ("filtered,f", po::value<std::string>(&_filtered_file_), "Path to output filtered file")
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

    // Initialise all filepaths
    try {
        const FilePaths filepaths(_input_file_, _working_dir_, _mapped_file_,
                                  _unmapped_file_, _all_file_, _filtered_file_,
                                  delete_wdir);
        bool write_filtered = !filepaths.filtered.string().empty();

        // Set option limits
        log_warning(BASEQUAL, "BASEQUAL", 0);
        log_warning(MINCOV, "MINCOV", 1);
        log_warning(MAPQUAL, "MAPQUAL", 0);

        // Print all option values
        std::cout << "MAPQUAL " << MAPQUAL << std::endl;
        std::cout << "BASEQUAL " << BASEQUAL << std::endl;
        std::cout << "MINCOV " << MINCOV << std::endl;
        std::cout << "Input file " << fs::system_complete(filepaths.inputfile).string() << std::endl;
        std::cout << "Mapped file " << fs::system_complete(filepaths.halfmapped).string() << std::endl;
        std::cout << "Unmapped file " << fs::system_complete(filepaths.halfunmapped).string() << std::endl;
        std::cout << "All file " << fs::system_complete(filepaths.bothunmapped).string() << std::endl;
        std::cout << "Filtered file " << fs::system_complete(filepaths.filtered).string() << std::endl;
        std::cout << "Write filtered?: " << (write_filtered ? "true" : "false") << std::endl;


        ClosingBamReader reader(filepaths.inputfile);
        const SamHeader header = reader.GetConstSamHeader();
        const RefVector references = reader.GetReferenceData();

        // 1: Extract all reads with at least 1 mate unmapped
        auto nfiltered = initial_extraction(reader, filepaths, BASEQUAL, MAPQUAL, write_filtered);

        // 2: Filter files to pair up mates
        std::vector<std::future<int>> filter_results;
        filter_results.push_back(
                std::async(filter_bam, filepaths.tmp_unmapped, filepaths.tmp_mapped,
                           filepaths.working_dir, filepaths.tmp_mapped_filtered, 1000000));
        filter_results.push_back(
                std::async(filter_bam, filepaths.tmp_both_1, filepaths.tmp_both_2,
                           filepaths.working_dir, filepaths.tmp_both_2_filtered, 1000000));
        filter_results.push_back(
                std::async(filter_bam, filepaths.tmp_both_2, filepaths.tmp_both_1,
                           filepaths.working_dir, filepaths.tmp_both_1_filtered, 1000000));

        // 3: Pileup and find reads passing minimum coverage threshold
        std::cout << "Checking coverage of filtered reads" << std::endl;
        auto visitor = std::make_unique<CoverageVisitor>(MINCOV);
        {
            PileupEngine pileup(true);
            pileup.AddVisitor(visitor.get());
            filter_results[0].get();  // this result needs to be ready now
            ClosingBamReader pileup_reader(filepaths.tmp_mapped_filtered);
            pileup_reader.CreateIndex();
            std::cout << "Piling up " << pileup_reader.GetFilename() << std::endl;
            BamAlignment read;
            while (pileup_reader.GetNextAlignmentCore(read)) {
                pileup.AddAlignment(read);
            }
            pileup.Flush();
        }
        visitor->Finalise();

        // Quit if no qualifying reads are found
        if (visitor->regions.empty()) {
            std::stringstream msg;
            msg << "[" << time_now() << "] "
                      << "Error - No qualifying reads were found.";
            throw NoResultsException(msg.str());
        }

        // Make sure files are ready
        for (auto it = std::next(filter_results.begin()); it != filter_results.end(); it++) {
            (*it).get();
        }

        // 4: Find mapped reads passing coverage checks to finalise half-mapped
        auto n_halfmapped = write_overlaps(filepaths.tmp_mapped_filtered, filepaths.halfmapped, visitor->regions);


        // 5: One last filter to finalise half-unmapped
        auto n_halfunmapped = filter_bam(filepaths.halfmapped, filepaths.tmp_unmapped,
                   filepaths.working_dir, filepaths.halfunmapped, 1000000);


        // 6: Consolidate to finalise both-unmapped
        std::vector<fs::path> both_mapped_paths{filepaths.tmp_both_1_filtered, filepaths.tmp_both_2_filtered};
        ClosingBamMultiReader consolidate_unmapped_reader(both_mapped_paths);
        ClosingBamWriter consolidate_unmapped_writer(filepaths.bothunmapped, consolidate_unmapped_reader.GetHeader(),
                                                     consolidate_unmapped_reader.GetReferenceData());
        BamAlignment read;
        unsigned long n_both_unmapped = 0;
        while (consolidate_unmapped_reader.GetNextAlignment(read)) {
            consolidate_unmapped_writer.SaveAlignment(read);
            n_both_unmapped++;
        }

        // Done: Write a message to confirm where the output was written
        std::cout << "Finished.\n" << std::string(60, '-') << std::endl;
        std::cout << "Wrote " << n_halfmapped
                  <<" half-mapped reads tied to high coverage areas to "
                  << filepaths.halfmapped << std::endl;

        std::cout << "Wrote " << n_halfunmapped
                  << " half-unmapped reads tied to high coverage areas to "
                  << filepaths.halfunmapped << std::endl;

        std::cout << "Wrote " << n_both_unmapped << " both-unmapped reads to "
                  << filepaths.bothunmapped << std::endl;

        if (nfiltered > 0) {
            std::cout << "Wrote " << nfiltered << " filtered reads to "
                      << filepaths.filtered << std::endl;
        }

//        if (delete_wdir) {
//            std::cout << "Cleaning up " << filepaths.working_dir << std::endl;
//            fs::remove_all(filepaths.working_dir);
//        }
        return 0;
    }
    catch (FilePathException &e) {
        std::cerr << "Error constructing filepaths: " << e.what() << std::endl;
        exit(1);
    }
    catch (NoResultsException &e) {
        std::cerr << "No results: " << e.what() << std::endl;
        exit(2);
    }
}
