//
// Created by Kevin Gori on 25/02/2017.
//

#include <string>
#include <unordered_set>
#include "BamfileIO.h"
#include "Utils.h"

using namespace BamTools;

void filter_bam(boost::filesystem::path query, boost::filesystem::path subject, boost::filesystem::path tmpdir, boost::filesystem::path outfile, int at_a_time) {
    std::vector<boost::filesystem::path> tmpfiles;

    ClosingBamReader query_reader(query);
    std::unordered_set<std::string> cache;

    int batch = 1;
    int written = 0;
    BamAlignment query_read;
    while (query_reader.GetNextAlignment(query_read)) {
        int nreads = 1; // already read one read
        cache.insert(query_read.Name);

        while (nreads < at_a_time && query_reader.GetNextAlignment(query_read)) {
            nreads++;
            cache.insert(query_read.Name);
        }

        std::cout << "[filter_bam] - batch number " << batch
                  << " processing " << cache.size() << " reads"
                  << std::endl;

        // Open subject file
        ClosingBamReader subject_reader(subject);

        // Open writer for this iteration
        boost::filesystem::path tmpfilename = tmpdir / boost::filesystem::unique_path();
        tmpfiles.push_back(tmpfilename);
        ClosingBamWriter batch_writer(tmpfilename, subject_reader.GetConstSamHeader(),
                                      subject_reader.GetReferenceData());

        BamAlignment subject_read;
        while(subject_reader.GetNextAlignment(subject_read)) {
            auto search = cache.find(subject_read.Name);
            if (search != cache.end()) {
                batch_writer.SaveAlignment(subject_read);
                cache.erase(search);
            }
        }

        // Cleanup
        cache.clear();
        batch++;
    }

    std::cout << "[filter_bam] - combining tmp bam files" << std::endl;
    {
        ClosingBamMultiReader multireader(tmpfiles);
        ClosingBamWriter writer(outfile, multireader.GetHeader(), multireader.GetReferenceData());

        BamAlignment multiread;
        while (multireader.GetNextAlignment(multiread)) {
            written++;
            writer.SaveAlignment(multiread);
        }
    }
    for (auto &path : tmpfiles) {
        remove(path);
        boost::filesystem::remove(path.string() + ".bai");
    }
    std::cout << "[filter_bam] - wrote " << written << " filtered reads" << std::endl;
}