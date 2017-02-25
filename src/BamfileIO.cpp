//
// Created by Kevin Gori on 24/02/2017.
//

#include <unordered_set>
#include <iomanip>
#include <future>
#include <iostream>
#include "BamfileIO.h"


ClosingBamReader::ClosingBamReader(const boost::filesystem::path filename) {
    if (!this->Open(filename.string())) {
        std::cerr << "Couldn't open " << filename << " for reading" << std::endl;
    }
}

ClosingBamReader::~ClosingBamReader() {
    if (this->IsOpen()) this->Close();
}


ClosingBamWriter::ClosingBamWriter(const boost::filesystem::path filename, const BamTools::SamHeader &header, const BamTools::RefVector &refs) {
        if (!this->Open(filename.string(), header, refs)) {
            std::cerr << "Couldn't open " << filename << " for writing" << std::endl;
        }
        this->filename = filename.string();
    }

ClosingBamWriter::~ClosingBamWriter() {
    if (this->IsOpen()) {
        this->Close();
        ClosingBamReader reader(this->filename);
        reader.CreateIndex();
    }
}

const std::string & ClosingBamWriter::GetFilename() {
    return this->filename;
}


ClosingBamMultiReader::ClosingBamMultiReader(const std::vector<boost::filesystem::path> filenames) {
    std::vector<std::string> sfilenames;
    for (auto filename : filenames) {
        sfilenames.push_back(filename.string());
    }
    if (!this->Open(sfilenames)) {
        std::cerr << "Couldn't open files for reading" << std::endl;
    }
}

ClosingBamMultiReader::~ClosingBamMultiReader() {
    this->Close();
}
