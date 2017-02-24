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
        std::__1::cerr << "Couldn't open " << filename << " for reading" << std::__1::endl;
    }
}

ClosingBamReader::~ClosingBamReader() {
    if (this->IsOpen()) this->Close();
}


ClosingBamWriter::ClosingBamWriter(const boost::filesystem::path filename, const BamTools::SamHeader &header, const BamTools::RefVector &refs) {
        if (!this->Open(filename.string(), header, refs)) {
            std::__1::cerr << "Couldn't open " << filename << " for writing" << std::__1::endl;
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


ClosingBamMultiReader::ClosingBamMultiReader(const std::__1::vector<boost::filesystem::path> filenames) {
    std::__1::vector<std::__1::string> sfilenames;
    for (auto filename : filenames) {
        sfilenames.push_back(filename.string());
    }
    if (!this->Open(sfilenames)) {
        std::__1::cerr << "Couldn't open files for reading" << std::__1::endl;
    }
}

ClosingBamMultiReader::~ClosingBamMultiReader() {
    this->Close();
}
