//
// Created by Kevin Gori on 24/02/2017.
//
#include <boost/filesystem.hpp>
#include <api/BamReader.h>
#include <api/BamMultiReader.h>
#include <api/BamWriter.h>


#ifndef _BAMFILEIO_H
#define _BAMFILEIO_H

class ClosingBamReader : public BamTools::BamReader {
public:
    ClosingBamReader(const boost::filesystem::path filename);
    ~ClosingBamReader();
};

class ClosingBamWriter : public BamTools::BamWriter {
public:
    ClosingBamWriter(const boost::filesystem::path filename,
                     const BamTools::SamHeader &header,
                     const BamTools::RefVector &refs);
    ~ClosingBamWriter();
    const std::string & GetFilename();
private:
    std::string filename;
};

class ClosingBamMultiReader : public BamTools::BamMultiReader {
public:
    ClosingBamMultiReader(const std::vector<boost::filesystem::path> filenames);
    ~ClosingBamMultiReader();
};

#endif //_BAMFILEIO_H
