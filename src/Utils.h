//
// Created by Kevin Gori on 25/02/2017.
//
#include <boost/filesystem.hpp>

#ifndef _UTILS_H
#define _UTILS_H

int filter_bam(boost::filesystem::path query, boost::filesystem::path subject, boost::filesystem::path tmpdir,
                boost::filesystem::path outfile, int at_a_time=1000000);

#endif //_UTILS_H
