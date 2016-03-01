/**
** Copyright (c) 2007-2010 Illumina, Inc.
**
** This software is covered by the "Illumina Genome Analyzer Software
** License Agreement" and the "Illumina Source Code License Agreement",
** and certain third party copyright/licenses, and any user of this
** source file is bound by the terms therein (see accompanying files
** Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
** Illumina_Source_Code_License_Agreement.pdf and third party
** copyright/license notices).
**
** This file is part of the Consensus Assessment of Sequence And VAriation
** (CASAVA) software package.
**
** @file FastqReader.h
**
** @brief This class is responsible for parsing fastq files.
**
** @author Michael Stromberg
**/

#pragma once

#include <algorithm>
#include <boost/regex.hpp>
#include <cstdlib>
#include <iostream>
#include <stdint.h>
#include <string>
#include "common/CasavaRead.hh"
#include "common/LineReader.hh"

namespace casava {
namespace common {

class FastqReader : public LineReader {
public:
    // constructor
    FastqReader(void);
    // destructor
    ~FastqReader(void);
    // returns true if there is another read available
    bool GetNextRead(CasavaRead& cr,
                     const bool isProvideHeader = true,
                     const bool isProvideQualites = true);
    // set to true if bases should be parsed, false otherwise
    void ProvideBases(bool b);

private:
    // extracts the metadata from the FASTQ header
    static void ExtractHeaderData(CasavaRead& cr, const std::string& s, bool& useCasavaHeaderStyle);
    // converts the FASTQ BQ offset (33) to the Illumina BQ offset (64)
    static inline char FastqToIlluminaOffset(char c);
    // regex that captures information from a FASTQ header file
    static const boost::regex mCasavaHeaderRegex;
    static const boost::regex mExternalHeaderRegex;
    // toggles base parsing
    bool mProvideBases;
    // uses the read name convention used by CASAVA
    bool mHasCasavaHeaderStyle;
    // used for temp line extraction:
    std::string mLine;
};

// converts the FASTQ BQ offset (33) to the Illumina BQ offset (64)
inline char FastqReader::FastqToIlluminaOffset(char c) {
    return c + 31;
}

}
}
