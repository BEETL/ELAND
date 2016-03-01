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
** @file AlignmentReader.h
**
** @brief This class handles the simultaneous parsing of ELAND extended files
**        and compressed FASTQ files.
**
** @author Michael Stromberg
**/

#pragma once

#include <boost/exception/all.hpp>
#include <boost/format.hpp>
#include <cstdlib>
#include <iostream>
#include <stdint.h>
#include <string>
#include <vector>
#include "common/CasavaRead.hh"
#include "common/ElandExtendedReader.h"
#include "common/Exceptions.hh"
#include "common/FastqReader.hh"
#include "kagu/ConfigurationSettings.h"
#include "kagu/KaguDataTypes.h"

namespace casava {
namespace kagu {

class AlignmentReader {
public:
    // constructor
    AlignmentReader(void);
    // destructor
    ~AlignmentReader(void);
    // closes the file streams
    void Close();
    // gets the read length from the ELAND extended file
    uint32_t GetReadLength(void);
    // opens the associated alignment and base quality files
    void Open(const std::string& alignmentFilename, const Filenames_t& bqFilenames,
        uint32_t numTrimPrefixBases, uint32_t numTrimSuffixBases,
        ReferenceRenamingStrategy_t strategy);
    // returns the next read
    bool GetNextRead(casava::common::CasavaRead& cr);
    // returns true if the underlying file stream is open
    bool IsOpen(void) const;
    // set to true if base qualities should be parsed, false otherwise
    void ProvideBaseQualities(bool b);
    // rewinds the input files to the beginning
    void Rewind(void);

private:
    // toggled according to the status of the underlying file streams
    bool mIsOpen;
    // toggles base quality file parsing
    bool mProvideBQs;
    uint32_t mNumTrimPrefixBases;
    uint32_t mNumTrimSuffixBases;
    // our input streams
    casava::common::ElandExtendedReader mAlignmentReader;
    casava::common::FastqReader mBQReader;
    // our read length before applying any use bases modifications
    uint32_t mUntrimmedReadLength;
    // our base quality (fastq) filenames
    Filenames_t mBaseQualityFilenames;
    Filenames_t::const_iterator mBqFilenameIter;
};

}
}
