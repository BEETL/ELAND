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
 ** @file LineReader.hh
 **
 ** @brief Superclass handling all of the low level buffering for
 **        parsing compressed and uncompressed text files.
 **
 ** @author Michael Stromberg
 **/

#pragma once

#include <boost/exception/all.hpp>
#include <boost/format.hpp>
#include <boost/utility.hpp>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <zlib.h>
#include "common/CasavaRead.hh"
#include "common/Exceptions.hh"

#define SR_BUFFER_SIZE 1048576

namespace casava {
namespace common {

class LineReader : boost::noncopyable
{
public:
    // constructor
    LineReader(void);
    // closes the underlying file stream(s)
    void Close(void);
    // returns true if the underlying file stream is open
    bool IsOpen(void) const;
    // opens the underlying file stream
    void Open(const std::string& filename, uint32_t numTrimPrefixBases = 0, uint32_t numTrimSuffixBases = 0);
    // rewinds the underlying file stream
    void Rewind(void);

protected:
    // prevent destruction in this base
    ~LineReader(void);
    // extracts another line from our memory buffer
    bool GetNextLine(std::string& s);
    // toggled according to the status of the underlying file stream(s)
    bool mIsOpen;
    // toggles base quality trimming
    bool mPerformTrimming;
    uint32_t mNumTrimPrefixBases;
    uint32_t mNumTrimSuffixBases;

private:
    // our underlying input stream
    gzFile mInStream;
    std::string mFilename;
    // these variables manage our getline buffer
    std::string mBuffer;
    char* mStartBuffer;
    char* mCurrentBuffer;
    int mBytesRead;
};

}
}
