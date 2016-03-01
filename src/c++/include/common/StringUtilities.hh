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
 ** @file StringUtilities.h
 **
 ** @brief Contains helper string functions that are used in kagu.
 **
 ** @author Michael Stromberg
 **/

#pragma once

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include <string>
#include <vector>
#include <zlib.h>
#include "common/CasavaRead.hh"

#define WC_BUFFER_SIZE (16 * 1024)

namespace casava {
namespace common {

class StringUtilities {
public:
    // copies the data between two character pointers into the supplied std::string
    static inline void CopyString(std::string& s, const char* pStart, const char* pEnd);
    // returns true if we were able to extract the value from the key value pair
    static bool ExtractKeyValuePair(const std::string& s, const uint32_t offset, const std::string& key, std::string& val);
    // returns the number of lines in the specified filename
    static uint32_t GetNumLines(const std::string& filename);
    // returns the read name given the supplied CASAVA read data structure
    static std::string GetReadName(const casava::common::CasavaRead& cr);
    // retrieves the splice length from the supplied read name (RNA)
    static int32_t GetSpliceLength(const std::string& readName);
    // splits the supplied delimited std::string into a vector
    static void Split(const std::string& s, char delimiter, std::vector<std::string>& v);
};

// copies the data between two character pointers into the supplied std::string
inline void StringUtilities::CopyString(std::string& s, const char* pStart, const char* pEnd) {
    const uint32_t len = (uint32_t)(pEnd - pStart);
    s.resize(len);
    memcpy((void*)s.data(), pStart, len);
}

}
}
