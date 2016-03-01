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
 ** @file StringUtilities.cpp
 **
 ** @brief Contains helper string functions that are used in kagu.
 **
 ** @author Michael Stromberg
 **/

#include <boost/exception/all.hpp>
#include <boost/format.hpp>
#include "common/Exceptions.hh"
#include "common/StringUtilities.hh"

using namespace std;

namespace casava {
namespace common {

// returns true if we were able to extract the value from the key value pair
bool StringUtilities::ExtractKeyValuePair(const std::string& s, const uint32_t offset, const std::string& key, std::string& val) {

    // find the start of the value
    const char* pBuffer = s.c_str();
    const char* pEnd    = pBuffer + s.size();
    const char* pKey    = key.c_str();

    // start at the offset
    pBuffer += offset;
    if(pBuffer >= pEnd) return false;

    // look for our key
    const char* pStart = strstr(pBuffer, pKey);
    if(!pStart) return false;
    pStart += key.size();

    // find the first occurrence of whitespace
    const char* pStop = pStart + 1;
    while(pStop != pEnd) {
        if(isspace(*pStop)) break;
        ++pStop;
    }

    // copy our value
    StringUtilities::CopyString(val, pStart, pStop);

    return true;
}

// returns the number of lines in the specified filename
//
// N.B.: a number of different ways were attempted, but this one was the fastest:
// All tests were performed on a FASTA file containing human genome NCBI36.
// using C fgetc:                            173.673 seconds
// using C++ getline:                         59.599 seconds
// using C read with 16 MB buffer:             8.813 seconds
// using C++ ifstream::read with 16 MB buffer: 5.723 seconds
// using C fread with 16 MB buffer:            4.214 seconds
// using gzread with 1 MB buffer:              ?.??? seconds
uint32_t StringUtilities::GetNumLines(const string& filename) {

    gzFile in = gzopen(filename.c_str(), "rb");

    if(!in) {
        BOOST_THROW_EXCEPTION(IoException(errno, (boost::format("Unable to open the file (%s) to count the number of lines") % filename).str()));
    }

    string buffer;
    buffer.resize(WC_BUFFER_SIZE + 1);
    char* pBuffer = (char*)buffer.data();

    uint32_t numLines = 0;
    int bytes_read;

    while((bytes_read = gzread(in, pBuffer, WC_BUFFER_SIZE)) > 0) {

        if(bytes_read == -1) {
            BOOST_THROW_EXCEPTION(IoException(EINVAL, (boost::format("Unable to read data from %s") % filename).str()));
        }

        char *p = pBuffer;

        while((p = (char*)memchr(p, '\n', (pBuffer + bytes_read) - p))) {
            ++p;
            ++numLines;
        }
    }

    gzclose(in);

    return numLines;
}

// returns the read name given the supplied CASAVA read data structure
string StringUtilities::GetReadName(const CasavaRead& cr) {
    ostringstream sb;
    sb << cr.Machine << ':'
        << cr.RunNumber << ':'
        << cr.FlowcellID << ':'
        << cr.Lane << ':'
        << cr.Tile << ':'
        << cr.XCoord << ':'
        << cr.YCoord;
    return sb.str();
}

// retrieves the splice length from the supplied read name (RNA)
int32_t StringUtilities::GetSpliceLength(const string& readName) {

    const uint32_t BUFFER_SIZE = (uint32_t)readName.size();
    const char* pBuffer = readName.data();

    // initialize pointers
    char* pUnderscore1 = NULL;
    char* pUnderscore2 = NULL;

    bool foundError = false;

    pUnderscore1 = (char*)memchr(pBuffer, '_', BUFFER_SIZE);
    if(!pUnderscore1) foundError = true;

    if(!foundError) pUnderscore2 = (char*)memchr(pUnderscore1 + 1, '_', BUFFER_SIZE - (pUnderscore1 - pBuffer));
    if(!pUnderscore2) foundError = true;

    if(foundError) {
        BOOST_THROW_EXCEPTION(CasavaException(EINVAL, (boost::format("Could not extract the splice length from the following read name: [%s]") % readName).str()));
    }

    string spliceLenString;
    CopyString(spliceLenString, pUnderscore1 + 1, pUnderscore2);

    return atoi(spliceLenString.c_str());
}

// splits the supplied delimited string into a vector
//
// boost split  - 8000000 iterations = 33.422 seconds
// mosaik split - 8000000 iterations = 15.593 seconds
// split        - 8000000 iterations =  7.397 seconds
void StringUtilities::Split(const string& s, const char delimiter, vector<string>& v) {

    const uint32_t sLen = (uint32_t)s.size();
    const char* pStart   = s.data();
    const char* pOld     = pStart;
    const char* pEnd     = pStart + sLen - 1;
    const char* pCurrent = NULL;

    // count the columns
    uint32_t numColumns = 0;
    while((pCurrent = (char*)memchr(pOld, delimiter, (pStart + sLen) - pOld))) {
        ++numColumns;
        if(pCurrent == pEnd) break;
        pOld = pCurrent + 1;
    }

    uint32_t numRemaining = (uint32_t)(pEnd - pOld + 1);
    if(numRemaining > 0) ++numColumns;

    // resize the vector
    v.resize(numColumns);
    vector<string>::iterator sIter = v.begin();

    // assign the vector elements
    pOld = pStart;

    while((pCurrent = (char*)memchr(pOld, delimiter, (pStart + sLen) - pOld))) {
        CopyString(*sIter, pOld, pCurrent);
        ++sIter;
        if(pCurrent == pEnd) break;
        pOld = pCurrent + 1;
    }

    numRemaining = (uint32_t)(pEnd - pOld + 1);
    if(numRemaining > 0) CopyString(*sIter, pOld, pEnd + 1);
}

}
}
