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
 ** @file ElandExtendedReader.cpp
 **
 ** @brief This class is responsible for parsing ELAND extended files.
 **
 ** @author Michael Stromberg
 **/

#include "common/ElandExtendedReader.h"
#include "common/StringUtilities.hh"

using namespace std;

namespace ck = casava::kagu;

namespace casava {
namespace common {

// regex that captures information from a VMF file
const boost::regex ElandExtendedReader::mPositionsRegex("^(-?\\d+)(F|R)(\\S+)$");

// constructor
ElandExtendedReader::ElandExtendedReader(void)
    : mUseContigNames(false)
    , mUseReferenceNames(false)
    , mProvideReadName(false)
{}

// destructor
ElandExtendedReader::~ElandExtendedReader(void) {}

// populates the supplied casava read with the ELAND extended read name
void ElandExtendedReader::ExtractReadName(CasavaRead& cr, const char* pBegin, const char* pEnd) {

    bool foundError = false;
    const uint32_t BUFFER_SIZE = (uint32_t)(pEnd - pBegin);

    // initialize pointers
    char* pUnderscore = NULL;
    char* pColon1     = NULL;
    char* pColon2     = NULL;
    char* pColon3     = NULL;
    char* pColon4     = NULL;
    char* pHash       = NULL;
    char* pSlash      = NULL;

    // HWI-BRUNOP16X_0001:8:1:3302:1060#0/1
    pUnderscore = (char*)memchr(pBegin, '_', BUFFER_SIZE);
    if(!pUnderscore) foundError = true;

    if(!foundError) pColon1 = (char*)memchr(pUnderscore + 1, ':', BUFFER_SIZE - (pUnderscore - pBegin));
    if(!pColon1) foundError = true;

    if(!foundError) pColon2 = (char*)memchr(pColon1 + 1, ':', BUFFER_SIZE - (pColon1 - pBegin));
    if(!pColon2) foundError = true;

    if(!foundError) pColon3 = (char*)memchr(pColon2 + 1, ':', BUFFER_SIZE - (pColon2 - pBegin));
    if(!pColon3) foundError = true;

    if(!foundError) pColon4 = (char*)memchr(pColon3 + 1, ':', BUFFER_SIZE - (pColon3 - pBegin));
    if(!pColon4) foundError = true;

    if(!foundError) pHash = (char*)memchr(pColon4 + 1, '#', BUFFER_SIZE - (pColon4 - pBegin));
    if(!pHash) foundError = true;

    if(!foundError) pSlash = (char*)memchr(pHash + 1, '/', BUFFER_SIZE - (pHash - pBegin));
    if(!pSlash) foundError = true;

    if(foundError) {
        string readName;
        StringUtilities::CopyString(readName, pBegin, pEnd);
        BOOST_THROW_EXCEPTION(CasavaException(EINVAL, (boost::format("Read name extraction failed: [%s]") % readName).str()));
    }

    StringUtilities::CopyString(cr.Machine,    pBegin,          pUnderscore);
    StringUtilities::CopyString(cr.RunNumber,  pUnderscore + 1, pColon1);
    StringUtilities::CopyString(cr.Lane,       pColon1 + 1,     pColon2);
    StringUtilities::CopyString(cr.Tile,       pColon2 + 1,     pColon3);
    StringUtilities::CopyString(cr.XCoord,     pColon3 + 1,     pColon4);
    StringUtilities::CopyString(cr.YCoord,     pColon4 + 1,     pHash);
    StringUtilities::CopyString(cr.Index,      pHash + 1,       pSlash);
    StringUtilities::CopyString(cr.ReadNumber, pSlash + 1,      pEnd);
}

// returns true if there is another read available
bool ElandExtendedReader::GetNextRead(CasavaRead& cr) {

    // return false if our file stream is closed
    if(!mIsOpen) return false;

    // get the next line from the alignment and base quality files
    string line;
    if(!GetNextLine(line)) return false;

    // extract the VMF fields
    //
    // using memchr & CopyString instead of boost::regex sped up kagu 12 %
    const uint32_t BUFFER_SIZE = (uint32_t)line.size();
    const char* pBuffer = line.data();
    const char* pEnd    = pBuffer + BUFFER_SIZE;

    bool foundError = false;

    // initialize pointers
    char* pTab1 = NULL;
    char* pTab2 = NULL;
    char* pTab3 = NULL;

    pTab1 = (char*)memchr(pBuffer, '\t', BUFFER_SIZE);
    if(!pTab1) foundError = true;

    if(!foundError) pTab2 = (char*)memchr(pTab1 + 1, '\t', BUFFER_SIZE - (pTab1 - pBuffer));
    if(!pTab2) foundError = true;

    if(!foundError) pTab3 = (char*)memchr(pTab2 + 1, '\t', BUFFER_SIZE - (pTab2 - pBuffer));
    if(!pTab3) foundError = true;

    if(foundError) {
        BOOST_THROW_EXCEPTION(CasavaException(EINVAL, (boost::format("Tab-delimited splitting could not applied to the following line: [%s]") % line).str()));
    }

    if(mProvideReadName) ExtractReadName(cr, pBuffer + 1, pTab1);
    StringUtilities::CopyString(cr.Bases,     pTab1 + 1, pTab2);
    StringUtilities::CopyString(cr.Status,    pTab2 + 1, pTab3);
    StringUtilities::CopyString(cr.Positions, pTab3 + 1, pEnd);

    cr.IsNm    = false;
    cr.IsQc    = false;
    cr.IsTmm   = false;
    cr.MStatus = MS_Unknown;

    // determine if the read is actually aligned
    const bool hasTooManyMatches = (cr.Positions[0] == '-' ? true : false);
    const string::size_type colonPos = cr.Status.find(':');
    const bool isAligned = (!hasTooManyMatches && (colonPos != string::npos) ? true : false);

    // extract the seed errors
    if(colonPos != string::npos) {

        const string::size_type secondColon = cr.Status.find(':', colonPos + 1);

        if(secondColon == string::npos) {
            BOOST_THROW_EXCEPTION(CasavaException(EINVAL, (boost::format("Unable to find the second colon in the neighborhood string: [%s]") % cr.Status).str()));
        }

        // extract the three strings
        cr.SeedErrors[0] = atoi(cr.Status.substr(0, colonPos).c_str());
        cr.SeedErrors[1] = atoi(cr.Status.substr(colonPos + 1, secondColon - colonPos - 1).c_str());
        cr.SeedErrors[2] = atoi(cr.Status.substr(secondColon + 1).c_str());

        if(hasTooManyMatches) {
            cr.IsTmm   = true;
            cr.MStatus = MS_Repeat;
        }

    } else {

        const char statusChar = cr.Status[0];

        if(statusChar == 'N') {
            cr.IsNm    = true;
            cr.MStatus = MS_NM;
        } else if(statusChar == 'Q') {
            cr.IsQc    = true;
            cr.MStatus = MS_QC;
        }

        cr.SeedErrors[0] = 0;
        cr.SeedErrors[1] = 0;
        cr.SeedErrors[2] = 0;
    }

    // extract the positions
    if(isAligned) {

        // split and assign the comma delimited positions
        vector<string> posVector;
        StringUtilities::Split(cr.Positions, ',', posVector);

        if(posVector.size() == 1) {
            cr.MStatus = MS_SingleAlignmentFound;
        } else {
            cr.MStatus = MS_ManyAlignmentsFound;
        }

        cr.Alignments.resize(posVector.size());
        CasavaAlignments::iterator caIter = cr.Alignments.begin();

        string currentReferenceName;
        string currentContigName;
        boost::smatch positionsResults;
        vector<string>::iterator sIter;
        for(sIter = posVector.begin(); sIter != posVector.end(); ++sIter, ++caIter) {

            // find the colon delimiter (if present)
            string::size_type colonPos = sIter->find(':');
            if(colonPos != string::npos) {

                // find the forward slash contig name delimiter
                currentReferenceName = sIter->substr(0, colonPos);
                currentContigName.clear();

                string::size_type slashPos = currentReferenceName.find('/');

                if(slashPos != string::npos) {
                    if(mUseContigNames) {
                        currentReferenceName = currentReferenceName.substr(slashPos + 1);
                    } else if(mUseReferenceNames) {
                        currentReferenceName = currentReferenceName.substr(0, slashPos);
                    } else {
                        currentContigName    = currentReferenceName.substr(slashPos + 1);
                        currentReferenceName = currentReferenceName.substr(0, slashPos);
                    }
                }

                *sIter = sIter->substr(colonPos + 1);
            }

            if(!boost::regex_search(*sIter, positionsResults, mPositionsRegex)) {
                BOOST_THROW_EXCEPTION(CasavaException(EINVAL, (boost::format("Regular expression (mPositionsRegex) could not applied to the following line: [%s]") % *sIter).str()));
            }

            caIter->ReferenceName     = currentReferenceName;
            caIter->ContigName        = currentContigName;
            caIter->IsReverseStrand   = (positionsResults[2].str()[0] == 'R' ? true : false);
            caIter->MatchDescriptor   = positionsResults[3].str();
            caIter->ReferencePosition = atoi(positionsResults[1].str().c_str());
        }

    } else cr.Alignments.clear();

    return true;
}

// set to true if base qualities should be parsed, false otherwise
void ElandExtendedReader::ProvideReadName(bool b) {
    mProvideReadName = b;
}

// sets the desired reference renaming strategy (contig, reference, or both)
void ElandExtendedReader::SetReferenceRenamingStrategy(ck::ReferenceRenamingStrategy_t strategy) {
    if(strategy == ck::USE_CONTIG_NAME)         mUseContigNames    = true;
    else if(strategy == ck::USE_REFERENCE_NAME) mUseReferenceNames = true;
}

}
}
