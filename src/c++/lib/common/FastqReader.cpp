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
 ** @file FastqReader.cpp
 **
 ** @brief This class is responsible for parsing fastq files.
 **
 ** @author Michael Stromberg
 **/

#include "common/FastqReader.hh"

using namespace std;

namespace casava {
namespace common {

// regular expressions
const boost::regex FastqReader::mCasavaHeaderRegex("^@([^:]*):([^:]*):([^:]*):([^:]+):([^:]+):([^:]+):(\\S+)\\s+([^:]+):([^:]+):([^:]+):(\\S*)");
const boost::regex FastqReader::mExternalHeaderRegex("^@(\\S+)");

// constructor
FastqReader::FastqReader(void)
    : mProvideBases(true)
    , mHasCasavaHeaderStyle(true)
{}

// destructor
FastqReader::~FastqReader(void) {}

// extracts the metadata from the FASTQ header
void FastqReader::ExtractHeaderData(CasavaRead& cr, const string& s, bool& useCasavaHeaderStyle) {

    boost::smatch headerResults;

    // expected line: @EAS139:136:FC706VJ:2:5:996:13539 1:Y:22:ATCACG
    if(useCasavaHeaderStyle && boost::regex_search(s, headerResults, mCasavaHeaderRegex)) {

        cr.Machine    = headerResults[1].str();
        cr.RunNumber  = headerResults[2].str();
        cr.FlowcellID = headerResults[3].str();
        cr.Lane       = headerResults[4].str();
        cr.Tile       = headerResults[5].str();
        cr.XCoord     = headerResults[6].str();
        cr.YCoord     = headerResults[7].str();

        cr.ReadNumber    = headerResults[8].str();
        cr.FailedFilters = (headerResults[9].str()[0] == 'Y' ? true : false);
        cr.ControlID     = headerResults[10].str();
        cr.Index         = headerResults[11].str();

        if(cr.Index.empty()) cr.Index = "0";

    } else {

        useCasavaHeaderStyle = false;
        if(!boost::regex_search(s, headerResults, mExternalHeaderRegex)) {
            BOOST_THROW_EXCEPTION(CasavaException(EINVAL, (boost::format("Regular expression (mExternalHeaderRegex) could not applied to the following line: [%s]") % s).str()));
        }

        cr.Machine = headerResults[1].str();
        cr.RunNumber.clear();
        cr.FlowcellID.clear();
        cr.Lane.clear();
        cr.Tile.clear();
        cr.XCoord.clear();
        cr.YCoord.clear();

        cr.ReadNumber.clear();
        cr.FailedFilters = false;
        cr.ControlID.clear();
        cr.Index = "0";
    }
}

// returns true if there is another read available
bool FastqReader::GetNextRead(CasavaRead& cr,
                              const bool isProvideHeader,
                              const bool isProvideQualities) {

    // return false if our file stream is closed
    if(!mIsOpen) return false;

    // retrieve the header
    if(!GetNextLine(mLine)) return false;

    // sanity check
    if(mLine.at(0) != '@') {
        BOOST_THROW_EXCEPTION(CasavaException(EINVAL, (boost::format("Expected a '@' in the FASTQ header, found '%c'.") % mLine.at(0)).str()));
    }

    if(isProvideHeader) {
        // extract the header data
        ExtractHeaderData(cr, mLine, mHasCasavaHeaderStyle);
    }

    // retrieve the bases
    if(!GetNextLine(mLine)) {
        BOOST_THROW_EXCEPTION(IoException(EINVAL, (boost::format("Premature end of the FASTQ file while retrieving bases.")).str()));
    }

    const uint32_t numBases = (uint32_t)mLine.size();
    if(mProvideBases) cr.Bases = mLine;

    // skip the second header
    if(!GetNextLine(mLine)) {
        BOOST_THROW_EXCEPTION(IoException(EINVAL, (boost::format("Premature end of the FASTQ file while retrieving second header.")).str()));
    }

    // retrieve the qualities
    if(!GetNextLine(mLine)) {
        BOOST_THROW_EXCEPTION(IoException(EINVAL, (boost::format("Premature end of the FASTQ file while retrieving base qualities.")).str()));
    }

    if(! isProvideQualities) return true;

    const uint32_t numBaseQualities = (uint32_t)mLine.size();

    // sanity check
    if(mProvideBases && (numBases != numBaseQualities)) {
        BOOST_THROW_EXCEPTION(CasavaException(EINVAL, (boost::format("The number of base qualities (%u) do not match the number of bases (%u).") % numBaseQualities % numBases).str()));
    }

    // modify the base qualities
    cr.Qualities = (mPerformTrimming ? mLine.substr(mNumTrimPrefixBases, numBaseQualities - mNumTrimPrefixBases - mNumTrimSuffixBases) : mLine);
    transform(cr.Qualities.begin(), cr.Qualities.end(), cr.Qualities.begin(), FastqToIlluminaOffset);

    return true;
}

// set to true if bases should be parsed, false otherwise
void FastqReader::ProvideBases(bool b) {
    mProvideBases = b;
}

}
}
