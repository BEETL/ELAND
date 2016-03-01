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
 ** @file AlignmentReader.cpp
 **
 ** @brief This class handles the simultaneous parsing of ELAND extended files
 **        and compressed FASTQ files.
 **
 ** @author Michael Stromberg
 **/

#include "kagu/AlignmentReader.h"

namespace casava {
namespace kagu {

using namespace std;
namespace cc = casava::common;

// constructor
AlignmentReader::AlignmentReader(void)
    : mIsOpen(false)
    , mProvideBQs(false)
    , mUntrimmedReadLength(0)
{
    // we extract our bases from the ELAND extended file
    mBQReader.ProvideBases(false);
}

// destructor
AlignmentReader::~AlignmentReader(void) {}

// closes the file streams
void AlignmentReader::Close() {
    if(mIsOpen) {
        mIsOpen = false;
        mAlignmentReader.Close();
        mBQReader.Close();
    }
}

// gets the read length from the ELAND extended file
uint32_t AlignmentReader::GetReadLength(void) {
    return mUntrimmedReadLength;
}

// opens the associated alignment and base quality files
void AlignmentReader::Open(const string& alignmentFilename, const Filenames_t& bqFilenames,
    uint32_t numTrimPrefixBases, uint32_t numTrimSuffixBases, ReferenceRenamingStrategy_t strategy) {

        // open the alignment file
        mAlignmentReader.Open(alignmentFilename, 0, 0);
        mAlignmentReader.SetReferenceRenamingStrategy(strategy);

        // open the base qualities file
        mBaseQualityFilenames = bqFilenames;
        mBqFilenameIter = mBaseQualityFilenames.begin();
        mBQReader.Open(*mBqFilenameIter, numTrimPrefixBases, numTrimSuffixBases);

        // toggle the open state
        mIsOpen = true;

        // localize the trimming data
        mNumTrimPrefixBases = numTrimPrefixBases;
        mNumTrimSuffixBases = numTrimSuffixBases;

        // grab the read length
        cc::CasavaRead cr;
        mAlignmentReader.GetNextRead(cr);
        mUntrimmedReadLength = (uint32_t)cr.Bases.size();
        mAlignmentReader.Rewind();
}

// returns the next read
bool AlignmentReader::GetNextRead(cc::CasavaRead& cr) {

    // sanity check: make sure the file is open
    if(!mIsOpen) return false;

    // add the information from the alignment file
    if(!mAlignmentReader.GetNextRead(cr)) return false;

    // add the information from the fastq file
    if(mProvideBQs) {

        bool moreReadsAvailable = mBQReader.GetNextRead(cr);

        // try the next fastq file
        if(!moreReadsAvailable) {
            ++mBqFilenameIter;

            // no more fastq files
            if(mBqFilenameIter == mBaseQualityFilenames.end()) {
                BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("More entries are available in the ELAND extended file, but all of the entries in the fastq files have already been processed. Are we missing some fastq files?")).str()));
            }

            mBQReader.Close();
            mBQReader.Open(*mBqFilenameIter, mNumTrimPrefixBases, mNumTrimSuffixBases);

            moreReadsAvailable = mBQReader.GetNextRead(cr);
        }

        // not a good sign
        if(!moreReadsAvailable) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("The alignment reader was able to retrieve the next entry from the ELAND extended file, but not from the fastq file.")).str()));
        }

        // sanity check
        if(cr.Qualities.size() != cr.Bases.size()) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("The number of bases (%u) in the ELAND extended file is not equal to the number of bases (%u) in the fastq file. Please check your use bases parameters (--ub1 and --ub2).") % cr.Qualities.size() % cr.Bases.size()).str()));
        }
    }

    // reset our qualities
    cr.MateAlignmentQuality     = 0;
    cr.FragmentAlignmentQuality = 0;

    return true;
}

// returns true if the underlying file stream is open
bool AlignmentReader::IsOpen(void) const {
    return mIsOpen;
}

// set to true if base qualities should be parsed, false otherwise
void AlignmentReader::ProvideBaseQualities(bool b) {
    mProvideBQs = b;
}

// rewinds the input files to the beginning
void AlignmentReader::Rewind(void) {
    mAlignmentReader.Rewind();
    mBqFilenameIter = mBaseQualityFilenames.begin();
    mBQReader.Close();
    mBQReader.Open(*mBqFilenameIter, mNumTrimPrefixBases, mNumTrimSuffixBases);
}

}
}
