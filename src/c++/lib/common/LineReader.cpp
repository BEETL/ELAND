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
 ** @file LineReader.cpp
 **
 ** @brief Superclass handling all of the low level buffering for
 **        parsing compressed and uncompressed text files.
 **
 ** @author Michael Stromberg
 **/

#include "common/LineReader.hh"
#include "common/StringUtilities.hh"

using namespace std;

namespace casava {
namespace common {

// constructor
LineReader::LineReader(void)
    : mIsOpen(false)
    , mPerformTrimming(false)
    , mNumTrimPrefixBases(0)
    , mNumTrimSuffixBases(0)
{
    mBuffer.resize(SR_BUFFER_SIZE);
    mStartBuffer = (char*)mBuffer.data();
    mBytesRead = 0;
}

// destructor
LineReader::~LineReader(void) {
    Close();
}

// closes the underlying file stream(s)
void LineReader::Close(void) {
    if(mIsOpen) {
        gzclose(mInStream);
        mIsOpen = false;
    }
}

// extracts another line from our memory buffer
bool LineReader::GetNextLine(string& s) {

    // skip if the file is not currently open or if we don't have any data in the buffer
    if(!mIsOpen || (mBytesRead <= 0)) return false;

    char* p = NULL;
    if((p = (char*)memchr(mCurrentBuffer, '\n', (mStartBuffer + mBytesRead) - mCurrentBuffer))) {

        StringUtilities::CopyString(s, mCurrentBuffer, p);
        mCurrentBuffer = p + 1;

    } else {

        const int32_t remainingLen = (int32_t)(mStartBuffer + SR_BUFFER_SIZE - mCurrentBuffer);
        memmove(mStartBuffer, mCurrentBuffer, remainingLen);

        mBytesRead = gzread(mInStream, mStartBuffer + remainingLen, SR_BUFFER_SIZE - remainingLen);

        if(mBytesRead == -1) {
            BOOST_THROW_EXCEPTION(IoException(EINVAL, (boost::format("Unable to read data from %s") % mFilename).str()));
        }

        if(mBytesRead == 0) return false;

        mBytesRead += remainingLen;
        mCurrentBuffer = mStartBuffer;

        // if we can find any newlines after fetching new data, give up
        if((p = (char*)memchr(mCurrentBuffer, '\n', (mStartBuffer + mBytesRead) - mCurrentBuffer))) {
            StringUtilities::CopyString(s, mCurrentBuffer, p);
            mCurrentBuffer = p + 1;
        } else return false;
    }

    return true;
}

// returns true if the underlying file stream is open
bool LineReader::IsOpen(void) const {
    return mIsOpen;
}

// opens the underlying file stream
void LineReader::Open(const string& filename, uint32_t numTrimPrefixBases, uint32_t numTrimSuffixBases) {

    mFilename = filename;
    mInStream = gzopen(filename.c_str(), "rb");

    if(!mInStream) {
        BOOST_THROW_EXCEPTION(IoException(errno, (boost::format("Unable to open the file (%s) for reading") % mFilename).str()));
    }

    mIsOpen = true;

    // fill the buffer
    mBytesRead = gzread(mInStream, mStartBuffer, SR_BUFFER_SIZE);

    if(mBytesRead == -1) {
        BOOST_THROW_EXCEPTION(IoException(EINVAL, (boost::format("Unable to read data from %s") % mFilename).str()));
    }

    mCurrentBuffer = mStartBuffer;

    // localize the trimming data
    mNumTrimPrefixBases = numTrimPrefixBases;
    mNumTrimSuffixBases = numTrimSuffixBases;

    if((numTrimPrefixBases > 0) || (numTrimSuffixBases > 0)) mPerformTrimming = true;
}

// rewinds the underlying file stream
void LineReader::Rewind(void) {
    gzrewind(mInStream);
    mBytesRead = gzread(mInStream, mStartBuffer, SR_BUFFER_SIZE);

    if(mBytesRead == -1) {
        BOOST_THROW_EXCEPTION(IoException(EINVAL, (boost::format("Unable to read data from %s") % mFilename).str()));
    }

    mCurrentBuffer = mStartBuffer;
}

}
}
