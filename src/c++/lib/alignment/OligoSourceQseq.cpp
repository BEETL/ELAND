/**
 ** Copyright (c) 2007-2009 Illumina, Inc.
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
 ** \file OligoSourceQseq.cpp
 **
 ** \brief Commonly used functions and definitions
 **
 ** \author Tony Cox
 **/

//#include <iostream>
//#include <ctime>
//#include <dirent.h>
//#include <boost/foreach.hpp>
//#include <boost/format.hpp>

#include "alignment/OligoSourceQseq.hh"
#include "eland_ms/ElandConstants.hh"

namespace casava
{
namespace alignment
{

/*****************************************************************************/
// **
// ** Function definitions for OligoSourceQseq
// **
/*****************************************************************************/


/*****************************************************************************/
// ctor
OligoSourceQseq::OligoSourceQseq( const list<fs::path> qseqFileList,
                                  const string useBases)
        : pSource_(NULL)
        , qseqFileList_(qseqFileList)
        , bUseBases_(expandUseBases(useBases))
{
    qseqFile_ = qseqFileList_.begin();
    pSource_ = getOligoSource(qseqFile_->string().c_str());

    curSeq_ = 1;
    skippedSequences_ = 0;

    cerr << "Using qseq files as input..." << endl;
}

/*****************************************************************************/
// Returns reference to next Sequence (supersedes getNextOligo).
// isValid will be false if there are no sequences left.
const casava::common::Sequence& OligoSourceQseq::getNextSequenceSelect(bool& isValid,
                                                                       const bool,
                                                                       const bool)
{
    if (!(isValid = (qseqFile_ != qseqFileList_.end())))
    {
        return sequence_;
    }

    casava::common::Sequence& sequence(
            const_cast<casava::common::Sequence&>( pSource_->getNextSequence(isValid) )
    );
    sprintf(nameBuf_, ">%s_%04u:%u:%u:%d:%d#%s/%u",
            sequence.getSpot().getTile().getMachineName().c_str(),
            sequence.getSpot().getTile().getRunNumber(),
            sequence.getSpot().getTile().getLaneNumber(),
            sequence.getSpot().getTile().getTileNumber(),
            sequence.getSpot().getX(), sequence.getSpot().getY(),
            sequence.getIndex().c_str(), sequence.getReadNumber() );

    while (isValid)
    {
        curSeq_++;
        if( isNoMask_ || mask_[curSeq_-1])
        {
            transform(sequence);
            return sequence;
        }
        else
        {
            casava::common::Sequence& sequence_tmp(
                const_cast<casava::common::Sequence&>( pSource_->getNextSequence(isValid) )
                );
            sprintf(nameBuf_, ">%s_%04u:%u:%u:%d:%d#%s/%u",
                    sequence_tmp.getSpot().getTile().getMachineName().c_str(),
                    sequence_tmp.getSpot().getTile().getRunNumber(),
                    sequence_tmp.getSpot().getTile().getLaneNumber(),
                    sequence_tmp.getSpot().getTile().getTileNumber(),
                    sequence_tmp.getSpot().getX(), sequence_tmp.getSpot().getY(),
                    sequence_tmp.getIndex().c_str(), sequence_tmp.getReadNumber() );
            sequence = sequence_tmp;
        }
    }

    if( !isValid ) {
        // No more sequences? try next tile
        if (++qseqFile_ != qseqFileList_.end())
        {
            delete pSource_;
            pSource_ = getOligoSource(qseqFile_->string().c_str());
        }
        return getNextSequence(isValid);
    } // ~else
    // The calling code is going to test isValid to figure out that the return value is invalid
    // TODO: improve this interface!
    return sequence_;

} // OligoSourceQseq::getNextSequence()


// Returns reference to last Sequence fetched (supersedes getLastOligo).
// isValid will be false if no sequences have been successfully read.
const casava::common::Sequence& OligoSourceQseq::getLastSequence(bool& isValid) const
{
    return ((isValid = (qseqFile_ != qseqFileList_.end()))
            ? pSource_->getLastSequence(isValid)
            : sequence_);
} // ~OligoSourceQseq::getLastSequence()

// Returns pointer to ASCII sequence of next oligo, or null if at end
const char* OligoSourceQseq::getNextOligo(void)
{
    bool isValid(false);
    const casava::common::Sequence& sequence(getNextSequence(isValid));
    return (isValid ? sequence.getData().c_str() : NULL);
} // ~OligoSourceQseq::getNextOligo()

// Returns pointer to ASCII name of last oligo read
const char* OligoSourceQseq::getLastName(void)
{
    if (qseqFile_ == qseqFileList_.end())
        return NULL;
    return nameBuf_;
} // ~OligoSourceQseq::getLastName()

// Rewind - next oligo read will be first in list
void OligoSourceQseq::rewind(void)
{
    delete pSource_;
    qseqFile_ = qseqFileList_.begin();
    pSource_ = getOligoSource(qseqFile_->string().c_str());
    // reset counter
    curSeq_=1;
} // ~OligoSourceQseq::rewind()


/*****************************************************************************/
// Applies UseBases mask and converts '.' -> 'N'
void OligoSourceQseq::transform( casava::common::Sequence& sequence )
{
    string data = sequence.getData();
    string quality = sequence.getQuality();
    string d,q;
    if (data.length() != bUseBases_.size())
        cerr << "Tried to apply a " << bUseBases_.size() << "-cycle UseBases mask "
             << "to data actually containing " << data.length() << " bases." << endl;
    assert( data.length() == bUseBases_.size() );

    string::const_iterator iData = data.begin();
    string::const_iterator iQuality = quality.begin();

    for(vector<bool>::const_iterator iUse = bUseBases_.begin();
        iUse != bUseBases_.end();
        ++iUse, ++iData, ++iQuality)
    {
        unsigned char c = static_cast<unsigned char>(*iData);
        if (*iUse)
        {
            d.append( 1U, (c != '.') ? c : 'N' );
            q.append( 1U, static_cast<unsigned char>(*iQuality) );
        }
    }
    sequence.setData(d);
    sequence.setQuality(q);
    // hopefully, the length returned will be shorter, which prevents
    // a spurious second-call to this method
}

} //namespace alignment
} //namespace casava
