/*
PROJECT: IMPALA (Inexact Matching Program ALlowing Ambiguity)
MODULE: OligoSourceQseq.h
AUTHOR: A. J. Cox

 * Copyright (c) 2003-2006 Solexa
 *
 ** This software is covered by the "Illumina Genome Analyzer Software
 ** License Agreement" and the "Illumina Source Code License Agreement",
 ** and certain third party copyright/licenses, and any user of this
 ** source file is bound by the terms therein (see accompanying files
 ** Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
 ** Illumina_Source_Code_License_Agreement.pdf and third party
 ** copyright/license notices).
 */

#ifndef CASAVA_ALIGNMENT_OLIGO_SOURCE_QSEQ_HH
#define CASAVA_ALIGNMENT_OLIGO_SOURCE_QSEQ_HH

#include "GlobalUtilities.hh"

namespace casava
{
namespace alignment
{

/*****************************************************************************/
// Read oligos from a list of QSEQ files
class OligoSourceQseq : public OligoSource
{
  public:
    OligoSourceQseq( const list<fs::path> qseqFileList,
                     const string useBases);
    ~OligoSourceQseq() { delete pSource_; }

    // Returns reference to next Sequence (supersedes getNextOligo).
    // isValid will be false if there are no sequences left.
    virtual const casava::common::Sequence& getNextSequenceSelect(bool& isValid,
                                                                  const bool isProvideHeaders,
                                                                  const bool isProvideQualities);

    // Returns reference to last Sequence fetched (supersedes getLastOligo).
    // isValid will be false if there are no sequences left.
    virtual const casava::common::Sequence& getLastSequence(bool& isValid) const;

    // Returns pointer to ASCII sequence of next oligo, or null if at end
    virtual const char* getNextOligo( void );

    // Returns pointer to ASCII sequence of last oligo fetched
    virtual const char* getLastOligo( void ) const
    {
        return pSource_->getLastOligo();
    } // ~getLastOligo( void )

    // Returns pointer to ASCII name of last oligo read
    virtual const char* getLastName( void );

    // Rewind - all streams
    virtual void rewind ( void );

  private:
    OligoSource* pSource_;
    const list<fs::path> qseqFileList_;
    list<fs::path>::const_iterator qseqFile_;
    const vector<bool> bUseBases_;
    bool format_;
    casava::common::Sequence sequence_;
    char nameBuf_[maxLineLength];

    // member variables for the second tier
    int curSeq_;
    int skippedSequences_;

    void transform( casava::common::Sequence& sequence );
};

} //namespace alignment
} //namespace casava

#endif //CASAVA_ALIGNMENT_OLIGO_SOURCE_QSEQ_HH
