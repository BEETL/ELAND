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
 ** \file OligoSourceFastq.hh
 **
 ** \brief Reads BCL files and other associated files (filter, position).
 **
 ** \author Roman Petrovski
 **/


#ifndef CASAVA_ALIGNMENT_OLIGO_SOURCE_FASTQ_HH
#define CASAVA_ALIGNMENT_OLIGO_SOURCE_FASTQ_HH

#include "GlobalUtilities.hh"
#include "common/FastqReader.hh"

namespace casava
{
namespace alignment
{

namespace cc=casava::common;

/*****************************************************************************/
// Read oligos from a list of QSEQ files
class OligoSourceFastq : public OligoSource
{
public:
    OligoSourceFastq(const fs::path &inputDirectory,
                    const std::string &sample,
                    const std::string &barcode,
                    const unsigned int lane, const unsigned int read,
                    const std::vector<unsigned int> &clusterSets,
                    const std::string &useBases);

    // Returns reference to next Sequence (supersedes getNextOligo).
    // isValid will be false if there are no sequences left.
    virtual const casava::common::Sequence& getNextSequenceSelect(bool& isValid,
                                                                  const bool isProvideHeader,
                                                                  const bool isProvideQualities);

    // Returns reference to last Sequence fetched (supersedes getLastOligo).
    // isValid will be false if there are no sequences left.
    virtual const casava::common::Sequence& getLastSequence(bool& isValid) const;

    // Returns pointer to ASCII sequence of next oligo, or null if at end
    virtual const char* getNextOligo( void );

    // Returns pointer to ASCII sequence of last oligo fetched
    const char * getLastOligo() const;
    // Returns pointer to ASCII name of last oligo read
    virtual const char* getLastName( void );

    // Rewind - all streams
    virtual void rewind ( void );

  private:
    const list<fs::path> fastqFiles_;
    list<fs::path>::const_iterator fastqFilesIterator_;
    unsigned useBasesLength_;
    typedef std::vector<std::pair<unsigned,bool> > UBRegions_t;
    UBRegions_t ubRegions_;
    cc::Sequence sequence_;
    enum { nameBufSize_ = 4096 };
    char nameBuf_[nameBufSize_];

    cc::CasavaRead read_;
    cc::FastqReader reader_;

    // member variables for the second tier
    int curSeq_;
    int skippedSequences_;
    bool isNameBuf_;

    void transform( const cc::CasavaRead &read, cc::Sequence& sequence, const bool isProvideQualities);

    static const list<fs::path> makeInputFileList(const fs::path &inputDirectory,
                                                  const std::string &sample,
                                                  const std::string &barcode,
                                                  const unsigned int lane, const unsigned int read,
                                                  const std::vector<unsigned int> &clusterSets);

};

} //namespace alignment
} //namespace casava

#endif //CASAVA_ALIGNMENT_OLIGO_SOURCE_QSEQ_HH
