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
 ** (CASAVA) software package. **
 ** \file OligoSourceFastq.cpp
 **
 ** \brief Reads fastq files.
 **
 ** \author Roman Petrovski
 **/

#include <boost/foreach.hpp>

#include "alignment/OligoSourceFastq.hh"
#include "blt_util/parse_util.hh"
#include "eland_ms/ElandConstants.hh"

namespace casava
{
namespace alignment
{

const list<fs::path> OligoSourceFastq::makeInputFileList(const fs::path &inputDirectory,
                                                         const std::string &sample,
                                                         const std::string &barcode,
                                                         const unsigned int lane, const unsigned int read,
                                                         const std::vector<unsigned int> &clusterSets)
{
    list<fs::path> ret;
    BOOST_FOREACH(unsigned int cluster, clusterSets)
    {
        ret.push_back(inputDirectory / (boost::format(
                "%s_%s_L%03u_R%d_%03u.fastq.gz") % sample % barcode % lane % read % cluster).str());
    }
    return ret;
}
/*****************************************************************************/
// ctor
OligoSourceFastq::OligoSourceFastq(const fs::path &inputDirectory,
                                   const std::string &sample,
                                   const std::string &barcode,
                                   const unsigned int lane, const unsigned int read,
                                   const std::vector<unsigned int> &clusterSets,
                                   const std::string &useBases)
    : fastqFiles_(makeInputFileList(inputDirectory, sample, barcode, lane, read, clusterSets))
    , fastqFilesIterator_(fastqFiles_.begin())
    , curSeq_(1)
    , isNameBuf_(false)
{
    cerr << "Using fastq files as input..." << endl;

    // initialize use bases data:
    std::vector<bool> bUseBases(expandUseBases(useBases));

    useBasesLength_=bUseBases.size();

    bool state(false);
    unsigned length(0);
    for(unsigned i(0);i<useBasesLength_;++i){
        if(bUseBases[i] == state) {
            length++;
        } else {
           if(length) ubRegions_.push_back(std::make_pair(length,state));
           length=1;
           state=(! state);
        }
    }
    if(length) ubRegions_.push_back(std::make_pair(length,state));
}

//std::ostream & operator << (std::ostream &os, const cc::CasavaRead &read)
//{
//    return os << read.Machine << " " << read.RunNumber << " " << read.Lane << " "
//              << read.Tile << " " << read.XCoord << " " << read.YCoord << " " << read.Index << " " << read.ReadNumber;
//}
/*****************************************************************************/
// Returns reference to next Sequence (supersedes getNextOligo).
// isValid will be false if there are no sequences left.
const casava::common::Sequence& OligoSourceFastq::getNextSequenceSelect(bool& isValid,
                                                                        const bool isProvideHeader,
                                                                        const bool isProvideQualities)
{
    using casava::blt_util::parse_int_str;
    using casava::blt_util::parse_unsigned_str;

    isValid = false;
    while (fastqFiles_.end() != fastqFilesIterator_)
    {
        if (!reader_.IsOpen())
        {
            reader_.Open(fastqFilesIterator_->string());
        }
        while (reader_.GetNextRead(read_,isProvideHeader,isProvideQualities))
        {
            curSeq_++;
            if( isNoMask_ || mask_[curSeq_-1] )
            {
                if(isProvideHeader) {
                    const int ret = snprintf(nameBuf_,nameBufSize_,">%s_%04u:%u:%u:%d:%d#%s/%u",
                             read_.Machine.c_str(),
                             parse_unsigned_str(read_.RunNumber),
                             parse_unsigned_str(read_.Lane),
                             parse_unsigned_str(read_.Tile),
                             parse_int_str(read_.XCoord),
                             parse_int_str(read_.YCoord),
                             read_.Index.c_str(),
                             parse_unsigned_str(read_.ReadNumber));

                    if((ret<0) || (ret>=nameBufSize_)) {
                             BOOST_THROW_EXCEPTION(
                                 cc::CasavaException(
                                     EINVAL,(boost::format("ERROR: Failed to format header string from fastq record number %i in file: '%s'\n")
                                         % (curSeq_-1) % (fastqFilesIterator_->string())).str()));
                    }
                    isNameBuf_=true;
                } else {
                    isNameBuf_=false;
                }
                transform(read_, sequence_,isProvideQualities);
                isValid = true;
                return sequence_;
            }
        }
        reader_.Close();
        ++fastqFilesIterator_;
    }
    return sequence_;
}

// Returns reference to last Sequence fetched (supersedes getLastOligo).
// isValid will be false if no sequences have been successfully read.
const casava::common::Sequence& OligoSourceFastq::getLastSequence(bool& isValid) const
{
    isValid = reader_.IsOpen();
    return sequence_;
}

// Returns pointer to ASCII sequence of next oligo, or null if at end
const char* OligoSourceFastq::getNextOligo(void)
{
    bool isValid(false);
    const casava::common::Sequence& sequence(getNextSequence(isValid));
    return (isValid ? sequence.getData().c_str() : NULL);
}

const char * OligoSourceFastq::getLastOligo() const
{
    bool isValid = true;
    const casava::common::Sequence& sequence(getLastSequence(isValid));
    return (isValid ? sequence.getData().c_str() : 0);
}

// Returns pointer to ASCII name of last oligo read
const char* OligoSourceFastq::getLastName(void)
{
    if(reader_.IsOpen() and isNameBuf_) return nameBuf_;
    return NULL;
}

// Rewind - next oligo read will be first in list
void OligoSourceFastq::rewind(void)
{
    reader_.Close();
    fastqFilesIterator_ = fastqFiles_.begin();
    curSeq_=1;
}


/*****************************************************************************/
// Applies UseBases mask
void OligoSourceFastq::transform( const cc::CasavaRead &read,
                                  cc::Sequence& sequence,
                                  const bool isProvideQualities)
{
    if (read.Bases.length() != useBasesLength_)
    {
        BOOST_THROW_EXCEPTION(
                cc::InvalidParameterException(
                        (boost::format("Tried to apply a %d-cycle UseBases mask "
                              "to data actually containing %d bases.")
                            % useBasesLength_ % read.Bases.length()).str()));
    }

    UBRegions_t::const_iterator i(ubRegions_.begin());
    const UBRegions_t::const_iterator i_end(ubRegions_.end());

    sequence.getData().clear();
    const char* iData(read.Bases.c_str());

    for(;i!=i_end;++i) {
       if(i->second) {
         sequence.getData().append(iData,i->first);
       }
       iData += i->first;
    }

    if(! isProvideQualities) return;

    i=ubRegions_.begin();

    sequence.getQuality().clear();
    const char* iQuality(read.Qualities.c_str());

    for(;i!=i_end;++i) {
       if(i->second) sequence.getQuality().append(iQuality,i->first);
       iQuality += i->first;
    }
}

} //namespace alignment
} //namespace casava
