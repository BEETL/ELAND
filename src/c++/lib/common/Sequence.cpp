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
 ** \file Alignment.cpp
 **
 ** \brief Implementation of the I/O API for sequence data type (qseq files).
 **
 ** \author Come Raczy
 **/

#include <errno.h>
#include <boost/lexical_cast.hpp>

#include "common/Exceptions.hh"
#include "common/Sequence.hh"

namespace casava
{
namespace common
{

UnfilteredSequence::UnfilteredSequence(const Spot &spot, const std::string& index,
        unsigned int readNumber, const std::string &data,
        const std::string &quality) :
    spot_(spot), index_(index), readNumber_(readNumber), data_(data), quality_(
            quality)
{
}

UnfilteredSequence::UnfilteredSequence(const UnfilteredSequence &unfiltered) :
    spot_(unfiltered.spot_), index_(unfiltered.index_), readNumber_(
            unfiltered.readNumber_), data_(unfiltered.data_), quality_(
            unfiltered.quality_)
{
}

UnfilteredSequence &UnfilteredSequence::operator=(
        const UnfilteredSequence &unfiltered)
{
    if (&unfiltered != this)
    {
        spot_ = unfiltered.spot_;
        index_ = unfiltered.index_;
        readNumber_ = unfiltered.readNumber_;
        data_ = unfiltered.data_;
        quality_ = unfiltered.quality_;
    }
    return *this;
}

bool UnfilteredSequence::operator==(const UnfilteredSequence &unfiltered) const
{
    return &unfiltered == this || (spot_ == unfiltered.spot_ && index_
            == unfiltered.index_ && readNumber_ == unfiltered.readNumber_
            && data_ == unfiltered.data_ && quality_ == unfiltered.quality_);
}

bool UnfilteredSequence::operator!=(const UnfilteredSequence &unfiltered) const
{
    return !(unfiltered == *this);
}

/*****************************************************************************/

void UnfilteredSequence::mask(const std::vector<uint>& cycleIndexVec)
{
    std::string maskedSeqStr;
    std::string maskedQualStr;
    const unsigned int origNumCycles(data_.length());

    for (std::vector<uint>::const_iterator cycleCIter(cycleIndexVec.begin()); cycleCIter
            != cycleIndexVec.end(); ++cycleCIter)
    {
        const unsigned int cycleInd(*cycleCIter);

        if (cycleInd >= origNumCycles)
        {
            BOOST_THROW_EXCEPTION(CasavaException(
                EINVAL, "Specified cycle ("
                + boost::lexical_cast<std::string>(cycleInd + 1)
                + ") is out of range (1 - "
                + boost::lexical_cast<std::string>(origNumCycles) + ")."));
        }

        maskedSeqStr += data_[cycleInd];
        maskedQualStr += quality_[cycleInd];
    }

    data_ = maskedSeqStr;
    quality_ = maskedQualStr;
}

/*****************************************************************************/

std::ostream &operator<<(std::ostream &os, const UnfilteredSequence &unfiltered)
{
    os << unfiltered.spot_;
    os.put('\t');
    os.write(unfiltered.index_.c_str(),unfiltered.index_.size());
		os.put('\t');
    putUnsignedInteger(os, unfiltered.readNumber_);
    os.put('\t');
    os.write(unfiltered.data_.c_str(), unfiltered.data_.size());
    os.put('\t');
    os.write(unfiltered.quality_.c_str(), unfiltered.quality_.size());
    return os;
}

std::istream &operator>>(std::istream &is, UnfilteredSequence &unfiltered)
{
    if (is >> unfiltered.spot_)
    {
      is.ignore();
      std::getline(is, unfiltered.index_, '\t');
      getUnsignedInteger(is, unfiltered.readNumber_, true);
      std::getline(is, unfiltered.data_, '\t');
      std::stringbuf sb;
      is.get(sb, '\t'); // the quality string is never empty
      unfiltered.quality_ = sb.str();
    }
    return is;
}

Sequence::Sequence(const Spot &spot, const std::string& index,
        unsigned int readNumber, const std::string &data,
        const std::string &quality, bool passed) :
    UnfilteredSequence(spot, index, readNumber, data, quality), passed_(passed)
{
}

Sequence::Sequence(const UnfilteredSequence &unfiltered, bool passed) :
    UnfilteredSequence(unfiltered), passed_(passed)
{
}

Sequence::Sequence(const Sequence &sequence) :
    UnfilteredSequence(sequence), passed_(sequence.passed_)
{
}

Sequence &Sequence::operator=(const Sequence &sequence)
{
    if (&sequence != this)
    {
        UnfilteredSequence::operator=(sequence);
        passed_ = sequence.passed_;
    }
    return *this;
}

bool Sequence::operator==(const Sequence &sequence) const
{
    return &sequence == this || (UnfilteredSequence::operator==(sequence)
            && passed_ == sequence.passed_);
}

bool Sequence::operator!=(const Sequence &sequence) const
{
    return !(sequence == *this);
}

std::ostream &operator<<(std::ostream &os, const Sequence &sequence)
{
    os << static_cast<const UnfilteredSequence &> (sequence);
    os.put('\t');
    putBool<'1', '0'> (os, sequence.passed_);
    os.put('\n'); // the delimiter for the end of the sequence
    return os;
}

std::istream &operator>>(std::istream &is, Sequence &sequence)
{
    if (is >> static_cast<UnfilteredSequence &> (sequence))
    {
        is.ignore();
        getBool<'1', '0'> (is, sequence.passed_);
        is.ignore(); // the delimiter for the end of the sequence
    }
    return is;
}

} // namespace common
} // namespace casava
