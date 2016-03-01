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
 ** \file Sequence.hh
 **
 ** \brief I/O API for sequence data type (qseq files).
 **
				 ** Stream-based API based on the Reader and Writer concepts defined
 ** in the "FilteringStreams.hh" header. The API factors out an
 ** UnfilteredSequence to ease reuse, for instance for Alignment data.
 **
 ** \author Come Raczy
 **/

#ifndef CASAVA_COMMON_SEQUENCE_HH
#define CASAVA_COMMON_SEQUENCE_HH

#include <string>
#include <vector>

#include "common/Spot.hh"
#include "common/FilteringStreams.hh"

typedef unsigned int uint;

namespace casava
{
namespace common
{

/**
 ** Definition of an unfiltered sequence.
 **
 ** This class is introduced to simplify the implementation of other
 ** file formats such as the format of the "export" files.
 **/
class UnfilteredSequence
{
public:
    UnfilteredSequence(const Spot &spot = Spot(), const std::string& index = "",
            unsigned int readNumber = 0, const std::string &data = "",
            const std::string &quality = "");
    UnfilteredSequence(const UnfilteredSequence &unfiltered);
    UnfilteredSequence &operator=(const UnfilteredSequence &unfiltered);
    bool operator==(const UnfilteredSequence &unfiltered) const;
    bool operator!=(const UnfilteredSequence &unfiltered) const;
    Spot getSpot() const
    {
        return spot_;
    }
    const std::string &getMachineName() const
    {
        return spot_.getTile().getMachineName();
    }
    int getRunNumber() const
    {
        return spot_.getTile().getRunNumber();
    }
    int getLaneNumber() const
    {
        return spot_.getTile().getLaneNumber();
    }
    int getTileNumber() const
    {
        return spot_.getTileNumber();
    }
    int getX() const
    {
        return spot_.getX();
    }
    int getY() const
    {
        return spot_.getY();
    }
    const std::string& getIndex() const
    {
        return index_;
    }
    std::string& getIndex()
    {
        return index_;
    }
    unsigned int getReadNumber() const
    {
        return readNumber_;
    }
    const std::string& getData() const
    {
        return data_;
    }
    std::string& getData()
    {
        return data_;
    }
    unsigned int getLength() const
    {
        return data_.length();
    }
    const std::string& getQuality() const
    {
        return quality_;
    }
    std::string& getQuality()
    {
        return quality_;
    }
    void setMachineName(const std::string &machineName)
    {
        spot_.getTile().setMachineName(machineName);
    }
    void setRunNumber(const unsigned int runNumber)
    {
        spot_.getTile().setRunNumber(runNumber);
    }
    void setLaneNumber(const unsigned int laneNumber)
    {
        spot_.getTile().setLaneNumber(laneNumber);
    }
    void setTileNumber(const unsigned int tileNumber)
    {
        spot_.getTile().setTileNumber(tileNumber);
    }
    void setSpot(const Spot &spot)
    {
        spot_ = spot;
    }
    void setX(const int x)
    {
        spot_.setX(x);
    }
    void setY(const int y)
    {
        spot_.setY(y);
    }
    void setIndex(const std::string& index)
    {
        index_ = index;
    }
    void setReadNumber(unsigned int readNumber)
    {
        readNumber_ = readNumber;
    }
    void setData(const std::string &data)
    {
        data_ = data;
    }
    void setQuality(const std::string &quality)
    {
        quality_ = quality;
    }

    /**
     ** @brief Mask according to supplied vector of indices.
     **/
    void mask(const std::vector<uint>& cycleIndexVec);

private:
    Spot spot_;
    std::string index_;
    unsigned int readNumber_;
    std::string data_;
    std::string quality_;
    friend std::ostream &operator<<(std::ostream &os,
            const UnfilteredSequence &unfiltered);
    friend std::istream &operator>>(std::istream &is,
            UnfilteredSequence &unfiltered);
};

/**
 ** @brief Write a complete object of type Sequence, EXcluding the delimiter.
 **/
std::ostream
        &operator<<(std::ostream &os, const UnfilteredSequence &unfiltered);
/**
 ** @brief Read a complete object of type Sequence, EXcluding the delimiter.
 **/
std::istream &operator>>(std::istream &is, UnfilteredSequence &unfiltered);

/**
 ** Definition of a filtered sequence.
 **
 **/
class Sequence: public UnfilteredSequence
{
public:
    Sequence(const Spot &spot = Spot(), const std::string& index = "",
            unsigned int readNumber = 0, const std::string &data = "",
            const std::string &quality = "", bool passed = false);
    Sequence(const UnfilteredSequence &sequence, bool passed);
    Sequence(const Sequence &sequence);
    Sequence &operator=(const Sequence &sequence);
    bool operator==(const Sequence &sequence) const;
    bool operator!=(const Sequence &sequence) const;
    bool getPassed() const
    {
        return passed_;
    }
    void setPassed(bool passed)
    {
        passed_ = passed;
    }
private:
    bool passed_;
    friend std::ostream &operator<<(std::ostream &os, const Sequence &sequence);
    friend std::istream &operator>>(std::istream &is, Sequence &sequence);
};

/**
 ** @brief Write a complete object of type Sequence, including the delimiter.
 **/
std::ostream &operator<<(std::ostream &os, const Sequence &sequence);
/**
 ** @brief Read a complete object of type Sequence, including the delimiter.
 **/
std::istream &operator>>(std::istream &is, Sequence &sequence);

/**
 ** @brief An input stream specialized for sequences.
 **/
typedef Reader<Sequence> SequenceReader;

/**
 ** @brief An output stream specialized for sequences.
 **/
typedef Writer<Sequence> SequenceWriter;

} // namespace common
} // namespace casava

#endif // #ifndef CASAVA_COMMON_SEQUENCE_HH
