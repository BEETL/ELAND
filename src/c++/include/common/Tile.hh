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
 ** \file Tile.hh
 **
 ** \brief Declaration of the data type used to identify a tile
 **
 ** \author Come Raczy
 **/

#ifndef CASAVA_COMMON_TILE_HH
#define CASAVA_COMMON_TILE_HH

#include <string>
#include <iostream>
#include "common/FastIo.hh"

namespace casava
{
namespace common
{

/**
 ** Identification of a tile.
 **/
class Tile
{
public:
    Tile() :
        machineName_(""), runNumber_(0), laneNumber_(0), tileNumber_(0)
    {
    }
    Tile(const std::string &machineName, unsigned int runNumber,
            unsigned int laneNumber, unsigned int tileNumber) :
        machineName_(machineName), runNumber_(runNumber), laneNumber_(
                laneNumber), tileNumber_(tileNumber)
    {
    }
    Tile(const Tile &tile) :
        machineName_(tile.machineName_), runNumber_(tile.runNumber_),
                laneNumber_(tile.laneNumber_), tileNumber_(tile.tileNumber_)
    {
    }
    Tile &operator=(const Tile &tile)
    {
        if (&tile != this)
        {
            machineName_ = tile.machineName_;
            runNumber_ = tile.runNumber_;
            laneNumber_ = tile.laneNumber_;
            tileNumber_ = tile.tileNumber_;
        }
        return *this;
    }
    inline bool operator==(const Tile &tile) const
    {
        return &tile == this || (machineName_ == tile.machineName_
                && runNumber_ == tile.runNumber_ && laneNumber_
                == tile.laneNumber_ && tileNumber_ == tile.tileNumber_);
    }
    inline bool operator!=(const Tile &tile) const
    {

        return !(tile == *this);
    }
    const std::string &getMachineName() const
    {
        return machineName_;
    }
    unsigned int getRunNumber() const
    {
        return runNumber_;
    }
    unsigned int getLaneNumber() const
    {
        return laneNumber_;
    }
    unsigned int getTileNumber() const
    {
        return tileNumber_;
    }
    void setMachineName(const std::string &machineName)
    {
        machineName_ = machineName;
    }
    void setRunNumber(unsigned int runNumber)
    {
        runNumber_ = runNumber;
    }
    void setLaneNumber(unsigned int laneNumber)
    {
        laneNumber_ = laneNumber;
    }
    void setTileNumber(unsigned int tileNumber)
    {
        tileNumber_ = tileNumber;
    }
private:
    std::string machineName_;
    unsigned int runNumber_;
    unsigned int laneNumber_;
    unsigned int tileNumber_;
    friend std::ostream &operator<<(std::ostream &os, const Tile &tile);
    friend std::istream &operator>>(std::istream &is, Tile &tile);
};

inline std::ostream &operator<<(std::ostream &os, const Tile &tile)
{
    os.write(tile.machineName_.c_str(), tile.machineName_.size());
    os.put('\t');
    putUnsignedInteger(os, tile.runNumber_);
    os.put('\t');
    putUnsignedInteger(os, tile.laneNumber_);
    os.put('\t');
    putUnsignedInteger(os, tile.tileNumber_);
    return os;
}

inline std::istream &operator>>(std::istream &is, Tile &tile)
{
    std::stringbuf sb;
    if (is.get(sb, '\t'))
    {
        is.ignore();
        tile.machineName_ = sb.str();
        getUnsignedInteger(is, tile.runNumber_, true);
        getUnsignedInteger(is, tile.laneNumber_, true);
        getUnsignedInteger(is, tile.tileNumber_, false);
    }
    return is;
}

} // namespace common
} // namespace casava

#endif // #ifndef CASAVA_COMMON_TILE_HH
