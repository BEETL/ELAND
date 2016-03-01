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
 ** \file Spot.hh
 **
 ** \brief Declaration of the data type used for a location on a tile
 **
 ** \author Come Raczy
 **/

#ifndef CASAVA_COMMON_SPOT_HH
#define CASAVA_COMMON_SPOT_HH

#include <iostream>
#include "common/Tile.hh"

namespace casava
{
namespace common
{

/**
 ** A location on a tile.
 **
 **/
class Spot
{
public:
    Spot(const Tile &tile = Tile(), int x = 0, int y = 0) :
        tile_(tile), x_(x), y_(y)
    {
    }
    Spot(const Spot &spot) :
        tile_(spot.tile_), x_(spot.x_), y_(spot.y_)
    {
    }
    Spot &operator=(const Spot &spot)
    {
        if (&spot != this)
        {
            tile_ = spot.tile_;
            x_ = spot.x_;
            y_ = spot.y_;
        }
        return *this;
    }
    bool operator==(const Spot &spot) const
    {
        return &spot == this || (tile_ == spot.tile_ && x_ == spot.x_ && y_
                == spot.y_);
    }
    bool operator!=(const Spot &spot) const
    {
        return !(spot == *this);
    }
    Tile getTile() const
    {
        return tile_;
    }
    Tile &getTile()
    {
        return tile_;
    }
    int getTileNumber() const
    {
        return tile_.getTileNumber();
    }
    int getX() const
    {
        return x_;
    }
    ;
    int getY() const
    {
        return y_;
    }
    ;
    void setTile(const Tile &tile)
    {
        tile_ = tile;
    }
    void setX(int x)
    {
        x_ = x;
    }
    void setY(int y)
    {
        y_ = y;
    }
private:
    Tile tile_;
    int x_;
    int y_;
    friend std::ostream &operator<<(std::ostream &os, const Spot &spot);
    friend std::istream &operator>>(std::istream &is, Spot &spot);
};

inline std::ostream &operator<<(std::ostream &os, const Spot &spot)
{
    os << spot.tile_;
    os.put('\t');
    putInteger(os, spot.x_);
    os.put('\t');
    putInteger(os, spot.y_);
    return os;
}

inline std::istream &operator>>(std::istream &is, Spot &spot)
{
    if (is >> spot.tile_)
    {
        is.ignore();
        getInteger(is, spot.x_, true);
        getInteger(is, spot.y_, false);
    }
    return is;
}

} // namespace common
} // namespace casava

#endif // CASAVA_COMMON_SPOT_HH
