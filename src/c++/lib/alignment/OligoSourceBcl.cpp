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
 ** \file OligoSourceBcl.cpp
 **
 ** \brief Commonly used functions and definitions
 **
 ** \author Tony Cox
 **/

//#include <iostream>
//#include <ctime>
//#include <dirent.h>
#include <boost/foreach.hpp>
//#include <boost/format.hpp>

#include "alignment/OligoSourceBcl.hh"
#include "eland_ms/ElandConstants.hh"

namespace casava
{
namespace alignment
{

/*****************************************************************************/
// **
// ** Function definitions for OligoSourceBcl
// **
/*****************************************************************************/

OligoSourceBcl::OligoSourceBcl(const std::vector<fs::path> &bclDirectoryList,
                               const std::vector<fs::path> &barcodeDirectoryList,
                               const fs::path &positionsDirectory,
                               const fs::path &filterDirectory,
                               const boost::format &positionsFileNameFormat,
                               const std::string &machineName,
                               const unsigned int runNumber,
                               const unsigned int lane,
                               const std::vector<unsigned int> tileList,
                               const unsigned int readNumber)
    : format_(true)
    , bclDirectoryList_(bclDirectoryList)
    , barcodeDirectoryList_(barcodeDirectoryList)
    , positionsDirectory_(positionsDirectory)
    , filterDirectory_(filterDirectory)
    , positionsFileNameFormat_(positionsFileNameFormat)
    , lane_(lane)
    , tileList_(tileList)
    , currentTile_(tileList_.begin())
    , currentCluster_(0)
    , currentClusterInTile_(0)
    , bclReader_(0)
    , barcodeReader_(0)
    , positionsReader_(0)
    , filtersReader_(0)
{
    sequence_.setMachineName(machineName);
    sequence_.setRunNumber(runNumber);
    sequence_.setLaneNumber(lane);
    sequence_.setReadNumber(readNumber);
    sequence_.getIndex().reserve(barcodeDirectoryList_.size());
    sequence_.getIndex().push_back('0');
    sequence_.getData().reserve(bclDirectoryList_.size());
    sequence_.getQuality().reserve(bclDirectoryList_.size());
    if (tileList_.end() != currentTile_)
    {
        initializeNewTile();
    }
}

void OligoSourceBcl::rewind()
{
    sequenceName_.clear();
    currentTile_ = tileList_.begin();
    currentCluster_ = 0;
    if (tileList_.end() != currentTile_)
    {
        initializeNewTile();
    }
}

const casava::common::Sequence &OligoSourceBcl::getLastSequence(bool& isValid) const
{
    isValid = (0 == bclReader_ || 0 == filtersReader_ || 0 == positionsReader_);
    return sequence_;
}

const casava::common::Sequence &OligoSourceBcl::getNextSequenceSelect(bool& isValid,
                                                                      const bool,
                                                                      const bool)
{
    isValid = getCluster();
    while (bclReader_ && !isValid)
    {
        isValid = getCluster();
    }
    return sequence_;
}

bool OligoSourceBcl::getCluster()
{
    sequenceName_.clear();
    if (0 == bclReader_ || bclReader_->getClusterCount() <= currentClusterInTile_)
    {
        initializeNewTile(); // skips empty tiles
    }
    if (bclReader_)
    {
        bclReader_->getCluster(sequence_.getData(), sequence_.getQuality());
        if (barcodeReader_)
        {
            std::string barcodeQuality;
            barcodeReader_->getCluster(sequence_.getIndex(), barcodeQuality);
        }
        typedef casava::alignment::PositionsReader::Position Position;
        Position position = positionsReader_->getPosition();
        sequence_.setX(position.first);
        sequence_.setY(position.second);
        unsigned int filterValue(0);
        sequence_.setPassed(filtersReader_->get(filterValue));
        ++currentCluster_;
        ++currentClusterInTile_;
        return (isNoMask_ || mask_[currentCluster_]);
    }
    return false;
}

const char * OligoSourceBcl::getNextOligo()
{
    bool isValid = true;
    getNextSequence(isValid);
    return (isValid ? sequence_.getData().c_str() : 0);
}

const char * OligoSourceBcl::getLastOligo() const
{
    bool isValid = true;
    getLastSequence(isValid);
    return (isValid ? sequence_.getData().c_str() : 0);
}

const char *OligoSourceBcl::getLastName()
{
    if (sequenceName_.empty())
    {
        sequenceName_ = (boost::format(">%s_%04u:%u:%u:%d:%d#%s/%u")
                       % sequence_.getSpot().getTile().getMachineName()
                       % sequence_.getSpot().getTile().getRunNumber()
                       % sequence_.getSpot().getTile().getLaneNumber()
                       % sequence_.getSpot().getTile().getTileNumber()
                       % sequence_.getSpot().getX()
                       % sequence_.getSpot().getY()
                       % sequence_.getIndex()
                       % sequence_.getReadNumber()).str();
    }
    return sequenceName_.c_str();
}

void OligoSourceBcl::initializeNewTile()
{
    using namespace std;
    using boost::format;
    delete bclReader_;
    bclReader_ = 0;
    delete barcodeReader_;
    barcodeReader_ = 0;
    delete positionsReader_;
    positionsReader_ = 0;
    delete filtersReader_;
    filtersReader_ = 0;
    if (tileList_.end() != currentTile_)
    {
        sequence_.setTileNumber(*currentTile_);
        const std::string bclFileName = (boost::format("s_%d_%d.bcl") % lane_ % (*currentTile_)).str();
        std::vector<fs::path> bclFileList;
        BOOST_FOREACH(const fs::path &d, bclDirectoryList_) {bclFileList.push_back(d / bclFileName);}
        bclReader_ = new casava::alignment::BclReader(bclFileList);
        if (!barcodeDirectoryList_.empty())
        {
            std::vector<fs::path> barcodeFileList;
            BOOST_FOREACH(const fs::path &d, bclDirectoryList_) {barcodeFileList.push_back(d / bclFileName);}
            barcodeReader_ = new casava::alignment::BclReader(barcodeFileList);
        }
        const std::string positionsFileName = (boost::format(positionsFileNameFormat_) % lane_ % (*currentTile_)).str();
        positionsReader_ = casava::alignment::PositionsReader::create(positionsDirectory_ / positionsFileName, bclReader_->getClusterCount());
        const std::string filterFileName = (boost::format("s_%d_%04d.filter") % lane_ % (*currentTile_)).str();
        filtersReader_ = new casava::alignment::FiltersReader(filterDirectory_ / filterFileName,
                filterDirectory_ != bclDirectoryList_[0].parent_path());
        ++currentTile_;
        currentClusterInTile_ = 0;
    }
    // skip empty tiles
    if (bclReader_ && 0 == bclReader_->getClusterCount())
    {
        initializeNewTile();
    }
}

} //namespace alignment
} //namespace casava
