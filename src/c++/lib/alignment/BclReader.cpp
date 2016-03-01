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
 ** \file BclReader.cpp
 **
 ** \brief Reads BCL files and other associated files (filter, position).
 **
 ** \author Come Raczy
 **/

#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "alignment/BclReader.hh"
#include "common/FastIo.hh"
#include "common/Exceptions.hh"

namespace casava
{
namespace alignment
{

namespace cc = casava::common;

BclReader::BclReader(const std::vector<fs::path> &pathList, const bool ignoreMissingBcl)
    : pathList_(pathList)
    , ignoreMissingBcl_(ignoreMissingBcl)
    , clusterCount_(readClusterCount())
    , currentCluster_(0)
{

}

unsigned int BclReader::readClusterCount(const size_t index)
{
    std::istream &is = at(index);
    const unsigned int clusterCount = cc::readUnsignedInteger<4>(is);
    if (!is)
    {
        using boost::format;
        const fs::path &path = pathList_[index];
        const format message = format("Failed to read number of clusters from %s") % path;
        if (ignoreMissingBcl_)
        {
            std::cerr << "WARNING: " << message.str() << std::endl;
            return 0;
        } else {
            BOOST_THROW_EXCEPTION(cc::IoException(errno, message.str()));
        }
    }
    return clusterCount;
}

unsigned int BclReader::readClusterCount()
{
    if (pathList_.empty())
    {
        return 0;
    }
    // open the files
    BOOST_FOREACH(const fs::path &path, pathList_)
    {
        using boost::format;
        if (!fs::exists(path))
        {
            const format message = format("File %s does not exist") % path;
            if (ignoreMissingBcl_)
            {
                std::cerr << "WARNING: " << message.str() << std::endl;
            } else {
                BOOST_THROW_EXCEPTION(cc::IoException(errno, message.str()));
            }
        }
        push_back(new std::ifstream(path.string().c_str()));
        if (!back())
        {
            const format message = format("Couldn't open BCL file %s.") % path;
            if (ignoreMissingBcl_)
            {
                std::cerr << "WARNING: " << message.str() << std::endl;
            } else {
                BOOST_THROW_EXCEPTION(cc::IoException(errno, message.str()));
            }
        }
    }
    size_t current = 0;
    const unsigned int clusterCount = readClusterCount(current);
    while (size() != ++current)
    {
        const unsigned int n = readClusterCount(current);
        if (clusterCount != n)
        {
            using boost::format;
            const fs::path path = pathList_[current];
            const format message = format("Incorrect number of clusters in %s: expected %d: got %d") % path % clusterCount % n;
            if (ignoreMissingBcl_)
            {
                std::cerr << "WARNING: " << message.str() << std::endl;
            } else {
                BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, message.str()));
            }
        }
    }
    return clusterCount;
}

unsigned char BclReader::decode(const size_t value) const
{
    static const unsigned char BaseUpperCase[4] = { 'A', 'C', 'G', 'T'};
    assert(4 > value);
    return BaseUpperCase[value];
}

void BclReader::getCluster(std::string &bases, std::string &qualities)
{
    using boost::format;
    if (clusterCount_ <= currentCluster_)
    {
        const format message = format("Method 'getCluster' called more than %d times") % clusterCount_;
        BOOST_THROW_EXCEPTION(cc::PreConditionException(message.str()));
    }
    ++currentCluster_;
    bases.clear();
    qualities.clear();
    for(boost::ptr_vector<std::ifstream>::iterator is = begin(); end() != is; ++is)
    {
        int c = is->get();
        if (is->fail())
        {
            const format message = format("Failed to read BCL file %s.") % pathList_[is - begin()];
            if (!ignoreMissingBcl_)
            {
                BOOST_THROW_EXCEPTION(cc::IoException(errno, message.str()));
            }
            c = 0;
        }
        if (is->eof())
        {
            const format message = format("Unexpected EOF for BCL file %s") % pathList_[is - begin()];
            if (!ignoreMissingBcl_)
            {
                BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, message.str()));
            }
            c = 0;
        }
        const unsigned int quality = ((c & 0xfc) >> 2);
        bases.push_back(quality ? decode(c & 0x3) : 'N');
        qualities.push_back(quality ? quality + 64 : 66);
    }
}

/**
 ** \brief Pre-RTA 1.9 filter file reader. Format: first four bytes are number of clusters,
 ** followed by one byte per cluster filtering values.
 **/
class FiltersReaderImpl_8bits : public FiltersReader::FiltersReaderImpl, boost::noncopyable
{
private:
    virtual unsigned int doReadClusterCount(std::istream &is);
    virtual unsigned int doGetNextFilter(std::istream& is);
};

unsigned int FiltersReaderImpl_8bits::doReadClusterCount(std::istream &is)
{
    using boost::format;
    // 4 bytes little endian
    const unsigned int header = cc::readUnsignedInteger<4>(is);
    if (!is) {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "Unable to read filter file header."));
    }
    if (0x00000000 != header)
    {
        // this is version 0:
        return header;
    }   // else it is version 3+:
    const unsigned int filterVersion = cc::readUnsignedInteger<4>(is);
    if (!is) {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "Unable to read filter file file version."));
    }

    const unsigned int earliestVersion = 0x00000003;
    if (earliestVersion > filterVersion)
    {
        const format message = format("Unexpected version byte %#04x found in filters file header."
                                      " Expected version > %#04x.") % filterVersion % (earliestVersion - 1);
        BOOST_THROW_EXCEPTION(cc::UnsupportedVersionException(message.str()));
    }

    const unsigned int clusterCount = cc::readUnsignedInteger<4>(is);
    if (!is) {
        const format message = format("Unable to read cluster count from filter file (version %#04x).") % filterVersion;
        BOOST_THROW_EXCEPTION(cc::IoException(errno,message.str()));
    }
    return clusterCount;
}

unsigned int FiltersReaderImpl_8bits::doGetNextFilter(std::istream &is)
{
    // 1 byte per cluster
    const unsigned int filter = static_cast<unsigned int>(is.get());
    if (!is) {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "Failed to read filter value."));
    }
    return filter;
}


/**
 ** \brief RTA 1.9/1.10 filter file reader. Recognizes v1, v2, and v3 filter files depending on header contents
 **/
class FiltersReaderImpl_16bits : public FiltersReader::FiltersReaderImpl, boost::noncopyable
{
private:
    virtual unsigned int doReadClusterCount(std::istream &is);
    virtual unsigned int doGetNextFilter(std::istream& is);
};

unsigned int FiltersReaderImpl_16bits::doReadClusterCount(std::istream &is)
{
    using boost::format;
    const unsigned int header = cc::readUnsignedInteger<4>(is);
    if (!is) {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "Unable to read filter/control file header."));
    }
    if (0x00000000 != header)
    {
        // this is version 1:
        return header;
    }   // else it is version 2+:
    const unsigned int filterVersion = cc::readUnsignedInteger<4>(is);
    if (!is) {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "Unable to read filter/control file file version."));
    }

    const unsigned int earliestVersion = 0x00000002;
    if (earliestVersion > filterVersion)
    {
        const format message = format("Unexpected version byte %#04x found in filters file header."
                                      " Expected version > %#04x.") % filterVersion % (earliestVersion - 1);
        BOOST_THROW_EXCEPTION(cc::UnsupportedVersionException(message.str()));
    }

    const unsigned int clusterCount = cc::readUnsignedInteger<4>(is);
    if (!is) {
        const format message = format("Unable to read cluster count from filter/control file (version %#04x).") % filterVersion;
        BOOST_THROW_EXCEPTION(cc::IoException(errno,message.str()));
    }
    return clusterCount;
}


unsigned int FiltersReaderImpl_16bits::doGetNextFilter(std::istream &is)
{
    // 2 bytes per cluster
    const unsigned int filter = cc::readUnsignedInteger<2>(is);
    if (!is) {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "Failed to read filter value."));
    }
    return filter;
}


FiltersReader::FiltersReader(const fs::path &filePath, bool ctrlIncluded)
    : filePath_(filePath)
    , is_(filePath_.string().c_str(),  std::ios_base::in | std::ios_base::binary)
    , readerImplPtr_(ctrlIncluded
            ? boost::shared_ptr<FiltersReaderImpl>(new FiltersReaderImpl_16bits)
            : boost::shared_ptr<FiltersReaderImpl>(new FiltersReaderImpl_8bits))
    , clusterCount_( 0 == filePath_.string().length()
                   ? 0
                   : executeGuarded<unsigned int>(
                            boost::bind(&FiltersReaderImpl::doReadClusterCount, readerImplPtr_, boost::ref(is_)) )
                   )
    , currentCluster_(0)
{
}


PositionsReader *PositionsReader::create(const fs::path &filePath, unsigned int clusterCount)
{
    PositionsReader *result = 0;
    const std::string extension = fs::extension(filePath);
    if (".locs" == extension)
    {
        result = new PositionsReaderBinary(filePath);
    }
    else if (".clocs" == extension)
    {
        result = new PositionsReaderCompressed(filePath, clusterCount);
    }
    else if (".txt" == extension)
    {
        result = new PositionsReaderText(filePath, clusterCount);
    }
    else
    {
        using boost::format;
        const format message = format("Unknown format for the positions file %s: supported formats are 'locs', 'clocs' and '_pos.txt': %s") % filePath % extension;
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, message.str()));
    }
    return result;
}

PositionsReaderBinary::PositionsReaderBinary(const fs::path &filePath)
    : PositionsReader(filePath, std::ios_base::in | std::ios_base::binary)
    , clusterCount_(readClusterCount())
{
}

unsigned int PositionsReaderBinary::readClusterCount()
{
    // first 8 bytes are unused
    is_.seekg(8, std::ios_base::beg);
    // 4 bytes little endian
    const unsigned int clusterCount = cc::readUnsignedInteger<4>(is_);
    if (!is_)
    {
        using boost::format;
        const format message = format("Failed to read number of clusters from %s.") % filePath_;
        BOOST_THROW_EXCEPTION(cc::IoException(errno, message.str()));
    }
    return clusterCount;
}

PositionsReader::FloatPosition &PositionsReaderBinary::doGetFloatPosition(FloatPosition &whereTo)
{
    float *f[] = {&whereTo.first, &whereTo.second};
    for (unsigned int i = 0; 2 > i; ++i)
    {
        *(f[i]) = cc::readDecimalNumber<4>(is_);
        if (!is_)
        {
            using boost::format;
            static const char coordinates[] = {'X', 'Y'};
            const format message = format("Failed to read %c coordinate for cluster %d in file %s")
                % coordinates[i] % currentCluster_ % filePath_;
            BOOST_THROW_EXCEPTION(cc::IoException(errno, message.str()));
        }
    }
    return whereTo;
}

PositionsReader::Position PositionsReader::getPosition()
{
    const FloatPosition floatPosition = this->getFloatPosition();
    return PositionsReader::Position(static_cast<int>(roundf(1000.0f + 10.0f*floatPosition.first)),
                                      static_cast<int>(roundf(1000.0f + 10.0f*floatPosition.second)));
}

PositionsReader::FloatPosition &PositionsReader::getFloatPosition(FloatPosition &whereTo)
{
    if (this->getClusterCount() <= currentCluster_)
    {
        using boost::format;
        const format message = format("Reading more positions than available in positions file %s: %d")
            % filePath_ % getClusterCount();
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, message.str()));
    }
    ++currentCluster_;
    this->doGetFloatPosition(whereTo);
    return whereTo;
}

PositionsReaderCompressed::PositionsReaderCompressed(const fs::path &filePath, unsigned int clusterCount)
    : PositionsReader(filePath, std::ios_base::in | std::ios_base::binary)
    , clusterCount_(clusterCount), blocksCount_(0), currentBlock_(0), currentBlockUnreadClusters_(0)
{
    const unsigned char clocsVersion = cc::readUnsignedInteger<1>(is_);
    if (!is_) {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "Failed to read version from clocs file " + filePath_.string()));
    }
    if (0x01 != clocsVersion) {
        BOOST_THROW_EXCEPTION(cc::UnsupportedVersionException((boost::format("Unexpected version byte %#04x"
                " found in clocs file header. Expected: 0x01, File: %s ") % int(clocsVersion) % filePath_.string()).str()));
    }
    blocksCount_ = cc::readUnsignedInteger<4>(is_);
    if (!is_) {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "Failed to read number of tiles from clocs file " + filePath_.string()));
    }

    currentBlockUnreadClusters_ = cc::readUnsignedInteger<1>(is_);
    if (is_.fail()) {
        BOOST_THROW_EXCEPTION(cc::IoException(errno,
                (boost::format("Failed to read number of tile fist clusters from clocs file %s."
                               "total tiles: %d" ) % filePath_.string() % blocksCount_).str()));
    }

}

PositionsReader::FloatPosition &PositionsReaderCompressed::doGetFloatPosition(FloatPosition &whereTo)
{

    while (!currentBlockUnreadClusters_--)
    {
        if (blocksCount_ <= ++currentBlock_) {
            BOOST_THROW_EXCEPTION(cc::IoException(EINVAL,
                    (boost::format("Attempt to read clocs file past the expected number of tiles. "
                                   "File: %s, Current cluster: %d, Total blocks: %d")
                                 % filePath_.string() % currentCluster_ % blocksCount_).str()));
        }

        currentBlockUnreadClusters_ = cc::readUnsignedInteger<1>(is_);
        if (is_.fail()) {
            BOOST_THROW_EXCEPTION(cc::IoException(errno,
                    (boost::format("Failed to read number of tile clusters from clocs file %s."
                                   "Current/total tile: %d/%d" ) % filePath_.string() % currentBlock_ % blocksCount_).str()));
        }
    }

    unsigned char dx = cc::readUnsignedInteger<1>(is_);
    unsigned char dy = cc::readUnsignedInteger<1>(is_);

    if (is_.fail()) {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "Failed to read position from clocs file " + filePath_.string()));
    }

    whereTo.first = float(blockSize) * (currentBlock_ % blocksPerlLine) + dx / 10.0;
    whereTo.second = float(blockSize) * (currentBlock_ / blocksPerlLine) + dy / 10.0;

    return whereTo;
}

PositionsReaderText::PositionsReaderText(const fs::path &filePath, unsigned int clusterCount)
    : PositionsReader(filePath, std::ios_base::in)
    , clusterCount_(clusterCount)
{
}

PositionsReader::FloatPosition &PositionsReaderText::doGetFloatPosition(FloatPosition &whereTo)
{
    static const std::string whiteSpaces = " \t\n\r\v";
    lineBuffer_.clear();
    if (getline(is_, lineBuffer_))
    {
        while(std::find(whiteSpaces.begin(), whiteSpaces.end(), *lineBuffer_.rbegin()) != whiteSpaces.end())
        {
            lineBuffer_.erase(lineBuffer_.size() - 1, 1);
        }
        using boost::format;
        float *f[] = {&whereTo.first, &whereTo.second};
        char *end = 0;
        for (unsigned int i = 0; 2 > i; ++i)
        {
            const char *begin = (end ? end : lineBuffer_.c_str());
            *(f[i]) = strtof(begin, &end);
            if (begin == end)
            {
                is_.clear(std::ios_base::failbit);
                static const char coordinates[] = {'X', 'Y'};
                const format message = format("Failed to read %c coordinate for cluster %d in file %s: %s")
                    % coordinates[i] % currentCluster_ % filePath_ % lineBuffer_;
                BOOST_THROW_EXCEPTION(cc::IoException(errno, message.str()));
            }
        }
        if (*end)
        {
            is_.clear(std::ios_base::failbit);
            const format message = format("Unexpected characters after Y coordinate for cluster %d in file %s: %s")
                % currentCluster_ % filePath_ % lineBuffer_;
            BOOST_THROW_EXCEPTION(cc::IoException(errno, message.str()));
        }
    }
    else
    {
        if (is_.fail()) {
            BOOST_THROW_EXCEPTION(cc::IoException(errno, "Failed to read position from pos file " + filePath_.string()));
        }
    }
    return whereTo;
}

} // namespace alignment
} // namespace casava

