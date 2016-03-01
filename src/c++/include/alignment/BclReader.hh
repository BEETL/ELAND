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
 ** \file BclReader.hh
 **
 ** \brief Reads BCL files and other associated files (filter, position).
 **
 ** \author Come Raczy
 **/

#ifndef CASAVA_ALIGNMENT_BCL_READER_HH
#define CASAVA_ALIGNMENT_BCL_READER_HH

#include <string>
#include <ios>
#include <fstream>
#include <vector>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "common/Exceptions.hh"

namespace casava
{
namespace alignment
{

namespace fs = boost::filesystem;

/**
 ** \brief Parallel ifstream providing sequences and qualities cluster by cluster.
 **/
class BclReader: protected boost::ptr_vector<std::ifstream>
{
public:
    typedef std::pair<std::string, std::string> Cluster;
    typedef Cluster RecordType;
    BclReader(const std::vector<fs::path> &pathList, const bool ignoreMissingBcl = false);
    unsigned int getClusterCount() const {return clusterCount_;}
    RecordType &get(RecordType &whereTo){
        getCluster(whereTo.first, whereTo.second);
        return whereTo;
    }
    std::string getDescription() const {
        return pathList_.begin()->string();
    }
    void getCluster(std::string &bases, std::string &quality);
    void getCluster(std::string &bases);
private:
    BclReader();
    BclReader(const BclReader &);
    BclReader &operator=(const BclReader &);
    const std::vector<fs::path> pathList_;
    const bool ignoreMissingBcl_;
    const unsigned int clusterCount_;
    unsigned int currentCluster_;
    unsigned int readClusterCount();
    unsigned int readClusterCount(const size_t index);
    unsigned char decode(const size_t value) const;
};

/**
 ** \brief ifstream providing filter information cluster by cluster.
 **/
class FiltersReader : boost::noncopyable
{
    typedef boost::error_info<struct tag_errmsg, std::string> errmsg_info;
public:
    typedef unsigned int RecordType;
    FiltersReader() : clusterCount_(0) {}
    FiltersReader(const fs::path &filePath, bool ctrlIncluded);
    unsigned int getClusterCount() const {return clusterCount_;}
    template <typename RetT, typename UnaryFunctionT> RetT executeGuarded(UnaryFunctionT f)
    {
        try
        {
            return f();
        }
        catch(boost::exception &e)
        {
            e << errmsg_info(" File: " + getDescription()
                    + " Cluster: " + boost::lexical_cast<std::string>(currentCluster_));
            throw;
        }
    }

    const RecordType &get(RecordType &whereTo)
    {
        if (clusterCount_ <= currentCluster_)
        {
            using boost::format;
            const format message = format("Attempt to read past the expected number of clusters (%d). File: %s")
                % clusterCount_ % getDescription();
            BOOST_THROW_EXCEPTION(casava::common::IoException(EINVAL, message.str()));
        }

        ++currentCluster_;

        whereTo = executeGuarded<unsigned int>(
                boost::bind(&FiltersReaderImpl::doGetNextFilter, readerImplPtr_, boost::ref(is_)));
        return whereTo;
    }

    const std::string getDescription() const {
        return filePath_.string();
    }

    class FiltersReaderImpl
    {
    public:
        virtual ~FiltersReaderImpl(){;}
        virtual unsigned int doReadClusterCount(std::istream &is) = 0;
        virtual unsigned int doGetNextFilter(std::istream& is) = 0;
    };

private:
    const fs::path filePath_;
    std::ifstream is_;

    boost::shared_ptr<FiltersReaderImpl> readerImplPtr_;
    const unsigned int clusterCount_;
    unsigned int currentCluster_;

    boost::shared_ptr<FiltersReaderImpl> createReader();
};

/**
 ** \brief ifstream providing position information cluster by cluster.
 **/
class PositionsReader
{
public:
    typedef std::pair<float, float> FloatPosition;
    typedef FloatPosition RecordType;

    typedef std::pair<int, int> Position;

    PositionsReader(const fs::path &filePath, std::ios_base::openmode mode)
    : filePath_(filePath), is_(filePath_.string().c_str(), mode), currentCluster_(0) {}

    RecordType& get(RecordType &whereTo) {
        doGetFloatPosition(whereTo);
        return whereTo;
    }

    std::string getDescription() const {
        return filePath_.string();
    }

    FloatPosition getFloatPosition()
    {
        FloatPosition ret;
        return getFloatPosition(ret);
    }


    FloatPosition &getFloatPosition(FloatPosition &whereTo);
    Position getPosition();
    virtual ~PositionsReader() {}
    static PositionsReader *create(const fs::path &filePath, unsigned int clusterCount);
    virtual unsigned int getClusterCount() const = 0;
protected:
    const fs::path filePath_;
    std::ifstream is_;
    unsigned int currentCluster_;
private:
    virtual FloatPosition &doGetFloatPosition(FloatPosition &whereTo) = 0;
};

class PositionsReaderBinary: public PositionsReader
{
public:
    PositionsReaderBinary(const fs::path &filePath);
    unsigned int getClusterCount() const {return clusterCount_;}
private:
    const unsigned int clusterCount_;
    virtual FloatPosition &doGetFloatPosition(FloatPosition &whereTo);
    unsigned int readClusterCount();
};

class PositionsReaderCompressed: public PositionsReader
{
public:
    PositionsReaderCompressed(const fs::path &filePath, unsigned int clusterCount);
    unsigned int getClusterCount() const {return clusterCount_;}
private:
    virtual FloatPosition &doGetFloatPosition(FloatPosition &whereTo);
    const unsigned int clusterCount_;
    unsigned int blocksCount_;
    unsigned int currentBlock_;
    unsigned char currentBlockUnreadClusters_;

    // Constants must match those in $Illumina.RTA/Dev/Trunk/Src/Shared/ILMNcommon/Common/FloatPoint.cs
    static const int blockSize = 25;
    static const int imageWidth = 2048;
    static const int blocksPerlLine = (imageWidth + blockSize - 1) / blockSize;
    static const int imageHeight = 20000;
};

class PositionsReaderText: public PositionsReader
{
public:
    PositionsReaderText(const fs::path &filePath, unsigned int clusterCount);
    unsigned int getClusterCount() const {return clusterCount_;}
private:
    virtual FloatPosition &doGetFloatPosition(FloatPosition &whereTo);
    const unsigned int clusterCount_;
    std::string lineBuffer_;
};

} // namespace alignment
} // namespace casava

#endif // #ifndef CASAVA_ALIGNMENT_BCL_READER_HH
