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
 ** \file FilteringStreams.hh
 **
 ** \brief I/O API for sequence data type (qseq files).
 **
 ** Stream-based API based on the Reader and Writer concepts defined
 ** in the "FilteringStreams.hh" header. The API factors out an
 ** UnfilteredSequence to ease reuse, for instance for Alignment data.
 **
 ** \author Come Raczy
 **/

#ifndef CASAVA_COMMON_FILTERING_STREAMS_HH
#define CASAVA_COMMON_FILTERING_STREAMS_HH

#include <boost/iostreams/filtering_stream.hpp>
#include <string>

#include "common/Compression.hh"

namespace casava
{
namespace common
{

/**
 ** The base input stream file to use for all data types.
 **
 ** The appropriate algorithm to load the data will depend on the
 ** underlying type of the data: text or binary.
 **
 ** Usage:
 ** const ga::common::Filter &decompress =
 **     ga::common::CompressionFactory::get(std::string(argv[1]));
 ** ga::common::Ifstream ifs(std::string(argv[1]), decompress);
 ** std::string line;
 ** while (getline(ifs, line))
 ** {
 **     // process the line
 ** }
 **/
class Ifstream: public boost::iostreams::filtering_istream
{
public:
    Ifstream(const std::string &filePath, const Filter &decompress);
    bool is_open();
private:
    Ifstream();
    Ifstream(const Ifstream &);
    Ifstream &operator=(const Ifstream &);
};

/**
 ** The base output stream file to use for all data types.
 **
 ** The appropriate algorithm to write the data will depend on the
 ** underlying type of the data: text or binary.
 **/
class Ofstream: public boost::iostreams::filtering_ostream
{
public:
    Ofstream(const std::string &filePath, const Filter &compress,
            ios_base::openmode mode = ios_base::out);
    Ofstream(const Filter &compress);
    void open(const std::string &filePath, ios_base::openmode mode =
            ios_base::out);
    void close();
    bool is_open();
private:
    Ofstream();
    Ofstream(const Ofstream &);
    Ofstream &operator=(const Ofstream &);

    void preventCorruptGzipFiles();
};

/**
 ** \brief A Reader class specialized for one specific type of input.
 **
 ** \tparam T must have a default constructor and an input stream operator >>
 **/
template<typename T>
class Reader: public Ifstream
{
public:
    /**
     ** @brief Construct an object of the Reader class.
     **
     ** Opens the file and set the compression filter.
     **/
    Reader(const std::string &filePath, const Filter &decompress) :
        Ifstream(filePath, decompress)
    {
    }
    /**
     ** @brief Read the first value available from the stream and
     ** stores it at the specified location.
     **
     ** Note: this method hides all the variants of istream::get
     **
     ** @return *this
     **
     **/
    Reader &get(T &value)
    {
        *this >> value;
        return *this;
    }
    /**
     ** @brief Read all the values available from the stream and
     ** stores it int the vector valueList.
     **
     ** Failure to read can be tested with the member fail.
     **
     ** Note: this method hides all the variants of istream::gread
     **/
    Reader &read(std::vector<T> &valueList)
    {
        T value;
        while (get(value))
        {
            valueList.push_back(value);
        }
        return *this;
    }

    /**
     ** @brief Read a block of n values from the stream and
     ** stores it int the vector valueList.
     **
     ** If the End-of-File is reached before n values have been
     ** read, the vector will contain all the elements read until it,
     ** and the failbit and eofbit will be set (which can be checked
     ** with members fail and eof respectively).
     **
     ** Note: this method hides all the variants of istream::gread
     **/
    Reader &read(std::vector<T> &valueList, unsigned int n)
    {
        T value;
        while (0 != n && get(value))
        {
            valueList.push_back(value);
            --n;
        }
        return *this;
    }

};

/**
 ** \brief A Writer class specialized for one specific type of output.
 **
 ** \tparam T must have a default constructor and an output stream operator <<
 **/
template<typename T>
class Writer: public Ofstream
{
public:
    /**
     ** @brief Construct an object of the Writer class.
     **
     ** Opens the file and set the compression filter.
     **/
    Writer(const std::string &filePath, const Filter &compress) :
        Ofstream(filePath, compress)
    {
    }

    /**
     ** @brief Construct an object of the Writer class.
     ** without opening the file.
     **
     ** Use open before you write anything to this stream
     **/
    Writer(const Filter &compress) : Ofstream(compress) {}

    /**
     ** @brief Write one value into the stream.
     **
     ** Note: this method hides all the variants of istream::put
     **
     ** @return *this
     **
     **/
    Writer &put(const T &value)
    {
        *this << value;
        return *this;
    }
    /**
     ** @brief Write a vector of values into the stream.
     **
     ** Note: this method hides all the variants of istream::write
     **/
    Writer &write(const std::vector<T> &valueList)
    {
        for (typename std::vector<T>::const_iterator value = valueList.begin(); valueList.end()
                != value; ++value)
        {
            this->put(*value);
        }
        return *this;
    }
};

} // namespace common
} // namespace casava

#endif // #ifndef CASAVA_COMMON_FILTERING_STREAMS_HH
