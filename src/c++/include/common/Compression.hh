/**
 ** Copyright (c) 2007-2010 Illumina, Inc.
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
 ** \file Compression.hh
 **
 ** \brief Type-safe support for the different types of compression.
 **
 ** Filtered-based API.
 **
 ** \author Come Raczy
 **/

#ifndef CASAVA_COMMON_COMPRESSION_HH
#define CASAVA_COMMON_COMPRESSION_HH

#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/static_assert.hpp>

#include "config.h"
namespace boost {namespace iostreams {class gzip_params;}}
namespace boost {namespace iostreams {class bzip2_params;}}

namespace casava
{
namespace common
{

/**
 ** Base class for all types of filters. Implements Filter "none".
 **
 **/
class Filter
{
public:
    virtual void operator()(boost::iostreams::filtering_ostream &) const
    {
    }
    virtual void operator()(boost::iostreams::filtering_istream &) const
    {
    }
protected:
    Filter()
    {
    }
    Filter(const Filter &)
    {
    }
    Filter &operator=(const Filter &)
    {
        return *this;
    }
    virtual ~Filter()
    {
    }
};

class CompressionFactory;

/**
 ** Base class for all types of compression.
 **
 ** Specializations: derived classes must provide alternative
 ** implementations of the operator(). Ideally, these implementations
 ** should just push the right compression into the compressioning stream.
 **
 ** Usage:
 **   const Compression &compress = CompressionFactory::get("gzip");
 **   boost::iostreams::filtering_ostream out;
 **   compress(out)
 **   boost::iostreams::filtering_istream in;
 **   compress(in)
 **/
class Compression: public Filter
{
public:
    //virtual void operator()(boost::iostreams::filtering_ostream &) const {};
    //virtual void operator()(boost::iostreams::filtering_istream &) const {};
    virtual const std::string getFileNameExtension() const
    {
        return "";
    }
    virtual const std::string getName() const
    {
        return "none";
    }
protected:
    Compression()
    {
    }
    Compression(const Compression &) :
        Filter()
    {
    }
    Compression &operator=(const Compression &)
    {
        return *this;
    }
    virtual ~Compression()
    {
    }
    friend class CompressionFactory;
};

/**
 **
 **/
class UnsupportedCompressionException: public std::logic_error
{
public:
    UnsupportedCompressionException(const std::string &algorithm);
};

/**
 * @brief Base template for resolving compression parameter structures into
 * a corresponding filter object. See concrete specializations for the list
 * of supported compressions below.
 *
 * Consider using casava::common::makeCompressionFilter instead of this template.
 */
template<typename ParamsT>
class CompressionFilter
{
    BOOST_STATIC_ASSERT(sizeof(ParamsT) == 0 && "Only specializations of this class are allowed to instantiate");
};

// Note: we're forced to put the implementations in cpp because of the following very noisy warnings:
// /home/rpetrovski/workspace/builds/BF00836/opt/bootstrap/include/boost/iostreams/filter/gzip.hpp:674: warning: overflow in implicit constant conversion
#ifdef HAVE_ZLIB
template<>
class CompressionFilter<boost::iostreams::gzip_params> : public Compression
{
public:
    CompressionFilter(const boost::iostreams::gzip_params &params) : params_(params){}
private:
    const boost::iostreams::gzip_params &params_;
    void operator()(boost::iostreams::filtering_ostream &os) const;
    void operator()(boost::iostreams::filtering_istream &is) const;
    const std::string getFileNameExtension() const {return ".gz";}
    const std::string getName() const {return "gzip";}
};
#endif

#ifdef HAVE_BZIP2
template<>
class CompressionFilter<boost::iostreams::bzip2_params> : public Compression
{
public:
    CompressionFilter(const boost::iostreams::bzip2_params &params) : params_(params){}
private:
    const boost::iostreams::bzip2_params &params_;
    void operator()(boost::iostreams::filtering_ostream &os) const;
    void operator()(boost::iostreams::filtering_istream &is) const;
    const std::string getFileNameExtension() const {return ".bz2";}
    const std::string getName() const {return "bzip2";}
};
#endif

/**
 * @brief Compile-time resolution of boost::iostreams::*_params into
 * a corresponding compression filter
 */
template <typename ParamsT>
CompressionFilter<ParamsT> makeCompressionFilter(const ParamsT params){
    return CompressionFilter<ParamsT>(params);
}

/**
 ** Factory to get default instances of different types of supported compressions at runtime.
 **
 ** Usage:
 **     const Compression &compress;
 **     if (CompressionFactory::isSupported("gzip"))
 **         compress = CompressionFactory::get("gzip");
 **/
class CompressionFactory
{
public:
    static const std::vector<std::string> getCompressionList();
    static bool isSupported(const std::string &compression);
    static const Compression &get(const std::string &compression)
            throw(UnsupportedCompressionException);

    template<typename ParamsT>
    static const CompressionFilter<ParamsT> create(const ParamsT &params);

    static const Compression &none();
    typedef std::map<std::string, const Compression *> CompressionMap;
private:
    static const CompressionMap &getCompressionMap();
};

} // namespace common
} // namespace casava

#endif // #ifndef CASAVA_COMMON_COMPRESSION_HH
