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
 ** \file Exceptions.cpp
 **
 ** \brief Implementation of the common exception mechanism.
 **
 ** \author Come Raczy
 **/

#include <cstring>
#include <errno.h>
#include <boost/date_time.hpp>

#include "common/Exceptions.hh"

namespace casava
{
namespace common
{

ExceptionData::ExceptionData(int errorNumber, const std::string &message) : boost::exception(),
            errorNumber_(errorNumber), message_(message)
{
}

std::string ExceptionData::getContext() const
{
    const std::string now = boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time());
//    return
//        std::string(boost::get_error_info<boost::throw_file>(*this) ?
//                    *boost::get_error_info<boost::throw_file>(*this) : "unknown file") + ": " +
//        std::string(boost::get_error_info<boost::throw_function>(*this) ?
//                    *boost::get_error_info<boost::throw_function>(*this) : "unknown function") + ": " +
//        (boost::get_error_info<boost::throw_line>(*this) ?
//         "line " + boost::lexical_cast<std::string>(*boost::get_error_info<boost::throw_line>(*this)) : "unknown line") + ": " +
//        now + ": " + std::string(strerror(errorNumber_));
    return now + ": " + std::string(strerror(errorNumber_)) + ": " + boost::diagnostic_information(*this);
}

IoException::IoException(int errorNumber, const std::string &message)
    : std::ios_base::failure(message)
    , ExceptionData(errorNumber, message)
{
}

UnsupportedVersionException::UnsupportedVersionException(const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

InvalidParameterException::InvalidParameterException(const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

InvalidOptionException::InvalidOptionException(const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

PreConditionException::PreConditionException(const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

PostConditionException::PostConditionException(const std::string &message)
    : std::logic_error(message)
    , ExceptionData(EINVAL, message)
{
}

} // namespace common
} // namespace casava
