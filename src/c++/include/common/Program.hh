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
 ** \file Program.hh
 **
 ** Declaration of the skeleton of all c++ programs.
 **
 ** \author Come Raczy
 **/

#ifndef CASAVA_COMMON_PROGRAM_HH
#define CASAVA_COMMON_PROGRAM_HH

#include <string>
#include <iostream>
#include <cstdlib>
#include <boost/program_options.hpp>

#include "common/Exceptions.hh"

namespace casava
{
namespace common
{

namespace po = boost::program_options;

/**
 ** Encapsulation of the pocessing of the command line options.
 **
 ** TODO: add config file and environment options
 **/
class Options
{
public:
    enum Action
    {
        RUN, HELP, ABORT
    };
    Options();
    virtual ~Options()
    {
    }
    Action parse(int argc, char *argv[]);
    std::string usage() const;
protected:
    po::options_description namedOptions_;
    po::options_description unnamedOptions_;
    po::positional_options_description positionalOptions_;
private:
    virtual std::string usagePrefix() const = 0;
    virtual std::string usageSuffix() const
    {
        return "";
    }
    virtual void postProcess(po::variables_map &)
    {
    }
    virtual bool tolerateCLIsyntax(po::invalid_command_line_syntax &) const
    {
        return false;
    }
    virtual bool correctCLIsyntax(int /*argc*/, char** &/*argv*/)
    {
        return false;
    }
};

/**
 ** Unified behavior of all programs.
 **/
template<class O>
void run(void(*callback)(const O &), int argc, char *argv[])
{
    try
    {
        O options;
        const typename O::Action action = options.parse(argc, argv);
        if (O::RUN == action)
        {
            callback(options);
        }
        else if (O::HELP == action)
        {
            std::cout << options.usage() << std::endl;
        }
        else
        {
            std::clog << options.usage() << std::endl;
            exit(1);
        }
    }
    catch (const casava::common::ExceptionData &exception)
    {
        std::clog << "Error: " << exception.getContext() << ": " << exception.getMessage() << std::endl;
        exit (1);
    }
    catch (const std::runtime_error &e)
    {
        std::clog << "runtime error: " << e.what() << std::endl;
        exit(2);
    }
    catch (const std::logic_error &e)
    {
        std::clog << "logic error: " << e.what() << std::endl;
        exit(3);
    }
}

} // namespace common
} // namespace casava

#endif // #ifndef CASAVA_COMMON_PROGRAM_HH
