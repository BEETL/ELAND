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
 ** \file Program.cpp
 **
 ** Implementation of the skeleton of all c++ programs.
 **
 ** \author Come Raczy
 **/

#include "common/Program.hh"

namespace casava
{
namespace common
{

Options::Options() :
    namedOptions_("Command line options")
{
    namedOptions_.add_options()("help,h", "produce help message and exit");
}

Options::Action Options::parse(int argc, char *argv[])
{
    try
    {
        po::options_description allOptions("Allowed options");
        allOptions.add(namedOptions_).add(unnamedOptions_);
        po::variables_map vm;
        po::store(
                po::command_line_parser(argc, argv). options(allOptions).positional(
                        positionalOptions_).run(), vm);
        po::notify(vm);
        postProcess(vm);
        if (vm.count("help"))
        {
            return HELP;
        }
        else
        {
            return RUN;
        }
    }
    catch (po::invalid_command_line_syntax& e)
    {
        if (this->tolerateCLIsyntax(e) && this->correctCLIsyntax(argc,argv))
            return parse(argc, argv);
//        std::clog << "Failed to parse the options: "
//                  << po::invalid_command_line_syntax::error_message(e.kind()) << std::endl;
        return ABORT;
    }
    catch (const std::exception &e)
    {
        std::clog << "Failed to parse the options: " << e.what() << std::endl;
        return ABORT;
    }
}

std::string Options::usage() const
{
    std::ostringstream os;
    os << this->usagePrefix() << std::endl << std::endl;
    os << namedOptions_ << std::endl;
    return os.str();
}

} // namespace common
} // namespace casava
