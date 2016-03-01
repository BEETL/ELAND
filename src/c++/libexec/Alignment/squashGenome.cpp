/**
 ** Copyright (c) 2003-2006 Solexa Limited. All rights reserved.
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
 ** \file squashGenome.cpp
 **
 ** \brief Application wrapper for squash genome apis
 **
 ** \author Roman Petrovski
 **/

#include <dirent.h>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "alignment/SquashGenome.hh"
#include "common/Program.hh"

namespace casava
{
namespace alignment
{

namespace po = boost::program_options;
namespace fs = boost::filesystem;

class SquashGenomeOptions : public casava::common::Options
{
public:
    SquashGenomeOptions();
    bool validateChromNames;
    bool allowManyContigs;
    std::string chromNameSource;
    std::vector<boost::filesystem::path> filesToSquash;
    boost::filesystem::path fileToUnsquash;
    boost::filesystem::path squashDirectory;
    boost::filesystem::path squashedFileOrDirectory;
    int logLevel;
private:
    std::string usagePrefix() const {
        return "Usage: \n"
        "squashGenome [options] <fileToUnsquash>\n"
        "  - unsquash file to standard output\n"
        "squashGenome [options] <squashDirectory>\n"
        "  - print xml containing contig size information of squashed files to standard output\n"
        "squashGenome [options] <targetDirectory> <fileToSquash1> <fileToSquash2> ...\n"
        "  - squash files and place the results in <targetDirectory>";
    }
    void postProcess(boost::program_options::variables_map &vm);
};


SquashGenomeOptions::SquashGenomeOptions()
    : validateChromNames(true)
    , allowManyContigs(false)
    , logLevel(1)
{
    // short namespaces
    namedOptions_.add_options()
        ("allow-many-contigs", po::value< bool >(&allowManyContigs),
                               "Will fail on .fa files containing multiple contigs if disabled. "
                               "On by default if --chrom-name-source is contigName. (Squash only)")
        ("validate-names"    , po::value< bool >(&validateChromNames),
                               ("Will fail on .fa files with contig names containing the following "
                               "characters: " + getContigNameForbiddenCharacters() +
                               "\nOn by default if --chrom-name-source is contigName. (Squash only)").c_str())
        ("chrom-name-source" , po::value< std::string >(&chromNameSource)->default_value("fileName"),
                               "Valid options are: contigName or fileName. Required for validations. (Squash only)")
        ("verbose-level,v" ,   po::value< int >(&logLevel)->default_value(1),
                               "Valid options are: 0 - no logging, 1 - user-level information and critical messages. 2 and above - debug loggging")
                               ;

    unnamedOptions_.add_options()
        ("squashedFileOrDirectory", po::value< fs::path >(&squashedFileOrDirectory),
                  "Path to a file to unsquash or a directory containing squashed files")
        ("filesToSquash", po::value< std::vector<fs::path> >(&filesToSquash),
                  "List of .fa files to squash")
                ;

    positionalOptions_.add("squashedFileOrDirectory",1);
    positionalOptions_.add("filesToSquash",-1);
}

void SquashGenomeOptions::postProcess(boost::program_options::variables_map &vm)
{
    // do not process exceptions if "--help" was given
    if (vm.count("help"))  return;

    using casava::common::InvalidOptionException;

    if (squashedFileOrDirectory.empty())
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** no parameters given ***\n"));
    }

    if ( filesToSquash.empty() && !fs::is_directory(squashedFileOrDirectory) )
    {
        fileToUnsquash = squashedFileOrDirectory;
    }
    else
    {
        squashDirectory = squashedFileOrDirectory;
    }

    if ( "fileName" != chromNameSource && "contigName" !=  chromNameSource)
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** the --chrom-name-source must be fileName or contigName ***\n"));
    }

    if (!vm.count("validate-names") && "fileName" == chromNameSource)
    {
        validateChromNames = false;
    }

    if (!vm.count("allow-many-contigs") && "contigName" == chromNameSource)
    {
        allowManyContigs = true;
    }
}

void squashGenome(const SquashGenomeOptions &options)
{
    if (!options.fileToUnsquash.empty())
    {
        unsquash(options.fileToUnsquash.string().c_str(), options.logLevel);
    }
    else
    {
        if (options.filesToSquash.empty())
        {
            std::cerr << "INFO: Trying to open directory " << options.squashDirectory << " ...\n";
            DIR* pDir = opendir(options.squashDirectory.string().c_str());
            if (pDir)
            {
                std::cerr << "INFO: ... success, will output file sizes to XML\n";
                outputSizesToXML(options.squashDirectory.string().c_str(), pDir, options.logLevel);
                closedir( pDir );
            }
            else
            {
                std::cerr << "ERROR: ... could not open directory: " << options.squashDirectory << "\n";
                exit(2);
            }
        }
        else
        {
            BOOST_FOREACH(const fs::path &fileToSquash, options.filesToSquash)
            {
                squash(options.squashDirectory.string().c_str(), fileToSquash.string().c_str(),
                       options.validateChromNames, options.allowManyContigs, options.logLevel);
            }
        }
    }
}

}
}

int main(int argc, char *argv[])
{
    casava::common::run(casava::alignment::squashGenome, argc, argv);
}
