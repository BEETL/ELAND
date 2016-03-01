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
 ** \file ELAND_options.cpp
 **
 ** \brief Command line options for ELAND.
 **
 ** \author Mauricio Varea
 **/
#include <sstream>
#include <iterator>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include "eland_ms/ELAND_options_ms.hh"
#include "common/Exceptions.hh"

namespace casava
{
namespace eland_ms
{

    ElandOptions::ElandOptions()
      : dataFormat_("bcl")
      , repeatFile_()
      , ungapped_(false)
      , singleseed_(false)
      , debug_(false)
      , sensitive_(false)
      , useBases_()
      , lane_(0)  // no default
      , read_(0)  // no default
      , inputDirectory_(".")
      , oligoLength_(0)
    {
      msg.push_back("[=N0[,N1,N2]]\n");
      msg[0] += "Output multiple hits per read. ";
      msg[0] += "At most N0,N1,N2 exact, 1-mismatch, 2-mismatch hits per read.";

      msg.push_back("file or directory of files\n");
      msg[1] += "  file type deduced from first character of each file:\n";
      msg[1] += "  '>' - fasta format\n";
      msg[1] += "  '#' - single molecule array format\n";
      msg[1] += "  [AGCTNagctn] - raw sequence format\n";

      msg.push_back("directory of genome files\n");
      msg[2] += "  preprocessed to 2-bits-per-base format using squashGenome\n";

      msg.push_back("name of output file\n");
      msg[3] += "  if name ends in '.vmf', use verbose match format,\n";
      msg[3] += "  else use format required by assembly module\n";

      msg.push_back("list of tiles to process\n");
      msg[4] += "  (only used when reading qseq files)\n";

      namedOptions_.add_options()
          ("multi", po::value< std::string >(&multi_)->implicit_value("10"),
                    msg[0].c_str())
          ("repeat-file", po::value< fs::path >(&repeatFile_),
                    "if given, points to a file containing the list of repeats to exclude (must be ASCII and in alphabetical order)")
          ("ungapped", po::value< bool >(&ungapped_)->zero_tokens(),
                    "output ungapped alignments instead of gapped")
          ("singleseed", po::value< bool >(&singleseed_)->zero_tokens(),
                    "do not use multiple seeds per read")
          ("debug", po::value< bool >(&debug_)->zero_tokens(),
                    "write the multi files")
          ("sensitive", po::value< bool >(&sensitive_)->zero_tokens(),
                    "increase sensitivity")
          ("lane", po::value< unsigned int >(&lane_),
                    "lane number (only used when reading qseq or bcl files)")
          ("read", po::value< unsigned int >(&read_),
                    "read number (only used when reading qseq or bcl files)")
          ("tiles", po::value< std::string >(&tilesString_)->default_value(""),
                    "list of tiles (only used when reading qseq or bcl files)")
          ("sample", po::value< std::string >(&sample_)->default_value("Sample"),
                    "sample name (for use with --data-format=fastq")
          ("barcode", po::value< std::string >(&barcode_)->default_value("empty"),
                    "barcode (for use with --data-format=fastq")
          ("cluster-sets", po::value< std::string >(&clusterSetsString_)->default_value(""),
                    "list of decimal cluster set numbers (for use with --data-format=fastq")
          ("instrument-name", po::value< std::string >(&instrumentName_)->default_value("unknown-instrument"),
                    "instrument name to use for the identification of the sequences (bcl input only)")
          ("run-number", po::value< unsigned int >(&runNumber_)->default_value(0),
                    "run-number to use for the identification of the sequences (bcl input only)")
          ("data-format", po::value< std::string >(&dataFormat_)->default_value("bcl"),
                    "format of the input data (bcl, qseq, fastq, fasta)")
          ("oligo-file", po::value< fs::path >(&oligoFile_),
                    "file containing the data (only for fastq and fasta format)")
          ("base-calls-dir", po::value< fs::path >(&inputDirectory_)->default_value("."),
                    "path to Base Calls directory")
          ("filter-directory", po::value< fs::path >(&filterDirectory_),
                    "directory containing the filter files, if different from the base calls directory (only for bcl input)")
          ("positions-directory", po::value< fs::path >(&positionsDirectory_),
                    "directory containing the positions files, if different from the parent of the base calls directory (only for bcl input)")
          ("positions-format", po::value< std::string >(&positionsFormatString_)->default_value("locs"),
                    "format of the position files, either 'locs', 'clocs' or 'txt' (only for bcl input)")
          ("output-file", po::value< fs::path >(&outputFile_),
                    "full path to the output file")
          ("tmp-file-prefix", po::value< fs::path >(&tmpFilePrefix_),
                    "path (including the file name) to form the temporary file paths. If unspecified, eland will create unique files in system temporary folder.")
          ("genome-directory", po::value< fs::path >(&genomeDirectory_),
                    "directory containing the squashed reference files")
          ("cycles", po::value< std::string >(&cycleString_),
                    "list of cycles to align (only for bcl input)")
          ("qseq-mask", po::value< std::string >(&useBases_),
                    "conversion mask - 'Y' (or 'y'), 'N' (or 'n') for 'use' or 'discard' respectively (only used when reading qseq files)")
          ("oligo-length", po::value< unsigned int >(&oligoLength_), "Seed length. Valid range is [8-32]")
          ;

      //argsH_.resize(2);
      //unnamedOptions_.add_options()
          //("arg1", po::value< std::string >(&argsH_[0]),
          //          "First positional argument")
          //("arg2", po::value< std::string >(&argsH_[1]),
          //          "Second positional argument")
          //("reminder", po::value< std::vector<std::string> >(&argsT_),
          //          "Rest of positional arguments");

      //positionalOptions_.add("arg1",1);
      //positionalOptions_.add("arg2",1);
      //positionalOptions_.add("reminder",-1);
    }

    std::string ElandOptions::usagePrefix() const
    {
        std::string usage = "Usage: eland_ms_" + boost::lexical_cast<std::string>(oligoLength_);
        usage += " oligoFile genomeDirectory outputFile[.vmf] [options]\n";
        usage += "   or: eland_ms_" + boost::lexical_cast<std::string>(oligoLength_);
        usage += " --qseq-source genomeDirectory outputFile[.vmf] [options]";
        usage += " tile1 [tile2 [.. tileN]]\n\n";

        usage += "oligoFile - " + msg[1] + "\n";
        usage += "genomeDirectory - " + msg[2] + "\n";
        usage += "outputFile - " + msg[3] + "\n";
        usage += "tile{1..N} - " + msg[4];

        return usage;
    }

    void ElandOptions::postProcess(po::variables_map &vm)
    {
        // do not process exceptions if "--help" was given
        if (vm.count("help"))  return;

        using casava::common::InvalidOptionException;
        if (! vm.count("multi") )
        {
            BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Missing switch '--multi'. For the moment, this is a mandatory switch ***\n"));
        }

        size_t p=0,n;
        while ( (n=multi_.find(',',p)) != std::string::npos )
        {
            maxNumMatches_.push_back( boost::lexical_cast<int>(multi_.substr(p,n-p)) );
            p=n+1;
        }
        maxNumMatches_.push_back( boost::lexical_cast<int>(multi_.substr(p,n-p)) );

        if (maxNumMatches_.size() == 1) {
        	maxNumMatches_.resize( 3, maxNumMatches_[0] );
        } else if (maxNumMatches_.size() != 3) {
            BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Problem parsing '--multi' CLI argument. Please provide either 0, 1, or 3 values ***\n"));
        }

        if (32 < oligoLength_ || oligoLength_ < 8) {
          BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** Problem parsing '--oligo-length' CLI argument. Please provide value in range [8-32] ***\n"));
        }

        // positional arguments interpretation depends on qseq-mode of operation
        //std::vector<std::string> tmp;
        //if (vm.count("reminder"))
        //{
        //	tmp = vm["reminder"].as< std::vector<std::string> >();
        //}
        if ("qseq" == dataFormat_ || "bcl" == dataFormat_)
        {
            std::istringstream is(tilesString_);
            std::copy(std::istream_iterator<unsigned int>(is), std::istream_iterator<unsigned int>(), std::back_inserter(tiles_));
            if (!is.eof())
            {
                BOOST_THROW_EXCEPTION(InvalidOptionException(
                              (boost::format("\n   *** failed to parse the list of tiles: '%s' ***\n") % tilesString_).str()));
            }
            if (tiles_.empty())
            {
                BOOST_THROW_EXCEPTION(InvalidOptionException(
                              (boost::format("\n   *** at least one tile must be provided for %s input format ***\n") % dataFormat_).str()));
            }
        }
        if ("fastq" == dataFormat_)
        {
            boost::tokenizer<boost::char_separator<char> > tknzr(clusterSetsString_, boost::char_separator<char>(" \t,"));
            std::transform(tknzr.begin(), tknzr.end(), std::back_inserter(clusterSets_),
                           static_cast<unsigned int (*) (const std::string&)>(&boost::lexical_cast<unsigned int>));

            if ( clusterSets_.empty() ) {
                BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** cluster-sets not valid: please provide list of cluster set numbers ***\n"));
            }
        }
        if ("bcl" == dataFormat_)
        {
            std::map<std::string, boost::format> positionsFormatMap;
            positionsFormatMap["txt"] = boost::format("s_%u_%04u_pos.txt");
            positionsFormatMap["locs"] = boost::format("s_%u_%04u.locs");
            positionsFormatMap["clocs"] = boost::format("s_%u_%04u.clocs");
            if (positionsFormatMap.end() == positionsFormatMap.find(positionsFormatString_))
            {
                BOOST_THROW_EXCEPTION(InvalidOptionException(
                              (boost::format("\n   *** invalid positions format: %s: supported formats are 'txt', 'locs' and 'clocs' ***\n") % positionsFormatString_).str()));

            }
            positionsFormat_ = positionsFormatMap[positionsFormatString_];
            std::istringstream is(cycleString_);
            std::copy(std::istream_iterator<unsigned int>(is), std::istream_iterator<unsigned int>(), std::back_inserter(cycles_));
            if (cycles_.empty())
            {
                BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** no cycles have been specified ***\n"));
            }
            if (!is.eof())
            {
                BOOST_THROW_EXCEPTION(InvalidOptionException((boost::format("\n   *** invalid cycles list: %s ***\n") % cycleString_).str()));
            }
        }
        if ("qseq" == dataFormat_)
        {
            //outputFile_ = vm["arg1"].as< std::string >();
            //for (std::vector<std::string>::iterator it=tmp.begin(); it != tmp.end(); ++it)
            //{
            //    if ( it->find_first_not_of("0123456789") == std::string::npos )
            //    {
            //        tiles_.push_back( boost::lexical_cast<unsigned int>(*it) );
            //    } else {
            //        BOOST_THROW_EXCEPTION(InvalidOptionException(
            //                  (boost::format("\n   *** '%s' is not a valid tile number ***\n") % it->c_str() ).str()
            //              ));
            //    }
            //}

            // validate the input directory
            if ( inputDirectory_.empty() )
            {
                BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** the input directory can't be empty ***\n"));
            } else if (! fs::is_directory(inputDirectory_) ) {
                BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** the input directory must exist ***\n"));
            }

            // validate the lane
            if ( 1 > lane_ || lane_ > 8 )
            {
                BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** lane not valid: please provide an integer in the range '1 <= n <= 8' ***\n"));
            }

            //validate use-bases
            size_t p = 0;
            if ( !useBases_.length() )
            {
                // should default to [Y*] but, given we don't know the lenght beforehand,
                // we'll make it compulsory for the moment
                BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** missing --qseq-mask option ***\n"));
            } else if ( (p = useBases_.find_first_not_of("YyNn0123456789")) != std::string::npos ) {
                BOOST_THROW_EXCEPTION(InvalidOptionException( (boost::format("\n   *** '%c' is not a valid char in --qseq-mask ***\n") % useBases_[p] ).str() ));
            }

            // validate the repeat file
            if (!repeatFile_.empty() && !fs::exists(repeatFile_) )
            {
                BOOST_THROW_EXCEPTION(InvalidOptionException("\n   *** the --repeat-file provided does not exist ***\n"));
            }

        } else {
            //oligoFile_ = vm["arg1"].as< std::string >();
            //genomeDirectory_ = vm["arg2"].as< std::string >();
            //outputFile_ = tmp[0];
        }
    }

} // eland_ms
} // casava
