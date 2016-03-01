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
 ** \file ELAND_options.hh
 **
 ** \brief Command line options for ELAND.
 **
 ** \author Mauricio Varea
 **/

#ifndef CASAVA_ALIGNMENT_ELAND_OPTIONS_HH
#define CASAVA_ALIGNMENT_ELAND_OPTIONS_HH

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "common/Program.hh"

namespace fs = boost::filesystem;
namespace po = boost::program_options;

namespace casava
{
namespace eland_ms
{

      class ElandOptions : public casava::common::Options
      {
      public:
          ElandOptions();
          std::string dataFormat_;
          fs::path oligoFile_;
          fs::path genomeDirectory_;
          fs::path outputFile_;
          std::vector<unsigned int> maxNumMatches_;
          fs::path repeatFile_;
          bool ungapped_;
          bool singleseed_;
          bool debug_;
          bool sensitive_;
          std::string useBases_;
          std::vector<unsigned int> cycles_;
          unsigned int lane_;
          unsigned int read_;
          fs::path inputDirectory_;
          fs::path filterDirectory_;
          fs::path positionsDirectory_;
          boost::format positionsFormat_;
          std::vector<unsigned int> tiles_;
          std::string sample_;
          std::string barcode_;
          std::vector<unsigned int> clusterSets_;
          std::string instrumentName_;
          unsigned int runNumber_;
          fs::path tmpFilePrefix_;
          unsigned int oligoLength_;
      private:
          std::string usagePrefix() const;
          void postProcess(po::variables_map &);
          std::vector<std::string> msg;
          //std::vector<std::string> argsH_;
          //std::vector<std::string> argsT_;
          std::string multi_;
          std::string cycleString_;
          std::string tilesString_;
          std::string clusterSetsString_;
          std::string positionsFormatString_;
      };

} // eland_ms
} // casava

#endif /* CASAVA_ALIGNMENT_ELAND_OPTIONS_HH */
