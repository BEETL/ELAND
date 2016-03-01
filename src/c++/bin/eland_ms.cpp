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
 ** \file ELAND_ms.cpp
 **
 ** \brief Stands for Efficient Local Alignment of Nucleotide Data
 **
 ** Matches oligos to genome data allowing for substitution errors and
 ** for ambiguity codes
 **
 ** \author Mauricio Varea
 **/
#include "eland_ms/ELAND_options_ms.hh"
#include "eland_ms/ELAND_main_ms.hh"
#include <boost/format.hpp>

template <int MAX_OLIGO_LEN>
void run_eland(unsigned int len,
      const fs::path &oligoFile,
      const fs::path &genomeDirectory,
      const fs::path &outputFile,
      const std::vector<unsigned int> &multi,
      const fs::path &repeatFile,
      const bool &singleSeed,
      const bool &debug,
      const bool &ungap,
      const bool &sensitive,
      const std::string &dataFormat,
      const std::string &useBases,
      const std::vector<unsigned int> &cycles,
      const fs::path &inputDirectory,
      const fs::path &filterDirectory,
      const fs::path &positionsDirectory,
      const std::string &instrumentName,
      const unsigned int runNumber,
      const unsigned int &lane,
      const unsigned int &read,
      const fs::path &tmpFilePrefix,
      const std::vector<unsigned int> &tiles,
      const std::string &sample,
      const std::string &barcode,
      const std::vector<unsigned int> &clusterSets,
      const boost::format &positionsFileNameFormat)
{
  if (MAX_OLIGO_LEN == len){
    casava::eland_ms::ELAND<MAX_OLIGO_LEN> eland(
        oligoFile,
        genomeDirectory,
        outputFile,
        multi,
        repeatFile,
        singleSeed,
        debug,
        ungap,
        sensitive,
        dataFormat,
        useBases,
        cycles,
        inputDirectory,
        filterDirectory,
        positionsDirectory,
        instrumentName,
        runNumber,
        lane,
        read,
        tmpFilePrefix,
        tiles,
        sample,
        barcode,
        clusterSets,
        positionsFileNameFormat);
    eland.run();
  } else {
    run_eland<MAX_OLIGO_LEN-1>(len, oligoFile,
        genomeDirectory,
        outputFile,
        multi,
        repeatFile,
        singleSeed,
        debug,
        ungap,
        sensitive,
        dataFormat,
        useBases,
        cycles,
        inputDirectory,
        filterDirectory,
        positionsDirectory,
        instrumentName,
        runNumber,
        lane,
        read,
        tmpFilePrefix,
        tiles,
        sample,
        barcode,
        clusterSets,
        positionsFileNameFormat);
  }
}


template <>
void run_eland<7>(unsigned int len,
      const fs::path &/*oligoFile*/,
      const fs::path &/*genomeDirectory*/,
      const fs::path &/*outputFile*/,
      const std::vector<unsigned int> &/*multi*/,
      const fs::path &/*repeatFile*/,
      const bool &/*singleSeed*/,
      const bool &/*debug*/,
      const bool &/*ungap*/,
      const bool &/*sensitive*/,
      const std::string &/*dataFormat*/,
      const std::string &/*useBases*/,
      const std::vector<unsigned int> &/*cycles*/,
      const fs::path &/*inputDirectory*/,
      const fs::path &/*filterDirectory*/,
      const fs::path &/*positionsDirectory*/,
      const std::string &/*instrumentName*/,
      const unsigned int /*runNumber*/,
      const unsigned int &/*lane*/,
      const unsigned int &/*read*/,
      const fs::path &/*tmpFilePrefix*/,
      const std::vector<unsigned int> &/*tiles*/,
      const std::string &/*sample*/,
      const std::string &/*barcode*/,
      const std::vector<unsigned int> &/*clusterSets*/,
      const boost::format &/*positionsFileNameFormat*/)
{
  BOOST_THROW_EXCEPTION(cc::InvalidParameterException(
        (boost::format("Eland oligo length %u not supported") % len).str()));
}


void eland_ms(const casava::eland_ms::ElandOptions &options)
{
  run_eland<32>(options.oligoLength_,
      options.oligoFile_,
      options.genomeDirectory_,
      options.outputFile_,
      options.maxNumMatches_,
      options.repeatFile_,
      options.singleseed_,
      options.debug_,
      options.ungapped_,
      options.sensitive_,
      options.dataFormat_,
      options.useBases_,
      options.cycles_,
      options.inputDirectory_,
      options.filterDirectory_,
      options.positionsDirectory_,
      options.instrumentName_,
      options.runNumber_,
      options.lane_,
      options.read_,
      options.tmpFilePrefix_,
      options.tiles_,
      options.sample_,
      options.barcode_,
      options.clusterSets_,
      options.positionsFormat_);
}


int main(int argc, char *argv[])
{
    casava::common::run(eland_ms, argc, argv);
}

