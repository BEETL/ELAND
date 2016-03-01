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
 ** \file eland_ms/ELAND_main_ms.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ALIGNMENT_ELAND_MAIN_HH
#define CASAVA_ALIGNMENT_ELAND_MAIN_HH

#include <dirent.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include "alignment/OligoSourceBcl.hh"
#include "alignment/OligoSourceQseq.hh"
#include "alignment/OligoSourceFastq.hh"

#include "MatchPositionTranslator.hh"
#include "SuffixScoreTable.hh"
#include "MatchTable.hh"
#include "RepeatTable.hh"

#include "ElandThread.hh"
#include "ElandDefines.hh"

namespace casava
{
namespace eland_ms
{

namespace fs = boost::filesystem;
using casava::common::CasavaException;
namespace ca = casava::alignment;


template <int OLIGO_LEN>
class ELAND
{
public:

private:
    const unsigned int oligoLength;
    const std::list<fs::path> qseq_file_list;
    const std::string genome_dir;
    OligoSource* pOligos;
    short no_of_seeds;
    MatchTable* pResults;
    MatchTable* pResults_2;
    const bool do_ungapped;
    bool do_singleseed;
    const bool do_debug;
    const bool do_sensitive;
    Timer timer;

    OligoSource* getOligoSource(const std::string &dataFormat,
                                   const std::string &machineName,
                                   const unsigned int runNumber,
                                   const unsigned int lane,
                                   const unsigned int read,
                                   const std::vector<unsigned int> &tiles,
                                   const std::string &sample,
                                   const std::string &barcode,
                                   const std::vector<unsigned int> &clusterSets,
                                   const fs::path &inputDirectory,
                                   const fs::path &filterDirectory,
                                   const fs::path &positionsDirectory,
                                   const std::string &useBases,
                                   const std::vector<unsigned int> &cycles,
                                   const fs::path &oligoFile,
                                   const boost::format &positionsFileNameFormat) const
{
    using boost::format;
    using namespace std;
    if (dataFormat == "bcl")
    {
        const fs::path laneDirectory = inputDirectory / (format("L%03u") % lane).str();
        std::vector<fs::path> bclDirectoryList;
        BOOST_FOREACH(unsigned int cycle, cycles)
        {
            bclDirectoryList.push_back(laneDirectory / (format("C%u.1") % cycle).str());
        }
        const std::vector<unsigned int> barcodeCycles;
        std::vector<fs::path> barcodeDirectoryList;
        BOOST_FOREACH(unsigned int cycle, barcodeCycles)
        {
            barcodeDirectoryList.push_back(laneDirectory / (format("C%u.1") % cycle).str());
        }
        return new ca::OligoSourceBcl(bclDirectoryList,
                                  barcodeDirectoryList,
                                  positionsDirectory,
                                  filterDirectory,
                                  positionsFileNameFormat,
                                  machineName,
                                  runNumber,
                                  lane,
                                  tiles,
                                  read);
    }
    else if (dataFormat == "qseq")
    {
        std::list<fs::path> qseqFileList;
        BOOST_FOREACH(unsigned int tile, tiles)
        {
            qseqFileList.push_back(inputDirectory /
                                   (format("s_%u_%u_%04u_qseq.txt") % lane % read % tile).str());
        }
        return new ca::OligoSourceQseq(qseqFileList, useBases);
    }
    else if (dataFormat == "fastq")
    {
        return new ca::OligoSourceFastq(inputDirectory, sample, barcode, lane, read, clusterSets, useBases);
    }
    // fallback to the auto-detection
    return OligoSourceFile::getOligoSourceFile(oligoFile.string().c_str());
}
public:

ELAND(const fs::path &oligoFile,
             const fs::path &genomeDirectory,
             const fs::path &outputFile,
             const std::vector<unsigned int> &maxNumMatches,
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
    : oligoLength(OLIGO_LEN)
    , genome_dir(genomeDirectory.string())
    , pOligos(getOligoSource(dataFormat, instrumentName, runNumber, lane, read, tiles,
            sample, barcode, clusterSets,
            inputDirectory, filterDirectory, positionsDirectory, useBases, cycles,
            oligoFile, positionsFileNameFormat))
    , pResults(NULL)
    , do_ungapped(ungap)
    , do_singleseed(singleSeed)
    , do_debug(debug)
    , do_sensitive(sensitive)
{
    assert(oligoLength!=0);



  // run from the constructor, at the moment, to mimic legacy behaviour.
  // we should look into factoring this out
  presentation();

  cerr << "Will use "
       << oligoLength
       << " oligos for seeds."
       << endl << endl;

  cerr << "Will use at most " << MAX_HASH_BITS
       << " bits in hash table" << endl;
  cerr << "Can process at most " << maxNumOligos
       << " oligos per batch" << endl;

#ifdef ONE_ERROR_PER_OLIGO
  cerr << "Will find all exact and single error matches" << endl;
#else
  cerr << "Will find all matches having 2 errors or less" << endl;
#endif

#ifdef DONT_SEARCH_REVERSE_STRAND
  cerr << "WARNING: will search for oligos in forward strand only"<< endl
       << "To search both strands, comment out DONT_SEARCH_REVERSE_STRAND "
       << "and recompile." << endl;
#endif

  const char* pTest=pOligos->getNextOligo();
  const int read_length = pTest ? strlen(pTest) : 0;
  if (pTest==NULL)
  {
    cerr << "WARNING: there do not appear to be any sequences in file "
	 << oligoFile.string() << endl;
  } // ~if
  else if (strlen(pTest)<oligoLength)
  {
    cerr << "ERROR: first sequence of " << oligoFile.string() << " contains "
	 << strlen(pTest) << " bases, but "
	 << oligoLength << " bases are needed for the alignment."
	 << endl
     << "Please use a different value for --oligo-len."
	 << endl << endl;
    exit(-1);
  } // ~else if
  else if (strlen(pTest)>oligoLength)
  {
//    cerr << "WARNING: first sequence of " << oligoFile.string() << " contains "
//	 << strlen(pTest) << " bases, but only the first "
//	 << oligoLength << " bases will be used in the alignment."
//	 << endl
//	 << "To change this value, edit OLIGO_LENGTH and recompile."
//	 << endl << endl;
  } // ~else if

  // deriving the number of seeds
  no_of_seeds = 1;
  uint bases_not_covered = 0;
  if (NULL != pTest)
  {
    no_of_seeds = strlen(pTest) / oligoLength;

    if( (strlen(pTest)-oligoLength)< 4 && (do_singleseed == false) )
    {
        do_singleseed = true;
        cerr << "Not running in multiseed mode because reads are too short." << endl;
    }

    if( no_of_seeds > 4 )
    {
        no_of_seeds = 4;
    }

    bases_not_covered = strlen(pTest) - no_of_seeds * oligoLength;
  }

  cerr << no_of_seeds << " seeds of length " << oligoLength << " will be used." << endl;
  cerr << bases_not_covered << " bases will not be covered by any of the seeds." << endl;

  pOligos->rewind();

  cerr << "Will read oligos from file " << oligoFile.string() << endl;

  cerr << "Will perform " << (do_ungapped ? "un" : "") << "gapped alignment." << endl;

  if (do_singleseed)
      cerr << "Will use only one seed per read." << endl;

  // MULTI command line argument
  if (3 == maxNumMatches.size())
  {
      pResults = new MatchTableMulti(OLIGO_LEN,outputFile.string().c_str(),
                                     debug,
                                     maxNumMatches[0],
                                     maxNumMatches[1],
                                     maxNumMatches[2],
                                     tmpFilePrefix.empty() ? 0 : tmpFilePrefix.string().c_str());
      pResults->setSensitivity(do_sensitive);

      // build up a second match table for the second tier
      pResults_2 = new MatchTableMultiSquareSeed(OLIGO_LEN,"/dev/null",
                                                 false,
                                                 (maxNumMatches[0]*6),
                                                 (maxNumMatches[1]*6),
                                                 (maxNumMatches[2]*6),
                                                 tmpFilePrefix.empty() ? 0 : (tmpFilePrefix.string() + ".t2").c_str());
//      pResults_2->setNoOfSeeds(no_of_seeds);
      pResults_2->setSensitivity(do_sensitive);
      pResults_2->setNoOfSeeds(4);
      pResults_2->setReadLength(read_length);
  } else {
      BOOST_THROW_EXCEPTION( CasavaException(EINVAL,
            (boost::format("Cannot match Table with given --multi. Expected 3, got %d parameters.") % maxNumMatches.size()).str()
      ));
  }

  RepeatTable<OLIGO_LEN>* pRepeats(NULL);
  if (!repeatFile.empty())
      pRepeats = new RepeatTable<OLIGO_LEN>(repeatFile.string().c_str());

  if (pRepeats!=NULL)
  {
    cerr << "Scanning for repeats in list " << repeatFile.string() << ": " << timer << endl;
    pRepeats->checkOligos( *pOligos, *pResults );
    cerr << "Scanned repeats: " << timer << endl;
    delete pRepeats;
  } // ~if
}


~ELAND()
{
    delete pOligos;
    delete pResults;
}

void presentation()
{
  cout << endl
       << "------------------------------------------------------------"
       << endl;

  cout << "ELAND: Efficient Local Alignment of Nucleotide Data" << endl
       << "Copyright (c) 2003-2006 Solexa Limited. All rights reserved."
       << endl << "Author: Anthony J. Cox" << endl << endl;

  cout << "Publications incorporating data generated by the use of\n"
       << "this software or modified versions thereof should cite:\n"
       << "Anthony J. Cox.\n"
       << "Ultra high throughput alignment of short sequence tags.\n"
       << "In preparation." << endl << endl;

  cout << "------------------------------------------------------------"
       << endl << endl;
}


void run()
{
  cerr << "Starting run! Time now: " << timer.timeNow();

  SuffixScoreTable scoreTable(ElandConstants<OLIGO_LEN>::fragLengthA,
                              ElandConstants<OLIGO_LEN>::fragLengthB,
                              ElandConstants<OLIGO_LEN>::fragLengthC,
                              ElandConstants<OLIGO_LEN>::fragLengthD);
  HashTableDataStore<ElandConstants<OLIGO_LEN>::useSplitPrefix> htds1, htds2;
  // second tier
  HashTableDataStore<ElandConstants<OLIGO_LEN>::useSplitPrefix> htds1_2, htds2_2;

  const char suffixName[] = ".2bpb";
  cerr << "Trying to open directory " << genome_dir  << " ..." << endl;

  DIR* pDir;
  if ( ! ( pDir = opendir( genome_dir.c_str() ) ) )
  {
    cerr << " ... failed!" << endl;
    exit(-1);
  } // ~if

  char dirBuffer[80];
  dirent* dirEntry;

  // chromosome names are indexed starting at 1
  // (because in OligoInfo a chromosome num of 0 indicates no match found)
  vector<string> chromNames(1); // %%%%% (1);
  vector<MatchPosition> blockStarts(1);

  // variables for the second tier
  vector<string> chromNames_2(1); // %%%%% (1);
  vector<MatchPosition> blockStarts_2(1);


  const string directoryName((genome_dir)+'/');
  //  string fullChromName;

  // get chromosome names from directory then close it
  while( (dirEntry = readdir(pDir)) )
  {
    strcpy(dirBuffer, dirEntry->d_name);
    // keep only the files with suffix matching suffixName
    char* const pSuffixStart = dirBuffer + strlen(dirBuffer) - strlen(suffixName);
    if ( strlen(dirBuffer) >= strlen(suffixName) && 0 == strcmp(pSuffixStart, suffixName) )
    {
      *pSuffixStart='\0';
      chromNames.push_back((string)dirBuffer);
      chromNames_2.push_back((string)dirBuffer);
    } // ~if
  } // ~while
  closedir( pDir );

  // Add in fix to issue IMP-21 - sort chromosome names so as to always get
  // consistent results
  cerr << "Sorting chromosome names" << endl;
  sort (chromNames.begin(),chromNames.end());
  sort (chromNames_2.begin(),chromNames_2.end());


  // do pass 0
  {
    OligoHashTable<0, OLIGO_LEN> hashTable (oligoLength, htds1, htds2, scoreTable, *pResults);
    scanAll<0, OLIGO_LEN>( pOligos, directoryName, chromNames, blockStarts, hashTable, timer,true );
  } // ~scope of hashTable

#ifndef ONE_ERROR_PER_OLIGO
  // do pass 1
  {
    OligoHashTable<1, OLIGO_LEN> hashTable (oligoLength, htds1, htds2, scoreTable, *pResults);
    scanAll<1, OLIGO_LEN>( pOligos, directoryName, chromNames, blockStarts, hashTable, timer,true );
  } // ~scope of hashTable

  // do pass 2
  {
    OligoHashTable<2, OLIGO_LEN> hashTable (oligoLength, htds1, htds2, scoreTable, *pResults);
    scanAll<2, OLIGO_LEN>( pOligos, directoryName, chromNames, blockStarts, hashTable, timer,true );
  } // ~scope of hashTable
#endif

  // clear some space - pResults->print may need it for MatchTableMulti
  htds1.clear();
  htds2.clear();

  MatchPositionTranslator getMatchPos( chromNames, blockStarts, directoryName );


  if( !do_singleseed )
  {

      cerr << "Looking for unmapped reads... ";
      vector<bool> unmappedReads;
      if( pResults->getUnmappedReads( unmappedReads ) == true )
      {
          cerr << "done." << endl;
      }
      else
      {
          cerr << "failed." << endl;
      }

      cerr << "Setting oligo mask...";
      pOligos->setMask( unmappedReads );

      // rewind pOligos again
      pOligos->rewind();


      // do the multiseed stage
      cerr << "Performing multi-seed for reads not matched so far..." << endl;

      // do pass 0
      {
          OligoHashTable<0, OLIGO_LEN> hashTable(oligoLength, htds1_2, htds2_2, scoreTable, *pResults_2);
          scanAll<0, OLIGO_LEN>( pOligos, directoryName, chromNames_2, blockStarts_2, hashTable, timer,false );
      } // ~scope of hashTable

#ifndef ONE_ERROR_PER_OLIGO
      // do pass 1
      {
          OligoHashTable<1, OLIGO_LEN> hashTable(oligoLength, htds1_2, htds2_2, scoreTable, *pResults_2);
          scanAll<1, OLIGO_LEN>( pOligos, directoryName, chromNames_2, blockStarts_2, hashTable, timer,false );
      } // ~scope of hashTable

      // do pass 2
      {
          OligoHashTable<2, OLIGO_LEN> hashTable(oligoLength, htds1_2, htds2_2, scoreTable, *pResults_2);
          scanAll<2, OLIGO_LEN>( pOligos, directoryName, chromNames_2, blockStarts_2, hashTable, timer,false );
      } // ~scope of hashTable
#endif

      htds1_2.clear();
      htds2_2.clear();


      // reset pOligos, otherwise we only print a subset of all the reads
      pOligos->unSetMask();

      // you have to ensure that pResults_2 is of type MatchTableMulti
      cerr << "Merging results..." << endl;
      if( pResults->mergeTable( pResults_2,getMatchPos ) == false )
      {
          cerr << "Error retrieving match information from the second run, will use information only from singleseed run." << endl;
      }
      // get rid of pResults_2 do save memory
      delete pResults_2;
      cerr << "done." << endl;

  } //~if( !do_singleseed )


  cerr << "Outputting results: " << timer << endl;


  pResults->setNoOfSeeds(no_of_seeds);
  pResults->printSquash
    ( *pOligos , getMatchPos, chromNames, blockStarts, scoreTable, oligoLength,directoryName,!do_ungapped );


  cerr << "... done " << timer << endl;
  cerr << "Run complete! Time now: " << timer.timeNow();
} // ~run

};


} // eland_ms
} // casava

#endif /* CASAVA_ALIGNMENT_ELAND_MAIN_HH */
