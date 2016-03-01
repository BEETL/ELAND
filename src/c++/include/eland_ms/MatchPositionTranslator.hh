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
 ** \file eland_ms/MatchPositionTranslator.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 ** Contains interface functionality from ELAND - specifically, stuff to do
 ** with scoring and storage of alignments
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_MATCH_POSITION_TRANSLATOR_H
#define CASAVA_ELAND_MS_MATCH_POSITION_TRANSLATOR_H

#include "ElandConstants.hh"
#include "ContigNameFinder.hh"

namespace casava
{
namespace eland_ms
{


// class MatchPositionTranslator
// This converts a match position into chromosome name + position
// and optionally a contig name
class MatchPositionTranslator
{
public:
  MatchPositionTranslator( const vector<string>& chromNames,
                           const vector<MatchPosition>& blockStarts,
                           const string& directoryName ) :
    chromNames_(chromNames)
  {
    int thisChrom(0);

    for (vector<MatchPosition>::const_iterator i(blockStarts.begin());
         i!=blockStarts.end()-1;
         ++i)
    {
      MatchPosition thisChromStartBlock(*i), nextChromStartBlock(*(i+1));
      BOOST_ASSERT((thisChromStartBlock & blockPositionMask)==0 && "thisChromStartBlock must not have any position bits set");
      BOOST_ASSERT((nextChromStartBlock & blockPositionMask)==0 && "nextChromStartBlock must not have any position bits set");
      BOOST_ASSERT(thisChromStartBlock < nextChromStartBlock && "blockStarts starts must be ordered");
      for ( uint j = thisChromStartBlock; j != nextChromStartBlock; j+=blockSize)
      {
          subtractTable_[j>>blockShift]=*i;
          chromTable_[j>>blockShift]=thisChrom;
      }
      ++thisChrom;
    }

    // First entry of chromNames is null so need to miss it out
    // and add a null entry to getContigName_ to match
    vector<string>::const_iterator i(chromNames.begin());
    i++;
    getContigName_.push_back( new ContigNameFinderNull() );

    for (;i!=chromNames.end(); ++i)
    {
      getContigName_.push_back
      (ContigNameFinder::getContigNameFinder(directoryName, *i));
    }
  }

  ~MatchPositionTranslator()
  {
    for (vector<ContigNameFinder*>::iterator i(getContigName_.begin());
	 i!=getContigName_.end();++i)
      delete *i;
  } // ~dtor

  void operator()( const MatchPosition originalPos,
		 const char*& chromName,
		 const char*& contigName,
		 MatchPosition& outputPos)
  {

	MatchPosition thisBlock = originalPos >> blockShift;
	outputPos=originalPos-subtractTable_[thisBlock];
	int thisChrom = chromTable_[thisBlock];
	chromName=chromNames_[thisChrom].c_str();
	contigName=chromName;
	assert((uint)thisChrom<getContigName_.size());
	(*getContigName_[thisChrom])( contigName, outputPos );
  } // ~operator()
 private:
  typedef MatchPosition SubtractTable[numPossibleChars];
  typedef int ChromTable[numPossibleChars];

  const vector<string>& chromNames_;
  ChromTable chromTable_;
  SubtractTable subtractTable_;
  vector<ContigNameFinder*> getContigName_;
}; // class MatchPositionTranslator


} // namespace eland_ms
} // namespace casava



#endif // CASAVA_ELAND_MS_MATCH_POSITION_TRANSLATOR_H
