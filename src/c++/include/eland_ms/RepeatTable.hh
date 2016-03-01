/**
 ** copyright (c) 2007-2010 Illumina, Inc.
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
 ** \file eland_ms/RepeatTable.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_REPEAT_TABLE_H
#define CASAVA_ELAND_MS_REPEAT_TABLE_H

#include "QueryGenerator.hh"

namespace casava
{
namespace eland_ms
{

// class RepeatTable: reads in an ASCII table of repeats and filters them out
// of the incoming oligos.
template <int OLIGO_LEN> class RepeatTable
{
public:

  // struct IsLess: returns true if the unmasked bits of a are < b
  // allows binary searching for oligos containing Ns
  struct IsLess
  {
    IsLess( const Oligo& mask ) :
      mask0_(~mask.ui[0]), mask1_(~mask.ui[1]) {}
    bool operator()( const Oligo& a, const Oligo& b )
    {
      return( ((a.ui[1]&mask1_)<(b.ui[1]&mask1_))
      	      ||( ((a.ui[1]&mask1_)==(b.ui[1]&mask1_))
		  &&((a.ui[0]&mask0_)<(b.ui[0]&mask0_))) );
    } // ~operator()
    const Word mask0_;
    const Word mask1_;
  }; // ~struct isLess

//  RepeatTable( const char* repeatFileName );
  ~RepeatTable()
  {
    cerr << "deleting repeat table" << endl;
  } // ~dtor

  bool check( const Oligo& oligo, const Oligo& mask )
  {
    return (!(binary_search(oligos_.begin(),
			    oligos_.end(),
			    oligo,
			    IsLess(mask) )));
  } // ~check

//  void checkOligos( OligoSource& oligos, MatchTable<OLIGO_LEN>& results );

private:
  vector<Oligo> oligos_;
public:

RepeatTable( const char* repeatFileName )
{
  OligoSourceRaw repeats(repeatFileName);
  Oligo last;
  uint numRepeats(0), numEntries(0);
  QueryGenerator<OLIGO_LEN> makeOligos;
  const char* pRepeat;

  cerr << "creating repeat table" << endl;
  last.ui[1]=~0;


  oligos_.push_back( Oligo() );
  while (1)
  {
    ++numEntries;
    if ((pRepeat=repeats.getNextOligo())==NULL) break;
    makeOligos( pRepeat, oligos_.back() );
    if (oligos_.back()!=last)
    {
      ++numRepeats;
      if (numEntries!=1) assert((last<oligos_.back())||(last==oligos_.back()));
      last=oligos_.back();
      oligos_.push_back( Oligo() );
    } // ~if
    //   else cout << " REPEATED";
    //   cout << endl;
  } // ~while

  oligos_.pop_back();

  sort(oligos_.begin(), oligos_.end());
  cerr << "storing " << numRepeats << " repeats out of "
       << numEntries << " entries" << endl;
  assert(numRepeats==oligos_.size());
} // ~void RepeatTable:RepeatTable( const char* repeatFileName )

// return false if a repeat, true if OK to hash
void checkOligos( OligoSource& oligos, MatchTable& results )
{
  const char* pOligo;
  Oligo oligo, mask;
  uint numRepeats(0), numEntries(0);
  QueryGenerator<OLIGO_LEN> makeOligos;
  cerr << "checking oligos for repeats" << endl;
  while (1)
  {
      if ((pOligo=oligos.getNextOligoSelect(false,false))==NULL) break;
    ++numEntries;

    // bug fix: oligo & mask need to be zeroed before each encoding - TC 18.7.3
    // (can get away with this in most cases because oligos are not reused)
    mask.ui[0]=0;
    mask.ui[1]=0;
    oligo.ui[0]=0;
    oligo.ui[1]=0;

    // changed to allow oligos with Ns to be checked - TC 11.07.03
    makeOligos.encodeOligo( pOligo, oligo, mask );
    if (check( oligo, mask )==false)
    {
      ++numRepeats;
      results.resize(numEntries+1);

      //      results.matchPosition_[numEntries]=repeatMasked;
      results.setRepeatMasked__(numEntries);

    } // ~if
  } // ~while

  oligos.rewind();
  cerr << "found " << numRepeats << " repeats in " << numEntries << " oligos" << endl;
  results.resize(numEntries+1);

}

};

} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_REPEAT_TABLE_H
