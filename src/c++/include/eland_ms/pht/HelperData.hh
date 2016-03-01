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
 ** \file eland_ms/HelperData.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_PHT_HELPER_DATA_H
#define CASAVA_ELAND_MS_PHT_HELPER_DATA_H

#include <boost/static_assert.hpp>

namespace casava
{
namespace eland_ms
{



// Caches arguments for calls to MatchTable.addMatch()
//
// Helps alleviate apparently slow access to the match descriptor
// table during the first pass. Observed reduction in scan time is
// 30% for the default human alignment case on Xeon-5400ish cpus,
// however there is almost no effect on subsequent scans. This seems
// to have very little effect on the newer 5600's, presumably
// due to higher cpu<->ram bandwidth.
//
struct MatchCache {

    MatchCache(MatchTable& tab)
        : tab_(tab)
        , head_(0)
        , cache_(cacheSize)
    {}

    ~MatchCache() { processMatches(); }

    inline
    MatchStore&
    setNewMatch() {
        if(head_==cacheSize) { processMatches(); }
        return cache_[head_++];
    }

private:
    void
    processMatches() {
        tab_.addMatch(cache_.begin(),cache_.begin()+head_);
        head_=0;
    } 


    enum { cacheSize = 200 };
    
    MatchTable& tab_;
    unsigned head_;
    MatchCacheStore cache_;
};




typedef uint TablePointer;

// struct HashTableDataStore
// This is a wrapper for the data used by PartitionHashTable. Idea is
// this persists between passes so saves unnecessary allocation/dellocation
// of memory
template <bool useSplitPrefix>
struct HashTableDataStore
{
  void clear( void )
  {
    //    entryPointer_.clear();
    //    hashRem_.clear();

    vector<TablePointer> clear_entryPointer_;
    vector<TableEntry<useSplitPrefix> > clear_hashRem_;

    /*     entryPointer_.clear(); */
    /*     hashRem_.clear(); */
    entryPointer_.swap( clear_entryPointer_ );
    hashRem_.swap( clear_hashRem_ );

    entryPointer_.clear();
    hashRem_.clear();



  } // ~clear

  vector<TablePointer> entryPointer_;
  vector<TableEntry<useSplitPrefix> > hashRem_;
}; // ~struct HashTableDataStore



typedef std::map<Word,uint32_t> MaskMapType;



template <bool useSplitPrefix>
struct PHTHelperData
{
  PHTHelperData( HashTableDataStore<useSplitPrefix>& data, MatchTable& results ) :
    data_(data), results_(results) {}

  HashTableDataStore<useSplitPrefix>& data_;
  TablePointer* pCount_;
  int lowerFragSize_;
  Word lowerFragMask_;
  const FragmentErrorType* lowerFragScore_;
  const FragmentErrorType* upperFragScore_;
  MatchTable& results_;
  vector<Word> maskTable_;
  //   next two quantities only used if useSplitPrefix==true
  Word splitPrefixMask_;
  int splitPrefixShift_;
};

template<int useSplitPrefix> struct PHTHelperPreBase :
public PHTHelperData<useSplitPrefix>
{
  PHTHelperPreBase( HashTableDataStore<useSplitPrefix>& data, MatchTable& results ) :
      PHTHelperData<useSplitPrefix>( data, results ) {}

  BOOST_STATIC_ASSERT(useSplitPrefix==2); // allow speciallizations only

  void countKey( MaskMapType& maskMap, const Word key, const Word prefixMask, const Word suffixMask);

  void hashEntry( const MaskMapType& maskMap, Word key, const Word keyMask,
                  const Word entry, const Word entryMask, OligoNumber oligoNum );

  void
  setTopPrefix(const uint32_t tableSize);
};

template<class Child, int useSplitPrefix> struct PHTHelperBase :
public PHTHelperPreBase<useSplitPrefix>
{
  PHTHelperBase( HashTableDataStore<useSplitPrefix>& data, MatchTable& results ) :
      PHTHelperPreBase<useSplitPrefix>( data, results ) {}

  BOOST_STATIC_ASSERT(useSplitPrefix==2); // allow speciallizations only

  void check(MatchCache& cache, Word prefix, const Word suff, const MatchPosition sequencePos);
};

template<int PASS, bool isFwd, int useSplitPrefix> struct PHTHelper
{
    BOOST_STATIC_ASSERT(useSplitPrefix==2); // allow specialializations only
};


} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_PHT_HELPER_DATA_H
