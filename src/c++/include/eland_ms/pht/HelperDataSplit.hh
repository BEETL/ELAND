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
 ** \file eland_ms/HelperDataSplit.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_PHT_HELPER_DATA_SPLIT_H
#define CASAVA_ELAND_MS_PHT_HELPER_DATA_SPLIT_H

#include "HelperData.hh"

namespace casava
{
namespace eland_ms
{


// Specialization for split prefix mode
template<> struct PHTHelperPreBase<true> :
public PHTHelperData<true>
{
  PHTHelperPreBase( HashTableDataStore<true>& data, MatchTable& results )
      : PHTHelperData<true>( data, results )
      , topShiftIn_(0)
      , topShiftOut_(0)
      , topMask_(~0)
  {}

  void countKey( MaskMapType& maskMap, const Word key, const Word prefixMask, const Word suffixMask)
  {
    if (prefixMask==0)
    {
      maskMap[suffixMask] = 0;
      ++pCount_[key>>splitPrefixShift_];
    }
    //    pCount_[key>>splitPrefixShift_]+=(mask==0);
  }

  void hashEntry( const MaskMapType& maskMap, Word key, const Word keyMask, const Word entry,
                  const Word entryMask, OligoNumber oligoNum )
  {
    if (0 == keyMask)
    {
      //      oligoNum|=((key&splitPrefixMask_)<<MAX_HASH_BITS);
      PrefixType prefix((PrefixType)(key&splitPrefixMask_));
      key >>= splitPrefixShift_;
      data_.hashRem_[pCount_[key]].prefix=prefix;
      data_.hashRem_[pCount_[key]].suffix.ui=entry;
      data_.hashRem_[pCount_[key]].mask=(maskMap.find(entryMask))->second;
      data_.hashRem_[pCount_[key]].position=oligoNum;
      pCount_[key]++;
    }
  }

  // add first prefix into the top bits of the table entry pointers 
  //
  void
  setTopPrefix(const uint32_t tableSize)
  {
    // first determine how many bits of the prefix are available: 
    static const unsigned tpBits(sizeof(TablePointer)*numBitsPerByte);
    uint8_t zBits(splitPrefixShift_);
    { 
      const TablePointer maxVal(pCount_[tableSize]);
      for(TablePointer i(maxVal);zBits && (i>>(tpBits-zBits));zBits--) {}
    }

    // topShiftIn - the amount of prefix information we loose -- this has to be taken from the LSB:
    topShiftIn_=(splitPrefixShift_-zBits);

    // topShiftOut - the offset of the (possibly LSB trimmed) prefix into the TablePointer: 
    topShiftOut_=(tpBits-zBits);

    for (TablePointer i(0) ; i < tableSize ; i++ ) {
      if ((pCount_[i+1]==pCount_[i])) continue;
      pCount_[i] |= ((data_.hashRem_[pCount_[i]].prefix>>topShiftIn_)<<topShiftOut_);
    }

    topMask_=((static_cast<uint64_t>(1)<<topShiftOut_)-1);
  }

protected:
    uint8_t topShiftIn_;
    uint8_t topShiftOut_;
    uint32_t topMask_;
};


// Specialization for split prefix mode
template<class Child> struct PHTHelperBase<Child, true> :
public PHTHelperPreBase<true>
{
  PHTHelperBase( HashTableDataStore<true>& data, MatchTable& results ) :
    PHTHelperPreBase<true>( data, results )
  {}

  void check(MatchCache& cache, Word prefix, const Word suff, const MatchPosition sequencePos )
  {
    //    printWord( prefix, 12 ); cout << "****" << endl;
    register Word thisMatch,thisMask;
    register FragmentErrorType errorLow,errorHigh;
    const PrefixType thisSplitPrefix((PrefixType)(prefix&splitPrefixMask_));
    //    OligoNumber thisSplitPrefix((prefix&splitPrefixMask_)<<MAX_HASH_BITS);
    prefix>>=splitPrefixShift_;


    uint i(topMask_&pCount_[prefix]);
    const uint i_end(topMask_&pCount_[prefix+1]);

    if((i==i_end) ||
       (((i+1)==i_end) && ((static_cast<uint32_t>(thisSplitPrefix)>>topShiftIn_)!=(pCount_[prefix]>>topShiftOut_))) || 
       ((static_cast<uint32_t>(thisSplitPrefix)>>topShiftIn_)<(pCount_[prefix]>>topShiftOut_)) ) return;

    for (; i!=i_end; ++i)
    {
      //      if ((data_.hashRem_[i].position&splitPrefixMaskHigh)
      //	    ==thisSplitPrefix) break;
      if (data_.hashRem_[i].prefix==thisSplitPrefix){
        break;
      }
    } // ~for

    for (;i!=i_end;++i)
    {
      //if ((data_.hashRem_[i].position&splitPrefixMaskHigh)!=thisSplitPrefix)
      //break;
      if (data_.hashRem_[i].prefix!=thisSplitPrefix) {
        break;
      }

      thisMask=0;
      thisMatch=(suff^data_.hashRem_[i].suffix.ui);
      if(data_.hashRem_[i].mask){
          thisMask=(maskTable_[data_.hashRem_[i].mask]);
          thisMatch&=(~thisMask);
       }
      if (((errorLow=(lowerFragScore_[thisMatch&lowerFragMask_])) < moreThanTwoErrors__) && 
           ((errorHigh=(upperFragScore_[thisMatch>>(numBitsPerBase*lowerFragSize_)])) < moreThanTwoErrors__))
      {
          if (static_cast<const Child*>(this)->wantMatch(errorLow, errorHigh, thisMask))
          {
              cache.setNewMatch().set(data_.hashRem_[i].position,sequencePos,
                                      ((errorLow>oneError__)+(errorLow>noErrors__)
                                       +(errorHigh>oneError__)+(errorHigh>noErrors__)));
          }
      }
    }
  }

//  void check__( const vector<pair<Oligo, MatchPosition> >& h )
//    {}
};

} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_PHT_HELPER_DATA_SPLIT_H
