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
 ** \file eland_ms/HelperDataStandard.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_PHT_HELPER_DATA_STANDARD_H
#define CASAVA_ELAND_MS_PHT_HELPER_DATA_STANDARD_H

#include "HelperData.hh"

namespace casava
{
namespace eland_ms
{


// Specialization for "standard" mode (prefix not split)
template<> struct PHTHelperPreBase<false> :
public PHTHelperData<false>
{
  PHTHelperPreBase( HashTableDataStore<false>& data, MatchTable& results ) :
    PHTHelperData<false>( data, results )
  {}

  void countKey( MaskMapType& maskMap, const Word key, const Word prefixMask, const Word suffixMask)
  {
    if (prefixMask==0)
    {
      maskMap[suffixMask] = 0;
      ++pCount_[key];
    }
  }

  void hashEntry( const MaskMapType& maskMap, Word key, const Word keyMask, const Word entry,
                  const Word entryMask, OligoNumber oligoNum )
  {
    if (0 == keyMask)
    {
      data_.hashRem_[pCount_[key]].suffix.ui=entry;
      data_.hashRem_[pCount_[key]].mask=(maskMap.find(entryMask))->second;
      data_.hashRem_[pCount_[key]].position=oligoNum;
      pCount_[key]++;
    }
  }

  void
  setTopPrefix(const uint32_t /*tableSize*/) { }

};


// Specialization for "standard" mode (prefix not split)
template<class Child> struct PHTHelperBase<Child, false> :
public PHTHelperPreBase<false>
{
  PHTHelperBase( HashTableDataStore<false>& data, MatchTable& results ) :
    PHTHelperPreBase<false>( data, results )
  {}

  void check(MatchCache& cache, Word prefix, const Word suff, const MatchPosition sequencePos )
  {
    register Word thisMatch,thisMask;
    register FragmentErrorType errorLow,errorHigh;
    uint i(pCount_[prefix]);
    for (; i!=pCount_[prefix+1]; ++i)
    {
#ifdef DEBUG
      printWord(data_.hashRem_[i].suffix.ui, maxBasesPerWord );
      cout << " - ";
      printWord(maskTable_[data_.hashRem_[i].mask],maxBasesPerWord);
      cout << endl;
#endif

      thisMask=0;
      thisMatch=(suff^data_.hashRem_[i].suffix.ui);
      if(data_.hashRem_[i].mask){
          thisMask=(maskTable_[data_.hashRem_[i].mask]);
          thisMatch&=(~thisMask);
      }

      if (((errorLow=(lowerFragScore_[thisMatch&lowerFragMask_])) < moreThanTwoErrors__) &&
           ((errorHigh=(upperFragScore_[thisMatch>>(numBitsPerBase*lowerFragSize_)])) < moreThanTwoErrors__))
      {
#ifdef DEBUG
          numErrors=errorLow;
          numErrors&=(~(errorMask1|errorMask2));
          numErrors+=errorHigh;
          numErrors>>=(2*errorBits);
          cout << (data_.hashRem_[i].position&(~isReverseOligo)) << " "
               << ((data_.hashRem_[i].position&isReverseOligo)?'R':'F') << " "
               << numErrors << " " << sequencePos-blockSize << " ";
          printWord(prefix,10);  cout << " ";
          printWord(suff, 10);
#endif
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
//  {
//    register Word thisMatch(0), thisMask;
//
//#ifdef DEBUG
//    cout << sequencePos-blockSize << " - position " << endl;
//    printWord( prefix, maxBasesPerWord ); cout << " - prefix" << endl;
//    printWord( suff, maxBasesPerWord); cout << " - suffix" << endl;
//#endif
//
//    vector<pair<uint, uint> >p(h.size());
//    for (uint j(0);j<h.size();++j)
//    {
//      p[j].first=this->pCount_[h[j].first.ui[0]];
//      p[j].second=this->pCount_[h[j].first.ui[0]+1];
//    } // ~for j
//
//    register FragmentErrorType errorLow, errorHigh;
//
//    for (uint j(0);j<p.size();++j)
//    {
//      for (uint i(p[j].first); i!=p[j].second; ++i)
//      {
//        //    for (uint i(pCount_[prefix]);
//        //   i!=pCount_[prefix+1];
//        //   ++i)
//
//        //      cout << "XX " << b-a << endl;
//#ifdef DEBUG
//        printWord(data_.hashRem_[i].suffix.ui, maxBasesPerWord );
//        cout << " - ";
//        printWord(maskTable_[data_.hashRem_[i].mask],maxBasesPerWord);
//        cout << endl;
//#endif
//
//        thisMatch = (h[j].first.ui[1]^this->data_.hashRem_[i].suffix.ui);
//
//        thisMask = this->maskTable_[this->data_.hashRem_[i].mask];
//        //      thisMatch &= (~data_.hashRem_[i].mask.ui);
//        thisMatch &= (~thisMask);
//        //      printWord(thisMatch,10); cout << endl;
//
//        //      errorLow=lowerFragScore_[thisMatch&lowerFragMask_];
//        if ((errorLow=this->lowerFragScore_[thisMatch&this->lowerFragMask_]) < moreThanTwoErrors__)
//        {
//          thisMatch>>=(numBitsPerBase*this->lowerFragSize_); // eliminate temporary
//          if ((errorHigh=this->upperFragScore_[thisMatch]) < moreThanTwoErrors__)
//          {
//
//            // numErrors=(errorLow>>(2*errorBits))+(errorHigh>>(2*errorBits));
//#ifdef DEBUG
//            numErrors=errorLow;
//            numErrors&=(~(errorMask1|errorMask2));
//            numErrors+=errorHigh;
//            numErrors>>=(2*errorBits);
//            cout << (data_.hashRem_[i].position&(~isReverseOligo)) << " "
//                 <<((data_.hashRem_[i].position&isReverseOligo)?'R':'F') << " "
//                 << numErrors << " " << sequencePos-blockSize << " ";
//            printWord(prefix,10);  cout << " ";
//            printWord(suff, 10);
//#endif
//            //            cout << "ck " << sequencePos-16777216 << " " << (errorLow&errorPosMask1__) << " "
//            //   << ((errorLow&errorPosMask2__)>>8) << " "
//            //   << (errorHigh&errorPosMask1__) << " "
//            //  << ((errorHigh&errorPosMask2__)>>8) << endl;
//
//
//            if (wantMatch(errorLow, errorHigh, thisMask))
//            {
//              //              cout << "wtd " << sequencePos-16777216 << " " << (errorLow&errorPosMask1__) << " "
//              //    << ((errorLow&errorPosMask2__)>>8) << " "
//              //  << (errorHigh&errorPosMask1__) << " "
//              //  << ((errorHigh&errorPosMask2__)>>8) << endl;
//              this->results_.addMatch( this->data_.hashRem_[i], thisMask, h[j].second, errorLow, errorHigh );
//            } // ~if
//
//          } // ~if
//        } // ~if
//      } // ~for	i
//    } // ~for j
//
//  }  // ~check__

};

} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_PHT_HELPER_DATA_STANDARD_H
