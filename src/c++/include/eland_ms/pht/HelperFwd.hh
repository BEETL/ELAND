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
 ** \file eland_ms/HelperFwd.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_PHT_HELPER_FWD_H
#define CASAVA_ELAND_MS_PHT_HELPER_FWD_H

#include "HelperDataStandard.hh"
#include "HelperDataSplit.hh"

namespace casava
{
namespace eland_ms
{

template <int useSplitPrefix> struct PHTHelper<0, true, useSplitPrefix> :
public PHTHelperBase<PHTHelper<0,true,useSplitPrefix>, useSplitPrefix>
{
  PHTHelper( HashTableDataStore<useSplitPrefix>& data, MatchTable& results ) :
      PHTHelperBase<PHTHelper<0,true,useSplitPrefix>,useSplitPrefix>( data, results )
  {}

// Functionality specific to the 6 partitions
  bool wantMatch( const FragmentErrorType errorLow, const FragmentErrorType errorHigh, const Word ) const
  {
    // pass 0 in forward direction: want any match with 2 differences
    // or less - only way this doesn't occur is if there are two
    // errors in one segment and at least one in the other

    //  return (!( ((errorLow>=twoErrors)&&(errorHigh>=oneError))
    //     ||((errorHigh>=twoErrors)&&(errorLow>=oneError))));

    return (!( ((errorLow>oneError__)&&(errorHigh>noErrors__))
          ||((errorHigh>oneError__)&&(errorLow>noErrors__))));

    //  return (numErrors<=2);

  }

#if 0
  Word getInfoN( const Word , const Word ) const
    { return noNonIncorps; }
#endif
};

template <int useSplitPrefix> struct PHTHelper<1, true, useSplitPrefix> :
        public PHTHelperBase<PHTHelper<1,true,useSplitPrefix>, useSplitPrefix>
  {
    PHTHelper( HashTableDataStore<useSplitPrefix>& data, MatchTable& results ) :
      PHTHelperBase<PHTHelper<1,true,useSplitPrefix>,useSplitPrefix>( data, results )
    {}

  bool wantMatch( const FragmentErrorType errorLow, const FragmentErrorType errorHigh, const Word mask ) const
  {
  // pass 1 in forward direction: want any match with 2 errors or less,
  // provided it has errors or Ns in 2 distinct segments (if all errors/Ns
  // are in a single segment it will be matched on pass 0.

  //  return (     (  ((errorLow!=0)||((mask&lowerFragMask_)!=0))
  //		  && ((errorHigh!=0)||((mask&(~lowerFragMask_))!=0)))
  //       && ( !( ((errorLow>=twoErrors)&&(errorHigh>=oneError))
  //       ||((errorHigh>=twoErrors)&&(errorLow>=oneError)))));

  return (     (  ((errorLow!=0)||((mask&this->lowerFragMask_)!=0))
		  && ((errorHigh!=0)||((mask&(~this->lowerFragMask_))!=0)))
	   && ( !( ((errorLow>oneError__)&&(errorHigh>noErrors__))
		   ||((errorHigh>oneError__)&&(errorLow>noErrors__)))));

  }

#if 0
  Word getInfoN( const Word suff, const Word mask ) const
    { return ( ((mask&0x3)!=0x3) ? noNonIncorps : (suff&0x3) ); }
#endif
};

template<int useSplitPrefix> struct PHTHelper<2, true, useSplitPrefix> :
    public PHTHelperBase<PHTHelper<2,true,useSplitPrefix>, useSplitPrefix>
  {
    PHTHelper( HashTableDataStore<useSplitPrefix>& data, MatchTable& results ) :
      PHTHelperBase<PHTHelper<2,true,useSplitPrefix>,useSplitPrefix>( data, results )
    {}

  bool wantMatch
( const FragmentErrorType errorLow,
  const FragmentErrorType errorHigh,
  const Word mask ) const
{
  // pass 2 in forward direction: want any match with 2 errors or less,
  // provided it has errors or Ns in 2 distinct segments (if all errors/Ns
  // are in a single segment it will be matched on pass 0.

  //  return (     (  ((errorLow!=0)||((mask&lowerFragMask_)!=0))
  //	  && ((errorHigh!=0)||((mask&(~lowerFragMask_))!=0)))
  //   && ( !( ((errorLow>=twoErrors)&&(errorHigh>=oneError))
  //	   ||((errorHigh>=twoErrors)&&(errorLow>=oneError)))));

  return (     (  ((errorLow!=0)||((mask&this->lowerFragMask_)!=0))
		  && ((errorHigh!=0)||((mask&(~this->lowerFragMask_))!=0)))
	   && ( !( ((errorLow>oneError__)&&(errorHigh>noErrors__))
		   ||((errorHigh>oneError__)&&(errorLow>noErrors__)))));

}

#if 0
  Word getInfoN( const Word , const Word ) const
    { return noNonIncorps; }
#endif
};


} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_PHT_HELPER_FWD_H
