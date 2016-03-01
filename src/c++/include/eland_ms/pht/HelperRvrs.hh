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
 ** \file eland_ms/HelperRvrs.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_PHT_HELPER_RVRS_H
#define CASAVA_ELAND_MS_PHT_HELPER_RVRS_H

#include "HelperDataStandard.hh"
#include "HelperDataSplit.hh"

namespace casava
{
namespace eland_ms
{

template<int useSplitPrefix> struct PHTHelper<0, false, useSplitPrefix> :
    public PHTHelperBase<PHTHelper<0,false, useSplitPrefix>, useSplitPrefix>
{
  PHTHelper( HashTableDataStore<useSplitPrefix>& data, MatchTable& results ) :
    PHTHelperBase<PHTHelper<0,false, useSplitPrefix>,useSplitPrefix>( data, results )
  {}

  bool wantMatch( const FragmentErrorType errorLow, const FragmentErrorType errorHigh, const Word mask ) const
  {
    // pass 0 in reverse direction: want any match with 1 or 2 errors,
    // only want zero error matches if mask is non zero

    //  return ( ((errorLow!=0)||(errorHigh!=0)||(mask!=0))
    //   &&( !( ((errorLow>=twoErrors)&&(errorHigh>=oneError))
    //     ||((errorHigh>=twoErrors)&&(errorLow>=oneError)))));

    return ( ((errorLow!=0)||(errorHigh!=0)||(mask!=0))
	   &&( !( ((errorLow>oneError__)&&(errorHigh>noErrors__))
          ||((errorHigh>oneError__)&&(errorLow>noErrors__)))));

    //  return ( ((numErrors==1)||(numErrors==2))
    //	   ||((numErrors==0)&(mask!=0)) );
  }

#if 0
  Word getInfoN( const Word suff, const Word mask ) const
    { return ( ((mask&0x3)!=0x3) ? noNonIncorps : (suff&0x3) ); }
#endif
};

template<int useSplitPrefix> struct PHTHelper<1, false, useSplitPrefix> :
    public PHTHelperBase<PHTHelper<1,false, useSplitPrefix>, useSplitPrefix>
  {
    PHTHelper( HashTableDataStore<useSplitPrefix>& data, MatchTable& results ) :
      PHTHelperBase<PHTHelper<1,false, useSplitPrefix>, useSplitPrefix>( data, results )
    {}

  bool wantMatch( const FragmentErrorType errorLow, const FragmentErrorType errorHigh, const Word mask ) const
  {
  // pass 1 in reverse direction: want any match with 2 errors or less,
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

template<int useSplitPrefix> struct PHTHelper<2, false, useSplitPrefix> :
    public PHTHelperBase<PHTHelper<2,false, useSplitPrefix>, useSplitPrefix>
  {
    PHTHelper( HashTableDataStore<useSplitPrefix>& data, MatchTable& results ) :
      PHTHelperBase<PHTHelper<2,false, useSplitPrefix>, useSplitPrefix>( data, results )
    {}

  bool wantMatch
( const FragmentErrorType errorLow,
  const FragmentErrorType errorHigh,
  const Word mask ) const
{
  // pass 2 in reverse direction: want any match with 2 errors or less,
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
  Word getInfoN( const Word suff, const Word mask ) const
    { return ( ((mask&0x3)!=0x3) ? noNonIncorps : (suff&0x3) );}
#endif
};


} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_PHT_HELPER_RVRS_H
