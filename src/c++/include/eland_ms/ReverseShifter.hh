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
 ** \file eland_ms/ReverseShifter.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_REVERSE_SHIFTER_H
#define CASAVA_ELAND_MS_REVERSE_SHIFTER_H

#include <boost/static_assert.hpp>

namespace casava
{
namespace eland_ms
{

template<int OLIGO_LEN> struct ReverseShifterBase
{
  ReverseShifterBase() :
    prefixMask_( (Word)0x3 ),
    suffixMask_( ((Word)0x3)<<(numBitsPerBase*(ElandConstants<OLIGO_LEN>::suffixLength-1)) ) {}

  const Word prefixMask_;
  const Word suffixMask_;
};

template< int DIFF, int OLIGO_LEN> struct ReverseShifter
{
    BOOST_STATIC_ASSERT(OLIGO_LEN == 0);
};

template<int OLIGO_LEN> class ReverseShifter<0, OLIGO_LEN> : private ReverseShifterBase<OLIGO_LEN>
{
public:
  void operator()( const Oligo& , Oligo& rc ) const
  {
    // prefix and suffix are same size, don't need to swap bits between elements
    rc.ui[0]>>=((maxBasesPerWord-ElandConstants<OLIGO_LEN>::suffixLength)*numBitsPerBase);
    rc.ui[1]>>=((maxBasesPerWord-ElandConstants<OLIGO_LEN>::prefixLength)*numBitsPerBase);
  }
};

template<int OLIGO_LEN> class ReverseShifter<1, OLIGO_LEN> : public ReverseShifterBase<OLIGO_LEN>
{
public:
  void operator()( const Oligo& ol, Oligo& rc ) const
  {
    //  cout << "move pref to suff" << endl;
    // move 1 base from prefix to suffix
    rc.ui[0]<<=numBitsPerBase;
    rc.ui[0]>>=((maxBasesPerWord-ElandConstants<OLIGO_LEN>::suffixLength)*numBitsPerBase);
    rc.ui[1]>>=((maxBasesPerWord-ElandConstants<OLIGO_LEN>::prefixLength)*numBitsPerBase);
    rc.ui[1]^=(ol.ui[1]&ReverseShifterBase<OLIGO_LEN>::prefixMask_);
    //  rc.ui[1]^=prefixMask_;
  }
};

template<int OLIGO_LEN> class ReverseShifter<-1, OLIGO_LEN> : public ReverseShifterBase<OLIGO_LEN>
{
public:
  void operator()( const Oligo& ol, Oligo& rc ) const
  {
    //  cout << "move suff to pref \n" ;
    // move 1 base from suffix to prefix
    rc.ui[0]>>=numBitsPerBase;
    rc.ui[0]>>=((maxBasesPerWord-ElandConstants<OLIGO_LEN>::suffixLength)*numBitsPerBase);
    rc.ui[1]>>=((maxBasesPerWord-ElandConstants<OLIGO_LEN>::prefixLength)*numBitsPerBase);
    rc.ui[0]|=(ol.ui[0]&ReverseShifterBase<OLIGO_LEN>::suffixMask_);
    rc.ui[0]^=ReverseShifterBase<OLIGO_LEN>::suffixMask_;
  }
};

} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_REVERSE_SHIFTER_H
