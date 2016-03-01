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
 ** \file eland_ms/Hasher.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_HASHER_H
#define CASAVA_ELAND_MS_HASHER_H

#include <boost/static_assert.hpp>
#include "ElandConstants.hh"
#include "SuffixScoreTable.hh"

namespace casava
{
namespace eland_ms
{

// Hash tables work as follows
// i) Using 2 bits per base encoding convert prefix of oligo or genome fragment
// into a number
// ii) Use this number to index into entryPointer_ - look up table of
// TablePointers
// iii) Value of corresponding TablePointer gives you the position in hashRem_
// of the first TableEntry for that prefix

// NB this won't work for n=32
#define hashsize(n) ((uint)1<<(n))
#define hashmask(n) (hashsize(n)-1)


// no need to compile this part 75 times! :
class HasherCore {
protected:
   void  do_interspersed( const Oligo& ol,Oligo& out ) const
   {
       static const Word MASK_A(0x33333333);  // 00110011001100110011001100110011
       static const Word MASK_B(0xcccccccc);  // 11001100110011001100110011001100

       const Word tmpA_lower = (ol.ui[0] & MASK_A);
       const Word tmpA_upper = (ol.ui[1] & MASK_A)<<2;
       out.ui[0] = tmpA_lower | tmpA_upper;

       const Word tmpB_lower = (ol.ui[0] & MASK_B);
       const Word tmpB_upper = (ol.ui[1] & MASK_B)>>2;
       out.ui[1] = tmpB_lower | tmpB_upper;
   }
};


// class Hasher: scrambles an oligo into a hash value
template<int OLIGO_LEN> class HasherBase : public HasherCore
{
public:
  HasherBase();

  // These are logical AND masks for the four fragments
  const Word fragMaskA_;
  const Word fragMaskB_;
  const Word fragMaskC_;
  const Word fragMaskD_;

  //  const int hashShift1_; - moved to ConvertOligo
  const int hashShift2_;

//  int getLowerFragSizePart1( void ) const;
//  Word getLowerFragMaskPart1( void ) const;
//  const FragmentErrorType*  getLowerFragScorePart1
//  ( const SuffixScoreTable& table ) const;
//  const FragmentErrorType*  getUpperFragScorePart1
//  ( const SuffixScoreTable& table ) const;
//  int getLengthPart1( void ) const;
//
//  int getLowerFragSizePart2( void ) const;
//  Word getLowerFragMaskPart2( void ) const;
//  const FragmentErrorType*  getLowerFragScorePart2
//  ( const SuffixScoreTable& table ) const;
//  const FragmentErrorType*  getUpperFragScorePart2
//  ( const SuffixScoreTable& table ) const;
//  int getLengthPart2( void ) const;
//
//    // do the interspersed hashing
//    inline void do_interspersed( const Oligo& ol,Oligo& out ) const;


private:
}; // ~class Hasher

template<int OLIGO_LEN> HasherBase<OLIGO_LEN>::HasherBase() :
  fragMaskA_(hashmask(numBitsPerBase*ElandConstants<OLIGO_LEN>::fragLengthA)),
  fragMaskB_(hashmask(numBitsPerBase*ElandConstants<OLIGO_LEN>::fragLengthB)
	     <<(numBitsPerBase*ElandConstants<OLIGO_LEN>::fragLengthA)),
  fragMaskC_(hashmask(numBitsPerBase*ElandConstants<OLIGO_LEN>::fragLengthC)),
  fragMaskD_(hashmask(numBitsPerBase*ElandConstants<OLIGO_LEN>::fragLengthD)
	     <<(numBitsPerBase*ElandConstants<OLIGO_LEN>::fragLengthC)),
  hashShift2_(numBitsPerBase*((OLIGO_LEN%ElandConstants<OLIGO_LEN>::fragmentsPerOligo)==2))
{
  //  cerr << fragLengthA << " " << fragLengthB << " "
  //      << fragLengthC << " " << fragLengthD << endl;

  BOOST_STATIC_ASSERT(ElandConstants<OLIGO_LEN>::fragLengthA
                    + ElandConstants<OLIGO_LEN>::fragLengthB
                    + ElandConstants<OLIGO_LEN>::fragLengthC
                    + ElandConstants<OLIGO_LEN>::fragLengthD ==OLIGO_LEN);
  assert((fragMaskA_&fragMaskB_)==0);
  assert((fragMaskC_&fragMaskD_)==0);

} // ~template<> Hasher<0>::Hasher( int keyBits, int oligoBases ) :

template<int PASS, int OLIGO_LEN> class Hasher
{
  BOOST_STATIC_ASSERT(OLIGO_LEN == 0);
};

template<int OLIGO_LEN> class Hasher<0, OLIGO_LEN> : private HasherBase<OLIGO_LEN>
{
public:
  // Parameters for pass 0: partition1=AB, partition2=CD
  int getLowerFragSizePart1( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthA; }
  Word getLowerFragMaskPart1( void ) const
  { return HasherBase<OLIGO_LEN>::fragMaskA_; }
  const FragmentErrorType* getLowerFragScorePart1
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragA_[0]; }
  const FragmentErrorType* getUpperFragScorePart1
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragB_[0]; }
  int getLengthPart1( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthA+ElandConstants<OLIGO_LEN>::fragLengthB; }
  int getLowerFragSizePart2( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthC; }
  Word getLowerFragMaskPart2( void ) const
  { return HasherBase<OLIGO_LEN>::fragMaskC_; }
  const FragmentErrorType* getLowerFragScorePart2
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragC_[0]; }
  const FragmentErrorType* getUpperFragScorePart2
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragD_[0]; }
  int getLengthPart2( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthC+ElandConstants<OLIGO_LEN>::fragLengthD; }

  void operator()( const Oligo& ol, Oligo& out ) const;
};

template<int OLIGO_LEN> class Hasher<1, OLIGO_LEN> : private HasherBase<OLIGO_LEN>
{
public:
  // Parameters for pass 1: partition1=CB, partition2=AD
  int getLowerFragSizePart1( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthC; }
   Word getLowerFragMaskPart1( void ) const
  { return HasherBase<OLIGO_LEN>::fragMaskC_; }
  const FragmentErrorType*  getLowerFragScorePart1
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragC_[0]; }
  const FragmentErrorType*  getUpperFragScorePart1
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragB_[0]; }
   int getLengthPart1( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthC+ElandConstants<OLIGO_LEN>::fragLengthB; }

  int getLowerFragSizePart2( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthA; }
   Word getLowerFragMaskPart2( void ) const
  { return HasherBase<OLIGO_LEN>::fragMaskA_; }
  const FragmentErrorType*  getLowerFragScorePart2
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragA_[0]; }
  const FragmentErrorType*  getUpperFragScorePart2
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragD_[0]; }
   int getLengthPart2( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthA+ElandConstants<OLIGO_LEN>::fragLengthD; }

   void operator()( const Oligo& ol, Oligo& out ) const;
};

template<int OLIGO_LEN> class Hasher<2, OLIGO_LEN> : private HasherBase<OLIGO_LEN>
{
public:
  // Parameters for pass 2: partition1=DB, partition2=CA
  int getLowerFragSizePart1( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthD; }
   Word getLowerFragMaskPart1( void ) const
  { return (HasherBase<OLIGO_LEN>::fragMaskD_>>(numBitsPerBase*ElandConstants<OLIGO_LEN>::fragLengthC)); }
  const FragmentErrorType*  getLowerFragScorePart1
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragD_[0]; }
  const FragmentErrorType*  getUpperFragScorePart1
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragB_[0]; }
   int getLengthPart1( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthD+ElandConstants<OLIGO_LEN>::fragLengthB; }

  // ordering of A and C swapped so that A is in the lower bits of the
  // fragment
  int getLowerFragSizePart2( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthA; }
   Word getLowerFragMaskPart2( void ) const
  { return HasherBase<OLIGO_LEN>::fragMaskA_; }
  const FragmentErrorType*  getLowerFragScorePart2
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragA_[0]; }
  const FragmentErrorType*  getUpperFragScorePart2
  ( const SuffixScoreTable& table ) const
  { return &table.scoreFragC_[0]; }
   int getLengthPart2( void ) const
  { return ElandConstants<OLIGO_LEN>::fragLengthC+ElandConstants<OLIGO_LEN>::fragLengthA; }

   void operator()( const Oligo& ol, Oligo& out ) const;
};


// hash0: just convert an oligo into two partitions without scrambling,
// e.g. for a 24-mer
// LSB    ----------------->    MSBlsb    ----------------->    msb
// ol[0]                           ol[1]
// 2 2 2 2 1 1 1 1 1 1 1 1 1 1
// 3.2.1.0.9.8.7.6.5.4.3.2.1.0.9.8.7.6.5.4.3.2.1.0.xxxxxxxxxxxxxxxx
// ->
// out[0]                          out[1]
// 2 2 2 2 1 1 1 1 1 1 1 1         1 1
// 3.2.1.0.9.8.7.6.5.4.3.2.xxxxxxxx1.0.9.8.7.6.5.4.3.2.1.0.xxxxxxxx
// A---------->B---------->        C---------->D---------->

template<int OLIGO_LEN> void Hasher<0, OLIGO_LEN>::operator()( const Oligo& ol, Oligo& out ) const
{
  BOOST_STATIC_ASSERT(32 != OLIGO_LEN); // Threre is a specialized version for this
  //  printWord(ol.ui[1],16);
  //  cout << '-';
  //  printWord(ol.ui[0],16);
  // cout << " -> ";

    out.ui[0]=ol.ui[0];
    out.ui[1]=ol.ui[1];

  //  printWord(out.ui[1],16);
  //  cout << '-';
  //  printWord(out.ui[0],16);
  //  cout << endl;

  //  cout << "done for this oligo" << endl << endl;

}
template<> void Hasher<0, 32>::operator()( const Oligo& ol, Oligo& out ) const;

// hash1: convert an oligo into two partitions and scramble as follows,
// e.g. for a 24-mer
// LSB    ----------------->    MSBlsb    ----------------->    msb
// ol[0]                           ol[1]
// 2 2 2 2 1 1 1 1 1 1 1 1 1 1
// 3.2.1.0.9.8.7.6.5.4.3.2.1.0.9.8.7.6.5.4.3.2.1.0.xxxxxxxxxxxxxxxx
// ->
// out[0]                          out[1]
// 2 2 2 2 1 1 1 1                 1 1 1 1 1 1
// 3.2.1.0.9.8.1.0.9.8.7.6.xxxxxxxx7.6.5.4.3.2.5.4.3.2.1.0.xxxxxxxx
// A---------->C---------->        B---------->D---------->


template<int OLIGO_LEN> void Hasher<1, OLIGO_LEN>::operator()( const Oligo& ol, Oligo& out ) const
{
  BOOST_STATIC_ASSERT(32 != OLIGO_LEN); // Threre is a specialized version for this
  //  printWord(ol.ui[1],16);
  //  cout << '-';
  //  printWord(ol.ui[0],16);
  //  cout << " -> ";

  out.ui[0]=ol.ui[0]&HasherBase<OLIGO_LEN>::fragMaskB_; // out.ui[0] = 0B
  out.ui[0]|=(ol.ui[1]&HasherBase<OLIGO_LEN>::fragMaskC_);// out.ui[0] = CB
  out.ui[1]=ol.ui[1]&HasherBase<OLIGO_LEN>::fragMaskD_; // out.ui[1] = 0D
  //  out.ui[1]|=(ol.ui[1]&fragMaskA_); // out.ui[1] = AD - oops! AC 6.6.3
  out.ui[1]|=(ol.ui[0]&HasherBase<OLIGO_LEN>::fragMaskA_); // out.ui[1] = AD

  //  printWord(out.ui[1],16);
  //  cout << '.';
  //  printWord(out.ui[0],16);
  //  cout << endl;

}
template<> void Hasher<1, 32>::operator()( const Oligo& ol, Oligo& out ) const;

// hash2: convert an oligo into two partitions and scramble as follows,
// e.g. for a 24-mer
// LSB    ----------------->    MSBlsb    ----------------->    msb
// ol[0]                           ol[1]
// 2 2 2 2 1 1 1 1 1 1 1 1 1 1
// 3.2.1.0.9.8.7.6.5.4.3.2.1.0.9.8.7.6.5.4.3.2.1.0.xxxxxxxxxxxxxxxx
// ->
// out[0]                          out[1]
// 2 2 2 2 1 1                     1 1         1 1 1 1 1 1
// 3.2.1.0.9.8.5.4.3.2.1.0.xxxxxxxx1.0.9.8.7.6.7.6.5.4.3.2.xxxxxxxx
// A---------->D---------->        C---------->B---------->

template<int OLIGO_LEN> void Hasher<2, OLIGO_LEN>::operator()( const Oligo& ol, Oligo& out ) const
{
  BOOST_STATIC_ASSERT(32 != OLIGO_LEN); // Threre is a specialized version for this
  //  printWord(ol.ui[1],16);
  //  cout << '-';
  //  printWord(ol.ui[0],16);
  //  cout << " -> ";

   out.ui[0]=ol.ui[0]&HasherBase<OLIGO_LEN>::fragMaskB_; // out.ui[0] = 0B
   out.ui[0]<<=HasherBase<OLIGO_LEN>::hashShift2_; // shift B if necessary so D fits next to it
   out.ui[0]|=ol.ui[1]>>(numBitsPerBase*ElandConstants<OLIGO_LEN>::fragLengthC); // out.ui[0]=DB

   // ordering of A and C swapped in second fragment TC 03.07.03
   //  out.ui[1]&=fragMaskC_; // out.ui[1]=C0
   // out.ui[1]|=((ol.ui[0]&fragMaskA_)<<(numBitsPerBase*fragLengthC));
   out.ui[1]=ol.ui[1]&HasherBase<OLIGO_LEN>::fragMaskC_;
   out.ui[1]<<=(numBitsPerBase*ElandConstants<OLIGO_LEN>::fragLengthA);
   out.ui[1]|=(ol.ui[0]&HasherBase<OLIGO_LEN>::fragMaskA_);


  //  printWord(out.ui[1],16);
  //  cout << '-';
  //  printWord(out.ui[0],16);
  //  cout << endl;

}
template<> void Hasher<2, 32>::operator()( const Oligo& ol, Oligo& out ) const;

} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_HASHER_H
