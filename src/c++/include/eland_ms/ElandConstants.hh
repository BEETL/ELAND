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
 ** \file eland_ms/ElandConstants.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 ** Contains interface functionality from ELAND - specifically, stuff to do
 ** with scoring and storage of alignments
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_ELAND_CONSTANTS_H
#define CASAVA_ELAND_MS_ELAND_CONSTANTS_H


#include <boost/assert.hpp>
#include <boost/integer/static_min_max.hpp>
#include "alignment/GlobalUtilities.hh"
#include "ElandDefines.hh"

namespace casava
{
namespace eland_ms
{

// Program-wide compile-time constants
template <unsigned int OLIGO_LEN> class ElandConstants
{
public:
  enum{
  // hashShift is the number of bits to shift for a given OLIGO_LEN
  // so that the number of bases stored in each of the two unsigned ints
  // is either equal (if OLIGO_LEN is even) or else (if OLIGO_LEN is
  // odd) differs only by 1.
  // hashShift = ((2*numBitsPerWord/numBitsPerBase)-OLIGO_LEN+1)&(~1),

  // an oligo is divided into four fragments
  fragmentsPerOligo = 4,
  //  basesPerFragment = (OLIGO_LEN/fragmentsPerOligo),
  //  bitsPerFragment = 2*basesPerFragment,

  // oligo is partitioned into a prefix and suffix, each containing two
  // fragments
  fragmentsPerPartition = 2,

  // there are 6 ways of partitioning the 4 fragments into a prefix and
  // a suffix (not 24 as we are not concerned about the ordering of the
  // fragments within the prefix and suffix). Two of these partitionings are
  // searched on each pass through the genome, hence 3 passes required in all
  numPasses = 3,

  // if more than this many matches of a given quality have already been found
  // for an oligo then don't look for any more, e.g. if we have more than 3
  // matches with 2 errors, then only single error or exact match errors
  // will be looked for on subsequent passes
  maxNumBestMatches = 3,

  // length (in bits) of largest suffix that will be needed for this oligo length
  // length 20 => fragments = 5 5 5 5 => largest suffix =5+5=10 bases = 20 bits
  // length 21 => fragments = 5 6 5 5 => largest suffix =5+6=11 bases = 22 bits
  // length 22 => fragments = 5 6 5 6 => largest suffix =6+6=12 bases = 24 bits
  // length 23 => fragments = 6 5 6 6 => largest suffix =6+6=12 bases = 24 bits
  // length 24 => fragments = 6 6 6 6 => largest suffix =6+6=12 bases = 24 bits
  maxPrefixLength = OLIGO_LEN
  +((OLIGO_LEN%fragmentsPerOligo)!=0)
  +((OLIGO_LEN%fragmentsPerOligo)==2),

  maxOligoLength = 2*maxBasesPerWord,

  // Oligos are split into fragments as follows:
  // Oligo Length  Fragment Lengths: A B C D
  // 20                              5 5 5 5
  // 21                              5 6 5 5
  // 22                              5 6 5 6
  // 23                              6 5 6 6
  // 24                              6 6 6 6
  // 25                              6 7 6 6 ... etc.
  fragLengthA=((OLIGO_LEN/fragmentsPerOligo)
	       + ((OLIGO_LEN%fragmentsPerOligo)==3)),
  fragLengthB=((OLIGO_LEN/fragmentsPerOligo)
	       + (    ((OLIGO_LEN%fragmentsPerOligo)==1)
		      ||((OLIGO_LEN%fragmentsPerOligo)==2))),
  fragLengthC=((OLIGO_LEN/fragmentsPerOligo)
	       + ((OLIGO_LEN%fragmentsPerOligo)==3)),
  fragLengthD=((OLIGO_LEN/fragmentsPerOligo)
	       + (    ((OLIGO_LEN%fragmentsPerOligo)==2)
		      ||((OLIGO_LEN%fragmentsPerOligo)==3))),
  prefixLength=fragLengthC+fragLengthD,
  suffixLength=fragLengthA+fragLengthB
  };

  BOOST_STATIC_ASSERT(OLIGO_LEN<=maxOligoLength);

// Software operates in two modes:
// 1. Standard mode, each oligo is split into prefix and suffix, prefix used to
// make a look up table
// 2. "Split prefix" mode: prefix is too large to be used to make a look up table
// so "hide" some bits in the TableEntry

static const bool useSplitPrefix = (maxPrefixLength>MAX_HASH_BITS);
//const int maxNumOligos((1<<(MAX_HASH_BITS-1))-1);

// next two should be OligoNumbers really?
static const Word splitPrefixMaskHigh = ( ((Word)~0) << MAX_HASH_BITS );
static const Word splitPrefixMaskLow = (~splitPrefixMaskHigh);

static const Word prefixMask =
(((uint)prefixLength==(uint)maxBasesPerWord)
 ? 0
 : (((Word)~0)<< boost::static_signed_min<numBitsPerWord - 1, prefixLength*numBitsPerBase>::value) );
static const Word suffixMask =
(((uint)suffixLength==(uint)maxBasesPerWord)
 ? 0
 : (((Word)~0)<< boost::static_signed_min<numBitsPerWord - 1, suffixLength*numBitsPerBase>::value) );

};

static const unsigned int maxOligoNum = 67108864;

typedef ushort FragmentErrorType;

#ifdef OLD_FRAGMENT_ERROR_DEFS
const int errorBits(6);
const FragmentErrorType incrementDist
( ((FragmentErrorType)1) << (2*errorBits) );
const FragmentErrorType errorMask1( (((FragmentErrorType)1)<<errorBits)-1 );
const FragmentErrorType errorMask2( errorMask1<<errorBits );
const FragmentErrorType distanceMask(~(errorMask1|errorMask2));
const FragmentErrorType moreThanTwoErrors
( ((FragmentErrorType)3) << (2*errorBits) );
const FragmentErrorType twoErrors
( ((FragmentErrorType)2) << (2*errorBits) );
const FragmentErrorType oneError
( ((FragmentErrorType)1) << (2*errorBits) );
#endif

// Constants to facilitate storing two error descriptions in an
// unsigned short.
//
// Storage scheme is (bits numbered starting at from least significant):
// Bits 1-6: Position of 1st error
// (so this storage scheme only works for <=2^6=64 base reads)
// Bits 7-8: XOR result of 1st error
// (can combine this with original base to infer what the error was)
// Bits 9-14: Position of 2nd error (if present, else zero)
// Bits 15-16: XOR of 2nd error (if present, else zero)
//
// If more than 2 errors, all bits are set to 1
// If more than 1 error, value of FragmentErrorType is greater than 255
// If no errors at all, all bits are off!

const int errorPosBits__(6);
const int errorTypeBits__(2);
const int errorBits__(errorPosBits__+errorTypeBits__);
const FragmentErrorType moreThanTwoErrors__(~((FragmentErrorType)0));
const FragmentErrorType errorPosMask1__(0x3F);  // binary 00111111
const FragmentErrorType errorTypeMask1__(0xC0); // binary 11000000
const FragmentErrorType errorInfoMask1__(errorPosMask1__|errorTypeMask1__);

const FragmentErrorType errorPosMask2__(errorPosMask1__<<errorBits__);
const FragmentErrorType errorTypeMask2__(errorTypeMask1__<<errorBits__);
const FragmentErrorType errorInfoMask2__(errorPosMask2__|errorTypeMask2__);
const FragmentErrorType oneError__(errorPosMask1__|errorTypeMask1__);
const FragmentErrorType noErrors__(0);

// The way in which an oligo matches to the genome data is described by
// one instance of MatchPosition
// one instance of MatchDescriptor
// The match as described by these two quantities may be categorised as
// being in one of seven states:
//
// 1: NM - No Match for this oligo
// 2: UE - Unique Exact match for this oligo
// 3: U1 - Unique match for this oligo with a single base substitution error
// 4: U2 - Unique match for this oligo, having 2 base substitution errors
// 5: RE - Unique Exact match for this oligo
// 6: R1 - Unique match for this oligo with a single base substitution error
// 7: R2 - Unique match for this oligo, having 2 base substitution errors


// MatchPosition stores a position in the genome (chromosome number + position)
// as a single unsigned int
// Capacity of an unsigned int is split into blocks according to the
// value of the most significant byte
// 0 : no match
// 1 - 239 : reserved for sequence data
// 240 - 255 : repeat
// blockRepeat + x : this oligo is the same as oligo number x in the batch
// ~0 (all bits 1) : match not attempted - oligo failed quality control (too many Ns etc)


typedef uint MatchPosition;

const MatchPosition noMatch(0);
const int blockShift(24);
const MatchPosition blockMask(((MatchPosition)0xFF)<<blockShift);
const MatchPosition blockPositionMask(((MatchPosition)-1)>>(8 * sizeof(MatchPosition) - blockShift));
const MatchPosition blockRepeat(((MatchPosition)0xF0)<<blockShift);
const MatchPosition blockSize(((MatchPosition)1)<<blockShift);
const MatchPosition qualityFailed(~(MatchPosition)0);
const MatchPosition repeatMasked(qualityFailed-1);

typedef uint OligoNumber;

const OligoNumber isReverseOligo
( ((OligoNumber)1) << ((8*sizeof(OligoNumber))-1) );
//( ((OligoNumber)1) << (MAX_HASH_BITS-1) );


const OligoNumber seed_bits[4] = { (((OligoNumber)0) << ((8*sizeof(OligoNumber))-3)),
                                   (((OligoNumber)1) << ((8*sizeof(OligoNumber))-3)),
                                   (((OligoNumber)2) << ((8*sizeof(OligoNumber))-3)),
                                   (((OligoNumber)3) << ((8*sizeof(OligoNumber))-3)) };



const int maxNumOligos(((uint)1<<((8*sizeof(OligoNumber))-1))-1);


// These are equivalent uchar definitions for pulling the error info
// out of a uchar in a MatchDescriptor prior to printing
const uchar ucharErrorPosMask(0x3F);
const uchar ucharErrorTypeMask(0xC0);

// These codes are used to store information about the Ns used in a match
const Word noNonIncorps(0);
// Any Ns that are in the oligo (possibly none) are all non-detections

const Word nonIncorpFirstN(0x55555555);
// this is 101010... in binary
// first N is a non-incorporation, second (if present) is a non-detection

const Word nonIncorpSecondN(0xAAAAAAAA);
// this is 010101... in binary
// second N is a non-incorporation, first is a non-detection

const Word nonIncorpBothNs(0xFFFFFFFF);
// this is 111111... in binary
// both Ns are non-incorporations

#define ALIGN_DP_BAND 10

} // namespace eland_ms
} // namespace casava

#endif // #ifndef CASAVA_ELAND_MS_ELAND_CONSTANTS_H
