/*
PROJECT: IMPALA (Inexact Matching Program ALlowing Ambiguity)
MODULE: GlobalUtilities.h
AUTHOR: A. J. Cox

 * Copyright (c) 2003-2006 Solexa
 *
 ** This software is covered by the "Illumina Genome Analyzer Software
 ** License Agreement" and the "Illumina Source Code License Agreement",
 ** and certain third party copyright/licenses, and any user of this
 ** source file is bound by the terms therein (see accompanying files
 ** Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
 ** Illumina_Source_Code_License_Agreement.pdf and third party
 ** copyright/license notices).
 */

#ifndef CASAVA_ALIGNMENT_GLOBAL_UTILITIES_HH
#define CASAVA_ALIGNMENT_GLOBAL_UTILITIES_HH

#include <fstream>
#include <sys/resource.h>
#include <sys/time.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/ptr_container/ptr_list.hpp>

// inner_product requires numeric
#include <numeric>

#include "common/Sequence.hh"
#include "alignment/BclReader.hh"

// Identifier to embed in source and binary files

// GlobalUtilities.h
// Header file for commonly used functions and definitions

using namespace std;
namespace fs = boost::filesystem;
namespace cc = casava::common;

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short ushort;

static const char baseNames [] = "ACGT";

// Word is used to store 2 bit-per-base encoded DNA data
typedef unsigned int Word;

enum
{
  numDifferentBases = 4,
  numBitsPerByte = 8,
  numPossibleChars = 1<<numBitsPerByte,
  numBitsPerBase = 2,
  //  mapChunkSize = 16777216,
  mapChunkSize = 1<<28, // 268435456
  numBitsPerWord = sizeof(Word)*numBitsPerByte,
  maxBasesPerWord = numBitsPerWord/numBitsPerBase,
  lineLength = 320,
  // maxSeqSize is the maximum length of sequence, used to create
  // buffers of appropriate size. NB IMPALA has its own (smaller) internal
  // sequence limit maxOligoLength, determined by the bit capacity of
  // unsigned integers.
  maxSeqSize = 256,
  // maxLineLength is the maximum number of characters assumed to be seen
  // in a line of an ASCII file. Used for buffer sizing
  maxLineLength= 8192
}; // ~enum


typedef Word TranslationTable[numPossibleChars];
typedef uchar TranslationTableChar[numPossibleChars];

static const Word nv(0xFF);

// whichBase is used to translate ASCII chars into 2 bit encodings
// 0 = A
// 1 = C
// 2 = G
// 3 = T
// AC 9.5.3 - 'N' was encoded as an A but this causes squashGenome to
// treat all bases as valid! Hence changing back to nv. Need to modify
// 'N' handling bit of impala to work round this
// whichBase is used by IMPALA
static const TranslationTable whichBase =
{
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv, /* next is 'A' */
00,nv,01,nv,nv,nv,02,nv,nv,nv,nv,nv,nv,nv,nv, /* next is 'P' */
nv,nv,nv,nv,03,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv, /* next is 'a' */
00,nv,01,nv,nv,nv,02,nv,nv,nv,nv,nv,nv,nv,nv, /* next is 'p' */
nv,nv,nv,nv,03,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,
nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv,nv
};


// nc and the following TranslationTableChar
// definitions are used by phageAlign

static const unsigned char nc('?');
//typedef unsigned char TranslationTable[256];

static const TranslationTableChar reverseCharASCII =
{
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,'.',nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc, /* next is 'A' */
'T',nc,'G',nc,nc,nc,'C',nc,nc,nc,nc,nc,nc,'N',nc, /* next is 'P' */
nc,nc,nc,nc,'A',nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc, /* next is 'a' */
'T',nc,'G',nc,nc,nc,'C',nc,nc,nc,nc,nc,nc,'N',nc, /* next is 'p' */
nc,nc,nc,nc,'A',nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc
};

// Map the real bases (ACGT) to 0-3 and Nn (plus legacy `.') to 4.
static const TranslationTable baseCodes =
{
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,04,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc, /* next is 'A' */
00,nc,01,nc,nc,nc,02,nc,nc,nc,nc,nc,nc,04,nc, /* next is 'P' */
nc,nc,nc,nc,03,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc, /* next is 'a' */
00,nc,01,nc,nc,nc,02,nc,nc,nc,nc,nc,nc,04,nc, /* next is 'p' */
nc,nc,nc,nc,03,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc
};

// Just ACGT.
static const TranslationTable realBaseCodes =
{
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc, /* next is 'A' */
00,nc,01,nc,nc,nc,02,nc,nc,nc,nc,nc,nc,nc,nc, /* next is 'P' */
nc,nc,nc,nc,03,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc, /* next is 'a' */
00,nc,01,nc,nc,nc,02,nc,nc,nc,nc,nc,nc,nc,nc, /* next is 'p' */
nc,nc,nc,nc,03,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,
nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc,nc
};

static const char baseChars[]="ACGT.";
//static const char* baseChars="ACGT.";

void resourceUsage( void );
void printWord( Word w, int l );
//char reverseComp( const char c);
//const char* getPrefixString( Word w, uint prefixLength );

// isBlank: return true if the character corresponds to
// a blank or uncalled base, i.e. 'n', 'N' or '.'
inline bool isBlank( const char c)
{
  return ((c=='.')||(c=='n')||(c=='N'));
} // ~isBlank


// countHeadNs: count the number of heading Ns in an oligo
//
inline
uint
countHeadNs(const char* pOligo,
            const uint oligoLength) {
  uint i(0);
  for (; i<oligoLength; ++i) {
    if(not isBlank(pOligo[i])) break;
  }
  return i;
} // ~countHeadNs


// countNs: count the number of heading, trailing and internal Ns in
// an oligo - used by adjustMatchPos in ELAND_outer.cpp and in Unsquash.cpp
inline void countNs( const char* pOligo,
		     const int oligoLength,
		     uint& headSize,
		     uint& tailSize,
		     uint& numInternalNs )
{
  headSize=0;
  tailSize=0;
  numInternalNs=0;
  bool inHead(true);
  for (int i(0); i<oligoLength; ++i)
  {
    //      if ((pOligo[i]=='n')||(pOligo[i]=='N'))
    if (isBlank(pOligo[i]))
    {
      headSize+=(inHead==true);
      ++tailSize;
      ++numInternalNs;
    } // ~if
    else
    {
      inHead=false;
      tailSize=0;
    } // ~else
  } // ~for
  numInternalNs-=headSize;
  numInternalNs-=tailSize;
} // ~countNs

/*****************************************************************************/
// Class: MappedReadOnlyPointer
// Memory map data in from a file and access it as you would a pointer
template <class T>
class MappedReadOnlyPointer
{
 public:
  MappedReadOnlyPointer( void ) : p_(NULL) {}
  void map( const char* name, const uint numElems=0 )
  {
    unmap();
    int fd;
    if ((fd=open(name, O_RDONLY, S_IRUSR))==-1)
    {
      cerr << "Error: could not open file " << name << endl;
      exit (1);
    } // ~if
    size_ = lseek( fd, 0, SEEK_END );
    if ((numElems!=0)&&(size_!= numElems*sizeof(T)))
    {
      cerr << "Expecting file of size " << numElems*sizeof(T)
	   << " bytes, got file of size " << size_ << endl;
      exit (1);
    } // ~if
    // mmap returns MAP_FAILED (-1) and not NULL on error
    if ((p_=(T*)mmap(0, size_, PROT_READ, MAP_SHARED, fd, 0))==(T*)MAP_FAILED)
    {
      cerr << "Could not memory map file " << name << endl;
      exit (1);
    } // ~if
    close (fd);

  } // ~void map( const string& name, const uint numElems )

  void unmap( void )
  {
    if (!p_) return;
    munmap( (void*)p_, size_ );
    p_ = NULL;
  } // ~void unmap( void )


  const T& operator*() const
  {
    return *p_;
  } // ~const T& operator*() const
  const T* operator->() const
  {
    return p_;
  } // ~const T* operator->() const

  const T& operator[]( int i ) const
  {
    return p_[i];
  } // ~const T& operator[]( int i ) const

  const T* operator+( const int& i ) const
  {
    return p_+i;
  } // ~const T* operator+( const int& i ) const

 private:
  uint size_;
  T* p_;

}; // ~template <class T> class MappedReadOnlyPointer

/*****************************************************************************/
// Class Name : Timer
// Description: Maintains info on actual and processing time
class Timer
{

  public:
  Timer( void );
  ostream& print( ostream& os );
  // timeNow: returns current date and time as an ASCII string
  const char* timeNow( void ) const;

  private:
  int numStamps;
  rusage thisUsage_;
  rusage lastUsage_;
  timeval thisTime_;
  timeval lastTime_;
  //  timeb thisTime_;
  //  timeb lastTime_;

}; // Timer

ostream& operator<<( ostream& os, Timer& timer );

/*****************************************************************************/

// reverseChar is used to compute reverse complements
// value of a uchar u is used as an index into reverseChar
// the value of reverseChar at that point gives the reverse complement of u
// e.g. reverseChar[AAAA]=reverseChar[0]=0xFF=TTTT
static const TranslationTableChar reverseChar =
{
  0xff, 0xbf, 0x7f, 0x3f, 0xef, 0xaf, 0x6f, 0x2f,
  0xdf, 0x9f, 0x5f, 0x1f, 0xcf, 0x8f, 0x4f, 0xf,
  0xfb, 0xbb, 0x7b, 0x3b, 0xeb, 0xab, 0x6b, 0x2b,
  0xdb, 0x9b, 0x5b, 0x1b, 0xcb, 0x8b, 0x4b, 0xb,
  0xf7, 0xb7, 0x77, 0x37, 0xe7, 0xa7, 0x67, 0x27,
  0xd7, 0x97, 0x57, 0x17, 0xc7, 0x87, 0x47, 0x7,
  0xf3, 0xb3, 0x73, 0x33, 0xe3, 0xa3, 0x63, 0x23,
  0xd3, 0x93, 0x53, 0x13, 0xc3, 0x83, 0x43, 0x3,
  0xfe, 0xbe, 0x7e, 0x3e, 0xee, 0xae, 0x6e, 0x2e,
  0xde, 0x9e, 0x5e, 0x1e, 0xce, 0x8e, 0x4e, 0xe,
  0xfa, 0xba, 0x7a, 0x3a, 0xea, 0xaa, 0x6a, 0x2a,
  0xda, 0x9a, 0x5a, 0x1a, 0xca, 0x8a, 0x4a, 0xa,
  0xf6, 0xb6, 0x76, 0x36, 0xe6, 0xa6, 0x66, 0x26,
  0xd6, 0x96, 0x56, 0x16, 0xc6, 0x86, 0x46, 0x6,
  0xf2, 0xb2, 0x72, 0x32, 0xe2, 0xa2, 0x62, 0x22,
  0xd2, 0x92, 0x52, 0x12, 0xc2, 0x82, 0x42, 0x2,
  0xfd, 0xbd, 0x7d, 0x3d, 0xed, 0xad, 0x6d, 0x2d,
  0xdd, 0x9d, 0x5d, 0x1d, 0xcd, 0x8d, 0x4d, 0xd,
  0xf9, 0xb9, 0x79, 0x39, 0xe9, 0xa9, 0x69, 0x29,
  0xd9, 0x99, 0x59, 0x19, 0xc9, 0x89, 0x49, 0x9,
  0xf5, 0xb5, 0x75, 0x35, 0xe5, 0xa5, 0x65, 0x25,
  0xd5, 0x95, 0x55, 0x15, 0xc5, 0x85, 0x45, 0x5,
  0xf1, 0xb1, 0x71, 0x31, 0xe1, 0xa1, 0x61, 0x21,
  0xd1, 0x91, 0x51, 0x11, 0xc1, 0x81, 0x41, 0x1,
  0xfc, 0xbc, 0x7c, 0x3c, 0xec, 0xac, 0x6c, 0x2c,
  0xdc, 0x9c, 0x5c, 0x1c, 0xcc, 0x8c, 0x4c, 0xc,
  0xf8, 0xb8, 0x78, 0x38, 0xe8, 0xa8, 0x68, 0x28,
  0xd8, 0x98, 0x58, 0x18, 0xc8, 0x88, 0x48, 0x8,
  0xf4, 0xb4, 0x74, 0x34, 0xe4, 0xa4, 0x64, 0x24,
  0xd4, 0x94, 0x54, 0x14, 0xc4, 0x84, 0x44, 0x4,
  0xf0, 0xb0, 0x70, 0x30, 0xe0, 0xa0, 0x60, 0x20,
  0xd0, 0x90, 0x50, 0x10, 0xc0, 0x80, 0x40, 0x0
}; // ~reverseChar

union Oligo
{
  Oligo() { ui[0]=0;ui[1]=0;}
  Word ui[2];
  uchar uc[2*sizeof(Word)];
  bool operator==( const Oligo& rhs ) const
  {
    return( (ui[0]==rhs.ui[0])&&(ui[1]==rhs.ui[1]));
  } // ~bool operator==( const Oligo& rhs )

  bool operator<( const Oligo& rhs ) const
  {
    return( (ui[1]<rhs.ui[1])
	    ||((ui[1]==rhs.ui[1])&&(ui[0]<rhs.ui[0])) );
  } // ~bool operator==( const Oligo& rhs )

  bool operator!=( const Oligo& rhs ) const
  {
    return( (ui[0]!=rhs.ui[0])||(ui[1]!=rhs.ui[1]));
  } // ~bool operator!=( const Oligo& rhs )

}; // ~Oligo

/*****************************************************************************/
// ExpandedTranslationTable
// Convert ASCII sequence file into 2 bits per base 4 bits at a time

class ExpandedTranslationTable
{
public:
  ExpandedTranslationTable(int oligoLength);
  void translate( const char* buf, Oligo& ol, Oligo& rc ) const;
  Word t_[numPossibleChars*numPossibleChars];
private:
  const int prefixLength_;
  const bool isOddLength_;
  const int reverseShift1_;
  const int reverseShift2_;
}; // ~class ExpandedTranslationTable

typedef short BaseScore;

/*****************************************************************************/
// ScoreSource: obtain score that should be assigned to a
// given sequenced base at a given cycle
class ScoreSource
{
 public:
  virtual ~ScoreSource() {}
  // Return score information: score to give a genomic base
  // when aligned against cycle-th base
  virtual BaseScore getScore( char base, uint cycle ) const=0;

  virtual BaseScore getScore( uint baseNum, uint cycle ) const=0;

}; // ~class ScoreSource

/*****************************************************************************/
// ScoreSourceFilter:
// wrap around another instance of ScoreSource and filter bases
// Memory management policy: is not responsible for deleting pRaw_.
// Based on OligoSourceFilter
class ScoreSourceFilter : public ScoreSource
{
 public:
  ScoreSourceFilter( const char* basesToUse ) :
    pRaw_(NULL)
  {
    for (uint i(0);i<strlen(basesToUse);i++)
    {
      if (toupper(basesToUse[i])=='Y')
        baseIndex_.push_back(i);
    } // ~for
  } // ~ctor

  virtual ~ScoreSourceFilter() {}

  // Link a raw data source to filter, deleting previous if necessary
  void link( ScoreSource* pRaw )
  {
    //    if (pRaw_!=NULL) delete pRaw_;
    pRaw_=pRaw;
  } // ~link


  // Return score information: score to give a genomic base
  // when aligned against cycle-th base
  virtual BaseScore getScore( char base, uint cycle ) const
  {
    return pRaw_->getScore(base, baseIndex_[cycle] );
  } // ~getScore

  virtual BaseScore getScore( uint baseNum, uint cycle ) const
  {
    return pRaw_->getScore(baseNum, baseIndex_[cycle] );
  } // ~getScore
 protected:
  ScoreSource* pRaw_;
  vector<uint> baseIndex_;

  }; // ~class ScoreSourceFilter


/*****************************************************************************/
// OligoSource: obtain oligos in ASCII format from some source
// (e.g. fasta file, raw sequence file, database, ... )
class OligoSource
{
 public:

// getOligoSource: given the name of a file of oligo data, returns a pointer
// to an instance of the appropriate subclass of OligoSource.
  static OligoSource* getOligoSource( const char* fileName );

  OligoSource()
      : isNoMask_(true) {}

  virtual ~OligoSource() {}

  // Returns reference to next Sequence (supersedes getNextOligo).
  // isValid will be false if there are no sequences left.
  virtual const casava::common::Sequence& getNextSequenceSelect(bool& isValid,
                                                                const bool isProvideHeader,
                                                                const bool isProvideQualities) = 0;

  const casava::common::Sequence& getNextSequence(bool& isValid) {
      return getNextSequenceSelect(isValid,true,true);
  }

  // Returns reference to last Sequence fetched (supersedes getLastOligo).
  // isValid will be false if no sequences have been successfully read.
  virtual const casava::common::Sequence& getLastSequence(bool& isValid) const = 0;

  // Returns pointer to ASCII sequence of next oligo, or null if at end
  virtual const char* getNextOligo( void ) =0;

  const char* getNextOligoSelect(const bool isProvideHeader,
                                 const bool isProvideQualities)
  {   
      bool isValid(false);
      const casava::common::Sequence& seq(getNextSequenceSelect(isValid,
                                                                isProvideHeader,
                                                                isProvideQualities));
      return (isValid ? seq.getData().c_str() : NULL);
  }

  // Returns pointer to ASCII sequence of last oligo fetched, or null at end
  virtual const char* getLastOligo( void ) const=0;


  // Returns pointer to ASCII name of last oligo read
  virtual const char* getLastName( void )
  {
    return NULL;
  } // ~getLastName

  // Rewind - next oligo read will be first in list
  virtual void rewind ( void ) =0;

  void setMask( const vector<bool>& mask ){
      isNoMask_ = false;
      mask_ = mask;
  }

  void unSetMask() {
      std::vector<bool> tmpMask;
      std::swap(tmpMask,mask_);
      tmpMask.clear();
      
      isNoMask_ = true;
  }

  // if we set mask_, then we have to know how many sequences we
  // skipped such that we have correct values for oligoNum in
  // buildTable, for instance.
  virtual int getNoSkippedSequences( void ) { return 1; }

protected:
  bool isNoMask_;
  std::vector<bool> mask_;
}; // ~class OligoSource


/*****************************************************************************/
// OligoSourceFilter
// wrap around another instance of OligoSource and filter bases
// Memory management policy: is not responsible for deleting pRaw_.
class OligoSourceFilter : public OligoSource
{
 public:
  OligoSourceFilter( const char* basesToUse )
      : pRaw_(NULL), sequenceIsValid_(false)
  {
    for (uint i(0);i<strlen(basesToUse);i++)
    {
      if (toupper(basesToUse[i])=='Y')
        baseIndex_.push_back(i);
    } // ~for
  } // ~ctor

  virtual ~OligoSourceFilter() {}

  // Link a raw data source to filter, deleting previous if necessary
  void link( OligoSource* pRaw )
  {
    //    if (pRaw_!=NULL) delete pRaw_;
    pRaw_=pRaw;
    sequenceIsValid_ = false;
  } // ~link

  // Returns reference to next Sequence (supersedes getNextOligo).
  // isValid will be false if there are no sequences left.
  virtual const casava::common::Sequence& getNextSequenceSelect(bool& isValid,
                                                                const bool isProvideHeader,
                                                                const bool isProvideQualities)
  {
      assert(pRaw_ != 0);
      maskedSequence_ = pRaw_->getNextSequenceSelect(sequenceIsValid_,
                                                     isProvideHeader,
                                                     isProvideQualities);
      isValid = sequenceIsValid_;

      if (isValid) {
          maskedSequence_.mask(baseIndex_);
      }

      return maskedSequence_;
  }

  // Returns reference to last Sequence fetched (supersedes getLastOligo).
  // isValid will be false if there are no sequences left.
  virtual const casava::common::Sequence& getLastSequence(bool& isValid) const
  {
      isValid = sequenceIsValid_;
      return maskedSequence_;
  }


  // Returns pointer to ASCII sequence of next oligo, or null if at end
  virtual const char* getNextOligo( void )
  {
      bool unusedBool(false);
      (void) getNextSequence(unusedBool);
      return (sequenceIsValid_ ? maskedSequence_.getData().c_str() : NULL);
  } // ~getNextOligo


  // Returns pointer to ASCII sequence of last oligo fetched
  virtual const char* getLastOligo( void ) const
  {
      return (sequenceIsValid_ ? maskedSequence_.getData().c_str() : NULL);
  } // ~getLastOligo


  // Returns pointer to ASCII name of last oligo read
  virtual const char* getLastName( void )
  {
    return pRaw_->getLastName();
  } // ~getLastName

  // Rewind - next oligo read will be first in list
  virtual void rewind ( void )
  {
    pRaw_->rewind();
  } // ~rewind

 protected:
  OligoSource* pRaw_;
  vector<uint> baseIndex_;
  casava::common::Sequence maskedSequence_;
  bool sequenceIsValid_;
}; // ~class OligoSourceFilter


/*****************************************************************************/
// OligoSourceFile: read oligos from a file
class OligoSourceFile : public OligoSource
{
 public:

  OligoSourceFile( const char* oligoFileName ) :
      pFile_(fopen(oligoFileName, "r"))
  {
    if (pFile_==NULL)
    {
      cerr << "Error in OligoSourceFile: could not open file "
	   << oligoFileName << endl;
      exit (1);
    } // ~if
  } // ~ctor

  ~OligoSourceFile()
  {
    fclose(pFile_);
  } // ~dtor

  // Returns reference to last Sequence fetched (supersedes getLastOligo).
  // isValid will be false if there are no sequences left.
  virtual const casava::common::Sequence& getLastSequence(bool& isValid) const
  {
      isValid = sequenceIsValid_;
      return sequence_;
  }

  // Returns pointer to ASCII sequence of next oligo, or null if at end
  virtual const char* getNextOligo( void )
  {
      bool unusedBool(false);
      (void) getNextSequence(unusedBool);
      return (sequenceIsValid_ ? sequence_.getData().c_str() : NULL);
  } // ~getNextOligo

  // Returns pointer to ASCII sequence of last oligo fetched
  virtual const char* getLastOligo( void ) const
  {
      return (sequenceIsValid_ ? sequence_.getData().c_str() : NULL);
  }

// getOligoSourceFile: given the name of a file of oligo data, returns a
// pointer to an instance of the appropriate subclass of OligoSourceFile.
// If can read a Sequence : OligoSourceGoat
// Otherwise : Examines first character of file:-
// If a '>':
// assumes a fasta file and returns an OligoSourceFasta
// If valid sequence character or blank:
// assumes raw sequence and  returns OligoSourceRaw
// If space or dash :
// assumes quality value format oligo file and returns OligoSourceScore
  static OligoSourceFile* getOligoSourceFile( const char* fileName );

  virtual void rewind( void )
  {
    curSeq_ = 1;
    ::rewind(pFile_);
  } // ~rewind

  virtual int getNoSkippedSequences( void ) { return skippedSequences_; }

 protected:
  FILE* pFile_;

  int curSeq_;
  int skippedSequences_;
  casava::common::Sequence sequence_;
  bool sequenceIsValid_;

  char nameBuf_[maxLineLength];
}; // ~class OligoSourceFile


/*****************************************************************************/
// Read oligos from a raw sequence file, i.e. an ASCII file with one oligo
// per line. NB performs no length checking on the oligos
class OligoSourceRaw : public OligoSourceFile
{
 public:
  OligoSourceRaw( const char* oligoFileName ) :
    OligoSourceFile(oligoFileName), oligoNum_(0)
    {
      strcpy(nameBuf_,oligoFileName);
      pNum_=strchr(nameBuf_,'\0');
      assert(pNum_!=NULL);
      *pNum_++='-';
    } // ~ctor

  // Returns reference to next Sequence (supersedes getNextOligo).
  // isValid will be false if there are no sequences left.
  virtual const casava::common::Sequence& getNextSequenceSelect(bool& isValid,
                                                                const bool,
                                                                const bool)
  {
      ++oligoNum_;
      isValid = sequenceIsValid_ = (fscanf(pFile_, "%s\n", oligoBuf_) != EOF);

      if (sequenceIsValid_) {
          sequence_.setData(oligoBuf_);
      }

      return sequence_;
  }

  virtual const char* getLastName( void )
  {
    sprintf( pNum_, "%d", oligoNum_ );
    return nameBuf_;
  } // ~getLastName( void )

  virtual void rewind( void )
  {
    oligoNum_=0;
    OligoSourceFile::rewind();
  } // ~rewind


 protected:
  char oligoBuf_[maxLineLength];
  int oligoNum_;
  char* pNum_;

}; // ~class OligoSourceRaw


/*****************************************************************************/
// Read oligos from a raw sequence file, i.e. an ASCII file with one oligo
// per line. NB performs no length checking on the oligos
class OligoSourceFasta : public OligoSourceFile
{
 public:

  OligoSourceFasta( const char* oligoFileName );

  // Returns reference to next Sequence (supersedes getNextOligo).
  // isValid will be false if there are no sequences left.
  virtual const casava::common::Sequence& getNextSequenceSelect(bool& isValid,
                                                                const bool,
                                                                const bool)
  {
    // reset the skipped sequences counter to 0
    skippedSequences_ = 0;
    //    cerr << "curSeq = " << curSeq_ << "\t" << "mask_[" << curSeq_ << "] = " << mask_[curSeq_] << endl;

    while(1)
    {
        sequenceIsValid_ = (fscanf(pFile_, "%s\n", nameBuf_) != EOF);

        if (sequenceIsValid_) {
            sequenceIsValid_ = (fscanf(pFile_, "%s\n", oligoBuf_) != EOF);


            // loop until we found the next available sequence
            if (sequenceIsValid_) {
                curSeq_++;
                if( isNoMask_ || mask_[curSeq_-1])
                {
                    sequence_.setData(oligoBuf_);

                    break;
                }
                else
                {
                    skippedSequences_++;
                }
            } else {
                cerr << "Error in OligoSourceFasta: failed to read sequence"
                << endl;
                // Existing behaviour.
                exit(1);
            }
        } else {
            cerr << "Error in OligoSourceFasta: found name but no sequence"
            << endl;

            break;
        }
    }

    isValid = sequenceIsValid_;
    return sequence_;
  }

  // Returns pointer to ASCII name of last oligo read
  virtual const char* getLastName( void )
  {
    return nameBuf_;
  } // ~getLastName( void )

  virtual int getNoSkippedSequences( void ) { return skippedSequences_; }

protected:
  char oligoBuf_[maxLineLength];

}; // ~class OligoSourceFasta : public OligoSourceFile


/*****************************************************************************/
// Read oligos from a Goat format file
// containing sequence interspersed with comment lines that start with a #
class OligoSourceGoat : public OligoSourceFile
{
 public:

  OligoSourceGoat( const char* oligoFileName );
  ~OligoSourceGoat();

  // Returns reference to next Sequence (supersedes getNextOligo).
  // isValid will be false if there are no sequences left.
  virtual const casava::common::Sequence& getNextSequenceSelect(bool& isValid,
                                                                const bool,
                                                                const bool);

  // Returns pointer to ASCII name of last oligo read
  virtual const char* getLastName( void );

  // Rewind - next oligo read will be first in list
  virtual void rewind ( void );

protected:
  std::ifstream qseq_file_;
}; // ~class OligoSourceGoat : public OligoSourceFile


/*****************************************************************************/
// Read oligos from a set of quality values. Inherits from OligoSourceRaw
// because
class OligoSourceScore: public OligoSourceRaw, public ScoreSource
{
 public:
  OligoSourceScore( const char* oligoFileName ) :
    OligoSourceRaw(oligoFileName) {}

  // Functions inherited from OligoSourceFile

  // Returns reference to next Sequence (supersedes getNextOligo).
  // isValid will be false if there are no sequences left
  virtual const casava::common::Sequence& getNextSequenceSelect(bool& isValid,
                                                                const bool isProvideHeader,
                                                                const bool isProvideQualities);

  // Functions inherited from ScoreSource

  // Returns quality info
  virtual BaseScore getScore( char base, uint cycle ) const;
  virtual BaseScore getScore( uint baseNum, uint cycle ) const;

 private:
  BaseScore scoreTable_[numDifferentBases][maxSeqSize];

}; // ~class OligoSourceScore


/*****************************************************************************/
// Read oligos from a directory containing other OligoSources
class OligoSourceDirectory : public OligoSource
{
 public:
  OligoSourceDirectory( const char* dirName );

  ~OligoSourceDirectory() { delete pSource_; }

  // Returns reference to next Sequence (supersedes getNextOligo).
  // isValid will be false if there are no sequences left.
  virtual const casava::common::Sequence& getNextSequenceSelect(bool& isValid,
                                                                const bool isProvideHeader,
                                                                const bool isProvideQualities);

  // Returns reference to last Sequence fetched (supersedes getLastOligo).
  // isValid will be false if there are no sequences left.
  virtual const casava::common::Sequence& getLastSequence(bool& isValid) const;

  // Returns pointer to ASCII sequence of next oligo, or null if at end
  virtual const char* getNextOligo( void );

  // Returns pointer to ASCII sequence of last oligo fetched, null at end
  virtual const char* getLastOligo( void ) const
  {
    return pSource_->getLastOligo();
  } // ~getLastOligo( void )

  // Returns pointer to ASCII name of last oligo read
  virtual const char* getLastName( void );

  // Rewind - next oligo read will be first in list
  virtual void rewind ( void );

 private:
  vector<string> fileNames_;
  vector<string>::iterator pName_;
  OligoSource* pSource_;
  casava::common::Sequence dummySequence_;
}; // ~class OligoSourceDirectory

// TBD: OligoSourceMySQL, OligoSourceFromChromosome


/*****************************************************************************/
// class ScoreSourceBasic: simplest possible scoring system, constant
// scores for match, mismatch and blank
class ScoreSourceBasic : public ScoreSource
{
 public:
  ScoreSourceBasic( const OligoSource& oligos,
                    BaseScore scoreMatch=1,
                    BaseScore scoreBlank=0,
                    BaseScore scoreMismatch=-1 ) :
    oligos_(oligos),
    scoreMatch_(scoreMatch),
    scoreBlank_(scoreBlank),
    scoreMismatch_(scoreMismatch) {}

  // Returns quality information
  virtual BaseScore getScore( uint baseNum, uint cycle ) const;
  virtual BaseScore getScore( char base, uint cycle ) const;

 private:
  const OligoSource& oligos_;
  const BaseScore scoreMatch_;
  const BaseScore scoreBlank_;
  const BaseScore scoreMismatch_;
}; // ~ScoreSourceBasic


/*****************************************************************************/
// class ScoreSourceCycle: scoring system as output by score.pl - separate
// substitution matrix for each cycle
class ScoreSourceCycle: public ScoreSource
{
 public:
  ScoreSourceCycle( const OligoSource& oligos, const char* scoreFile );

  virtual BaseScore getScore( uint baseNum, uint cycle ) const;
  virtual BaseScore getScore( char base, uint cycle ) const;


 protected:
  const OligoSource& oligos_;
  BaseScore
    scoreTable_[maxSeqSize][numDifferentBases+1][numDifferentBases];
}; // ~ScoreSourceCycle


/*****************************************************************************/
// ValidRegion describes a region of valid bases in a squashed file,
// corresponding to a region containing only A, G, C or T (i.e. no Ns or other
// funny business) in the original sequence.
//
// bases in a squashed file are numbered starting at zero
// (although match positions are output numbered starting at 1)
// start gives the first valid base in a region (0=first base in file)
// finish gives the last valid base of the region (ie not the first invalid
// base after the end of the region)
struct ValidRegion
{
  ValidRegion() : start(0), finish(0){}
  ValidRegion(uint start, uint finish) : start(start), finish(finish){}
  uint start;
  uint finish;
}; // struct ValidRegion

// FileReader: memory maps a chromosome file into memory and gives out
// characters byte by byte, ignoring the first line and any '\n' characters
class FileReader
{
public:
  FileReader( const char* fileName );
  ~FileReader();

  //  MatchPosition scan( OligoHashTable& hashTable, const MatchPosition currentBlock );

  const ValidRegion* getFirstValid( void ) const
  {
    return pValid_;
  } // ~getFirstValid( void ) const
  const ValidRegion* getLastValid( void ) const
  {
    return pLastValid_;
  } // ~getFirstValid( void ) const
  const Word* getSeqStart( void ) const
  {
    return pStart_;
  } // ~getSeqStart( void ) const

  // getLastValidBase - position of last valid base, where
  // first base in file has base zero
  uint getLastValidBase( void ) const
  {
    return (pLastValid_-1)->finish;
  } // ~getLastBase;

private:
  //  const int prefixLength_;
  //  const int oligoLength_;
  //  const Word prefixMask_;
  //  Word wordMask_[maxBasesPerWord];

  string seqFileName_;
  string vldFileName_;

  int seqFileSize_;
  int vldFileSize_;

  void *pMapSeq_;
  void *pMapVld_;

  // pointer to first Word in file
  Word* pStart_;
  // pointer to end of file, ie first Word not in file
  Word* pEnd_;
  // pointer to current Word in current valid region
  Word* pWord_;
  // pointer to last Word in current valid region
  Word* pLastInRegion_;
  // pointer to first Word after end of current valid region
  Word* pEndOfRegion_;

  //  int nextBase_;
  //  int basesPerWord_;
  //  int basesPerWordLast_;


  ValidRegion* pValid_;
  ValidRegion* pLastValid_;


}; // ~class FileReader


// MirrorSorter - sort elements pBegin[left] to pBegin[right] *inclusive*,
// uses quicksort routine from pg.87 of Kernighan and Ritchie
// replaces mirrorSort and Swapper stuff in previous versions
// Change is to make interface neater (arguably) and allows sorting criterion
// to involve V and W. Need to provide a LessThan class, for 2 arrays of
// ints this might be:
//struct LessThan
//{
//  LessThan( int* p1, int* p2 ) : p1_(p1), p2_(p2) {}
//  bool operator()( int i, int j ) const
//  {
//    //    return true;
//    return ( (p1_[i]<p1_[j])
//	     || ( (p1_[i]==p1_[j])&&(p2_[i]<p2_[j])));
//  } // ~op<
//  const int* p1_;
//  const int* p2_;
//};

template <class V, class W>
class MirrorSorter
{
public:
  MirrorSorter( V* pv, W* pw ) : pBeginV_(pv), pBeginW_(pw) {}
  void swap( int left, int right )
  {
    v_ = pBeginV_[left];
    pBeginV_[left] = pBeginV_[right];
    pBeginV_[right] = v_;

    w_ = pBeginW_[left];
    pBeginW_[left] = pBeginW_[right];
    pBeginW_[right] = w_;
  } // ~void operator()(V* begin, V* left, V* right )

  template <class LessThan>
  void operator()( int left, int right, const LessThan& lessThan )
  {
    int i, last;

    if (left>=right) return;
    swap( left, (left+right)/2);
    last = left;
    for (i=left+1; i <= right; i++ )
      if (lessThan(i,left)) swap ( ++last, i );
    //      if (pBegin[i]<pBegin[left]) swap ( pBegin, ++last, i );
    swap( left, last );
    this->operator()( left, last-1, lessThan );
    this->operator()( last+1, right, lessThan );
  }

  // provides STL-style interface,
  // i.e. sort from *pStart up to but not including *pEnd
  template <class LessThan>
  void operator()( V* pStart, V* pEnd, const LessThan& lessThan )
  {
    this->operator()(pStart-pBeginV_, pEnd-pBeginV_-1, lessThan );
  }


private:
  V v_;
  W w_;
  V* pBeginV_;
  W* pBeginW_;
}; // ~template <class V, class SWAPPER> void mirrorSort



// compute the Hamming distance between two strings
class Hamming
{

public:
  Hamming(void) {} // c'tor
  ~Hamming(void) {} //d'tor

  int operator()(string const& s1, string const& s2)
  {
    int difference = std::inner_product(s1.begin(),s1.end(),
					s2.begin(),
					0,
					std::plus<int>(),diff_
					);
    return difference;
  }

private:
  inline static int diff_(char s1,char s2)
    {
        return s1 == s2?0:1;
    }

};




#ifdef XXX
// mirrorSort and Swapper: sort a vector of V and a vector of W by value of V.
// formerly in IndexGenome.cpp

// Example usage:
//  Swapper<Word, MatchInfoType> swapPos(pos_.begin());
//  mirrorSort(suf_.begin(),
//	 suf_.begin()+pArray_[0][j],
//	 suf_.begin()+pAccum_[0][j],
//	 swapPos);

// Swapper: swaps two elements of a vector of V and "silently" swaps the
// corresponding elements of a vector of W

template <class V, class W> class Swapper
{
public:
  Swapper( W* pBeginW ) : pBeginW_(pBeginW) {}

  void operator()(V* pBegin, int left, int right )
  {
    v_ = pBegin[left];
    pBegin[left] = pBegin[right];
    pBegin[right] = v_;

    w_ = pBeginW_[left];
    pBeginW_[left] = pBeginW_[right];
    pBeginW_[right] = w_;

  } // ~void operator()(V* begin, V* left, V* right )


private:
  V v_;
  W w_;
  W* pBeginW_;
}; // ~template <class V, class W> class Swapper

// Swapper3: swaps three elements of a vector of V and "silently" swaps the
// corresponding elements of two vectors W and X

template <class V, class W, class X> class Swapper3
{
public:
  Swapper3( W* pBeginW, X* pBeginX ) : pBeginW_(pBeginW), pBeginX_(pBeginX) {}

  void operator()(V* pBegin, int left, int right )
  {
    v_ = pBegin[left];
    pBegin[left] = pBegin[right];
    pBegin[right] = v_;

    w_ = pBeginW_[left];
    pBeginW_[left] = pBeginW_[right];
    pBeginW_[right] = w_;

    x_ = pBeginX_[left];
    pBeginX_[left] = pBeginX_[right];
    pBeginX_[right] = x_;

  } // ~void operator()(V* begin, V* left, V* right )


private:
  V v_;
  W w_;
  X x_;
  W* pBeginW_;
  X* pBeginX_;
}; // ~template <class V, class W, class X> class Swapper3

// mirrorSort - sort elements pBegin[left] to pBegin[right] inclusive,
// uses quicksort routine from pg.87 of Kernighan and Ritchie
template <class V, class SWAPPER> void mirrorSort
( V* pBegin, int left, int right, SWAPPER& swap )
{
  int i, last;

  if (left>=right) return;
  swap( pBegin, left, (left+right)/2);
  last = left;
  for (i=left+1; i <= right; i++ )
    if (pBegin[i]<pBegin[left]) swap ( pBegin, ++last, i );
  swap( pBegin, left, last );
  mirrorSort( pBegin, left, last-1, swap );
  mirrorSort( pBegin, last+1, right, swap );
} // ~template <class V, class SWAPPER> void mirrorSort

// mirrorSort - provides STL-style interface to mirrorSort,
// i.e. sort from *pLeft up to but not including *pRight
template <class V, class SWAPPER> void mirrorSort
( V* pBegin, V* pLeft, V* pRight, SWAPPER& swap )
{
  mirrorSort(pBegin, pLeft-pBegin, pRight-pBegin-1, swap);
} // ~template <class V, class SWAPPER> void mirrorSort



#endif

const vector<bool> expandUseBases(const std::string &useBases);

#endif //CASAVA_ALIGNMENT_GLOBAL_UTILITIES_HH
