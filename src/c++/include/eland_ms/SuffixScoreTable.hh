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
 ** \file eland_ms/SuffixScoreTable.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 ** Contains interface functionality from ELAND - specifically, stuff to do
 ** with scoring and storage of alignments
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_SUFFIX_SCORE_TABLE_H
#define CASAVA_ELAND_MS_SUFFIX_SCORE_TABLE_H

#include "ElandConstants.hh"

namespace casava
{
namespace eland_ms
{

// class SuffixScoreTable
// This stores the scoring information for fast inexact matching of suffixes
// suffixScore_ - gives match score based on XOR of two suffixes
// quickByte_ - filters out words based on examining single bytes
struct SuffixScoreTable
{
SuffixScoreTable(const int lenA,
                 const int lenB,
                 const int lenC,
                 const int lenD)
  {
    int offset(0);

    buildScoreTable( scoreFragD_, lenD, offset );
    offset+=lenD;
    buildScoreTable( scoreFragC_, lenC, offset );
    offset+=lenC;
    buildScoreTable( scoreFragB_, lenB, offset );
    offset+=lenB;
    buildScoreTable( scoreFragA_, lenA, offset );

  } // ~SuffixScoreTable::SuffixScoreTable( int numBits )

  //  typedef bool ByteTable[numPossibleChars];
  //  typedef ByteTable ByteTableSet[numPossibleChars];

  void buildScoreTable( vector<FragmentErrorType>& score, const int fragLength, const int offset );

  //  vector<SuffixErrorType> suffixScore_;
  vector<FragmentErrorType> scoreFragA_;
  vector<FragmentErrorType> scoreFragB_;
  vector<FragmentErrorType> scoreFragC_;
  vector<FragmentErrorType> scoreFragD_;

  //  ByteTableSet quickByte_;
};

} // namespace eland_ms
} // namespace casava



#endif // CASAVA_ELAND_MS_SUFFIX_SCORE_TABLE_H
