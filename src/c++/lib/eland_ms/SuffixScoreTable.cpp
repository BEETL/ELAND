/**
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
 ** \file eland_ms/SuffixScoreTable.cpp
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#include "eland_ms/SuffixScoreTable.hh"

namespace casava
{
namespace eland_ms
{

void
SuffixScoreTable::buildScoreTable
( vector<FragmentErrorType>& score, const int fragLength, const int offset )
{
  const Word numFrags( 1<<(fragLength<<1) );
    const int start_test(((Word)0x3) << (numBitsPerBase*(fragLength-1)));
  Word test, result;

  FragmentErrorType thisScore;
  int numErrors;
  score.resize(numFrags);
  for ( Word i(0); i<numFrags; i++ )
  {
    test = start_test;
    thisScore=0;
    numErrors=0;
    for (int j(1); j<=fragLength; j++, test>>=numBitsPerBase )
    {
      if ((result=(i&test))!=0)
      {
        ++numErrors;

        if (numErrors<=2)
        {
          if (numErrors==2)
            thisScore<<=errorBits__;

          thisScore|= (j+offset);

          thisScore
            |= ((result>>(numBitsPerBase*(fragLength-j)))<<errorPosBits__);


          //            |= ((result>>(numBitsPerBase*(fragLength-1)))<<errorPosBits__);
          //	thisScore+=incrementDist;
          //	if ((thisScore&errorMask1)==0)
          //  thisScore+=j+offset;
          //	else if ((thisScore&errorMask2)==0)
          //	  thisScore+=((j+offset)<<errorBits);
        } // ~if
        else thisScore=moreThanTwoErrors__;

      } // ~if
    } // ~for j
    //    if (thisScore>=moreThanTwoErrors) thisScore &= distanceMask;
    score[i]=thisScore;
    //    printWord(i, fragLength);
    //   cout << " " << (thisScore&errorPosMask1__) << " "<< ((thisScore&errorPosMask2__)>>8) << endl;

    //        cout << " " << (thisScore>>(2*errorBits))
    // << " " << (thisScore&errorMask1)
    // << " " << ((thisScore&errorMask2)>>errorBits) << endl;
  } // ~for i

}

} //namespace eland_ms
} //namespace casava

