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
 ** \file eland_ms/ContigNameFinder.cpp
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/
/*
PROJECT: ELAND (Efficient Local Alignment of Nucleotide Data)
MODULE: ELAND_outer.cpp
AUTHOR: A. J. Cox

Copyright (c) 2003-2006 Solexa Limited. All rights reserved.
This source file is covered by the "Solexa Public Source License"
agreement and bound by the terms therein.
*/

#include "alignment/ELAND_unsquash.h"
#include "eland_ms/ContigNameFinder.hh"

namespace casava
{
namespace eland_ms
{


// getContigNameFinder: given the name of directory of squashed sequences
// and a sequence within it, return the appropriate contig name finder
ContigNameFinder* ContigNameFinder::getContigNameFinder
( const string& directoryName, const string& chromName )
{
  string indexName=directoryName+chromName+(string)".idx";
  cerr << "Will look for index file " << indexName << endl;
  FILE* pIndex=fopen(indexName.c_str(),"r");
  if (pIndex==NULL)
  {
    cerr << "... not found, will not index" << endl;
    return new ContigNameFinderNull();
  } // ~if
  else
  {
    cerr << "... found, creating index" << endl;
    return new ContigNameFinderIndex( pIndex );
  } // ~else
  //fclose (pIndex);
}

} //namespace eland_ms
} //namespace casava
