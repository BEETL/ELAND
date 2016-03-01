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
 ** \file eland_ms/SeedMatch.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 ** Contains interface functionality from ELAND - specifically, stuff to do
 ** with scoring and storage of alignments
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_SEED_MATCH_H
#define CASAVA_ELAND_MS_SEED_MATCH_H

#include "ElandConstants.hh"

namespace casava
{
namespace eland_ms
{

struct SeedMatch
{
  SeedMatch( MatchPosition pos, uchar errors, bool reverse, uchar seed ) :
    pos_(pos), seeds_(0), errors_(errors), lastSeed_(seed), reverse_(reverse) {}
  MatchPosition pos_; // Position of seed 0 before extension. (this does not mean the match was created by seed 0)

  uchar seeds_   :3; //number of seeds that extended the match. 0 if it was not extended
  uchar errors_  :2; //number of errors. Valid numbers are 0,1,2
  uchar lastSeed_:2; //seed that matched. Last one in case the initial seed got extended.
                     //Valid are: (0-3 correspond to the multiseed seeds)
  uchar reverse_ :1; //1 indicates a reverse match

  // the operator is overloaded differently than MultiMatch; this one
  // is only used for sorting the matches according to their
  // seed/error ratio
  bool operator<( const SeedMatch& b ) const
  {
    return ( seeds_ == b.seeds_ ) ?  (errors_ > b.errors_ ) : ( seeds_ < b.seeds_ );
  }

} __attribute__ ((packed)); // ~struct SeedMatch

BOOST_STATIC_ASSERT(sizeof(SeedMatch)==sizeof(MatchPosition) + sizeof(uchar)); //only specializations are allowed to instantiate

} // namespace eland_ms
} // namespace casava



#endif // CASAVA_ELAND_MS_SEED_MATCH_H
