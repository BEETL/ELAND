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
 ** \file eland_ms/StatMachine.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 ** Contains interface functionality from ELAND - specifically, stuff to do
 ** with scoring and storage of alignments
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_STATE_MACHINE_H
#define CASAVA_ELAND_MS_STATE_MACHINE_H

#include "SeedMatch.hh"
#include "MatchPositionTranslator.hh"
#include "MatchDescriptor.hh"
#include "MultiMatch.hh"

namespace casava
{
namespace eland_ms
{

class StateMachine {

 public:
  StateMachine() {} // constructor
  ~StateMachine() {} // destructor

  // resize the state table; every read is associated with a certain
  // state according to the number of seeds that match
  void initialize( const OligoNumber&  s, const vector<int>& seedOffsets )
      { states_.resize( s ); matchType_.resize( s ); seedOffsets_ = seedOffsets; }

  // update the status of a read
  bool updateStatus( const OligoNumber& oligo );

  // Insert pass hit, if routine returns true, then call fetchHits, do
  // nothing otherwise
  bool insertSeedHit( const MatchPosition& thisPos, const MatchPosition& thisCode, int& update_cnt, MatchPositionTranslator& getMatchPos );

  // Extract the hits for a specific oligo
  list<MultiMatch> getHits( const OligoNumber& oligo, int maxItems, int seedNumberBase );

  // stitching together the reads we built our own matchType_ array
  // incorporate now the counts from MatchTableMulti::matchType_
  bool combineMatchDescriptor( const OligoNumber& oligo,const MatchDescriptor& descriptor );

  // how often did we query the list of positions
  void outputAccessStatistics(void);

  void clear( const OligoNumber& oligo );

  void clear( void );

  // keep track of your own MatchDescriptro vector
  vector< MatchDescriptor > matchType_;

 private:
  vector< vector<SeedMatch> > states_;
  vector<int> seedOffsets_;

  bool compare_matches_( const SeedMatch& a,const SeedMatch& b );

  //  vector<int> no_of_inserted_seeds_;
  //  vector<int> no_of_accesses_;

};


} // namespace eland_ms
} // namespace casava



#endif // CASAVA_ELAND_MS_STATE_MACHINE_H
