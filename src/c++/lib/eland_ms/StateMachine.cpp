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
 ** \file eland_ms/StateMachine.cpp
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#include "eland_ms/StateMachine.hh"

namespace casava
{
namespace eland_ms
{

// SEED_DEVIATION gives the number of indels that we allow per
// seed. We have have a 100bp read and SEED_DEVIATION is 2, then we
// allow 2 indels for seed 1 and 4 indels for seed 2
//
// if SEED_DEVIATION is set to 0, then we do not allow any deviation
#define SEED_DEVIATION 3


bool StateMachine::insertSeedHit( const MatchPosition& thisPos, const MatchPosition& thisCode, int& update_cnt,MatchPositionTranslator& getMatchPos)
{
  uint numErrors=(thisCode>>30)&0x3;
  uchar reverseFlag=(thisCode>>29)&0x1;
  uchar seedNo=(thisCode>>27)&0x3;
  uint thisOligo=thisCode&(((uint)~0)>>5);

  assert( thisOligo < states_.size() );

  // check if the position allows for the seed extension
  MatchPosition corrected_pos(0);
  const char* ext_chromName;
  const char* ext_contigName;

  getMatchPos( thisPos,
               ext_chromName,
               ext_contigName,
               corrected_pos );

  // correct for the offset within the
  uint adapted_position = 0;
  if( reverseFlag == 0 )
  {
      if( corrected_pos <= (MatchPosition)(seedOffsets_[seedNo]) )
      {
          return false;
      }
      adapted_position = thisPos - (MatchPosition)(seedOffsets_[seedNo]);
  }
  else
  {
      adapted_position = thisPos + (MatchPosition)(seedOffsets_[seedNo]);

      // additional check if the position points to a valid chromosome file
      MatchPosition corrected_pos(0);
      const char* ext_chromName;
      const char* ext_contigName;
      getMatchPos( adapted_position,
                   ext_chromName,
                   ext_contigName,
                   corrected_pos );

      if( string(ext_chromName)=="" )
      {
          return false;
      }
  }

  int existing_idx = -1;
  for( int i=states_[thisOligo].size()-1;i>=0;i-- )
  {
      if( ( abs( (int)(states_[thisOligo][i].pos_ - adapted_position) ) <= seedNo*SEED_DEVIATION ) &&
              states_[thisOligo][i].reverse_ == reverseFlag )
      {
          // also check if we have the hits are on the same strand
          // we found a hit that is being extended
          existing_idx = i;
          break;
      }
  }

  if( -1 == existing_idx )
  {
      // update matchType_.r[numError] correspondingly
      matchType_[thisOligo].r[numErrors]++;
      // insert the actual hit
      states_[thisOligo].push_back( SeedMatch(adapted_position, numErrors, reverseFlag, seedNo) );
//      std::cerr << "new seed match added " << thisOligo << " pos " << adapted_position
//               << " idx " << (int)states_[thisOligo].size()-1
//               << " errors " << (int)numErrors
//               << " new errors " << (int)states_[thisOligo].back().errors_
//               << " reverse " << (int)reverseFlag
//               << " seed " << (int)seedNo
//               << std::endl;
  }
  else
  {
      // update the hit that we found
      SeedMatch & tmp_type(states_[thisOligo][existing_idx]);

      // if the new seed hit has a lower number of errors, then update
      // the matchType_ entry, otherwise do nothing (if the number of
      // errors is the same or greater)
      if( numErrors < tmp_type.errors_ )
      {
//          std::cerr << "seed match extended " << thisOligo << " pos "  << adapted_position
//                   << " ori pos " << (int)states_[thisOligo][existing_idx].pos_
//                   << " idx " << (int)existing_idx
//                   << " errors " << (int)numErrors
//                   << " old errors " << (int)tmp_type.errors_
//                   << " reverse " << (int)reverseFlag
//                   << " seed " << (int)seedNo
//                   << std::endl;
          // the current seeds has a lower number of errors, update matchType_
          matchType_[thisOligo].r[numErrors]++;
          matchType_[thisOligo].r[tmp_type.errors_]--;
          tmp_type.errors_ = numErrors;
          tmp_type.lastSeed_ = seedNo;
//          tmp_type.pos_ = adapted_position;

      }
      else
      {
//          std::cerr << "seed match not extended " << thisOligo << " pos "  << adapted_position
//                   << " ori pos " << (int)states_[thisOligo][existing_idx].pos_
//                   << " idx " << (int)existing_idx
//                   << " errors " << (int)numErrors
//                   << " old errors " << (int)tmp_type.errors_
//                   << " reverse " << (int)reverseFlag
//                   << " seed " << (int)seedNo
//                   << std::endl;
          // keeping track of how many hits we used, subtract the
          // number of unused hits from the total number of hits ->
          // determine a consistent way of determining how many
      }
      ++tmp_type.seeds_;
      update_cnt++;
  }

  return true;
}


void StateMachine::clear( const OligoNumber& oligo )
{
  vector<SeedMatch> tmp;
  states_[oligo].swap(tmp);
}


void StateMachine::clear( void )
{
  vector< vector<SeedMatch> > tmp;
  states_.swap( tmp );
}


list<MultiMatch> StateMachine::getHits( const OligoNumber& oligo, int maxItems, int seedNumberBase)
{
  // mask the upper five bits (encode R/F/seed/errors)
  uint thisOligo=oligo&(((uint)~0)>>5);

  list<MultiMatch> l;

  // sort the hits that we found by some criterion, output the first x hits
  sort( states_[thisOligo].begin(),states_[thisOligo].end() );

  // exit if there's nothing to return
  if( states_[thisOligo].empty() )
	return l;

  // detect the maximal number of seed hits
  const unsigned int maxNumSeeds( states_[thisOligo].back().seeds_);
//  std::cerr << "max num seeds: " << maxNumSeeds << " oligo " << thisOligo << std::endl;

  for( vector<SeedMatch>::const_reverse_iterator rit = states_[thisOligo].rbegin();
      rit != states_[thisOligo].rend() && rit->seeds_ >= maxNumSeeds && maxItems;
      ++rit, --maxItems )
  {
//      std::cerr << " pushing " << rit->pos_
//                << " errors " << (int)rit->errors_ << " lastSeed " << (int)rit->lastSeed_ + seedNumberBase
//                << " reverse " << (int)rit->reverse_ << std::endl;
      // add 1 to our seed number as multiseed numbers begin from 1 (0 is reserved for singleseed matches)
      l.push_back( MultiMatch(rit->pos_,rit->errors_, rit->lastSeed_ + seedNumberBase, rit->reverse_) );
  }

  return l;
}


// combine StateMachine::matchType_[oligo] with descriptor
// update r[0],r[1],r[2]
bool StateMachine::combineMatchDescriptor( const OligoNumber& oligo,const MatchDescriptor& descriptor )
{
  for(uint i=0;i<2;i++ )
    {
      // if we have more than 10 matches, then add it up to StateMachine::matchType_
      if( descriptor.r[i] > 10 )
	{
	  matchType_[oligo].r[i] = ((matchType_[oligo].r[i]+descriptor.r[i])<=255)?(matchType_[oligo].r[i]+descriptor.r[i]):255;
	}
    }

  // the two error cases are a special case: If we have an exact
  // match, then we don't output them, but they are still counted in
  // the MultiMatchTable::addMatch method, but they are not inserted
  // into the StateMachine::insertSeedHit and therefore the
  // StateMachine.matchType_ vector is not aware of them.
  // as a workaround for now, combine both values

  // distinction if we have an exact match and if we don't have one
  if( matchType_[oligo].r[0] > 0 )
    {
      // we don't have to check whether the number of 2-error matches
      // is greater than 10, because if we have an exact hit, we don't
      // insert any 2-error matches to the state machine anyways!
      matchType_[oligo].r[2] = ((matchType_[oligo].r[2]+descriptor.r[2])<=255)?(matchType_[oligo].r[2]+descriptor.r[2]):255;
    }
  else
    {
      if( descriptor.r[2] > 10 )
	{
	  matchType_[oligo].r[2] = ((matchType_[oligo].r[2]+descriptor.r[2])<=255)?(matchType_[oligo].r[2]+descriptor.r[2]):255;
	}
    }


  return true;
}


bool StateMachine::compare_matches_( const SeedMatch& /*a*/,const SeedMatch& /*b*/ )
{


  return true;
}


void StateMachine::outputAccessStatistics(void)
{

}

} //namespace eland_ms
} //namespace casava
// end of ELAND_outer.cpp
