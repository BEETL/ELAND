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
 ** \file eland_ms/MultiMatch.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 ** Contains interface functionality from ELAND - specifically, stuff to do
 ** with scoring and storage of alignments
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_MULTI_MATCH_H
#define CASAVA_ELAND_MS_MULTI_MATCH_H

namespace casava
{
namespace eland_ms
{


struct MultiMatch
{
  MultiMatch( MatchPosition pos, uchar errors, uchar lastSeed, bool reverse ) :
    pos_(pos), errors_(errors), lastSeed_(lastSeed), reverse_(reverse){}
  MatchPosition pos_;

  uchar errors_  :2; //number of errors. Valid numbers are 0,1,2
  uchar lastSeed_:3; //seed that matched. Last one in case the initial seed got extended.
                     //  Valid are: 0 - singleseed, 1,2,3,4 multiseeds
  uchar reverse_ :1; //1 indicates a reverse match

  bool operator==( const MultiMatch& rhs ) const
  {
    return(pos_==rhs.pos_);
  } // ~bool operator==( const MultiMatch& rhs )

  bool operator<( const MultiMatch& rhs ) const
  {
    return(pos_<rhs.pos_);
  } // ~bool operator==( const MultiMatch& rhs )

  bool operator!=( const MultiMatch& rhs ) const
  {
    return(pos_!=rhs.pos_);
  } // ~bool operator!=( const MultiMatch& rhs )



} __attribute__ ((packed)); // ~struct MultiMatch

BOOST_STATIC_ASSERT(sizeof(MultiMatch)==sizeof(MatchPosition) + sizeof(uchar)); //only specializations are allowed to instantiate
//static bool lessThanMultiMatch
//( const MultiMatch& lhs, const MultiMatch& rhs )
//{
//  return (lhs.pos_<rhs.pos_);
//} // ~lessThanMultiMatch


} // namespace eland_ms
} // namespace casava



#endif // CASAVA_ELAND_MS_MULTI_MATCH_H
