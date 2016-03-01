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
 ** \file eland_ms/MatchDescriptor.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 ** Contains interface functionality from ELAND - specifically, stuff to do
 ** with scoring and storage of alignments
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_MATCH_DESCRIPTOR_H
#define CASAVA_ELAND_MS_MATCH_DESCRIPTOR_H
namespace casava
{
namespace eland_ms
{
// MatchDescriptor contains the remaining information necessary to describe
// a match between an oligo and the genome
//
// errorType:
//  bit 7: 1 = reverse complement match (ignored if MatchPosition is a repeat)
//  bit 2-5: pass at which match found, needed to convert error positions
//   in the suffix to positions in the oligo.
//  bit 0-1: number of errors in best match found so far (if zero, need to
//  check MatchPosition is non zero to distinguish from the case where no
//  match has yet been found)
//
// registers r[0], r[1], r[2]:
//  The roles of these registers vary according to the current state of the
//  match. The table below gives the contents of the three registers for each
//  of the seven states, and how the state responds upon receipt of each of
//  the three types of match.
//
//  Key:
//  'E'=exact match, '1'=match with 1 error, '2'=match with 2 errors
//  e0, e1 = positions of errors in suffix
//  #E, #1, #2 = counts of E,1,2 matches
//  xx = don't care
//  Copy= replace current match details, ->XX transition to state XX
//
// State match r0 r1 r2 | E               | 1               | 2
// ---------------------|-----------------|-----------------|----------------
// NM    noM    0  0  0 | Copy, r0++ ->UE | Copy, r1++ ->U1 | Copy, r2++ ->U2
// UE    mPos   1 #1 #2 |       r0++ ->RE |       r1++      |       r2++
// U1    mPos  e0  1 #2 | Copy, r0++ ->UE |       r1++ ->R1 |       r2++
// U2    mPos  e0 e1  1 | Copy, r0++ ->UE | Copy, r1++ ->U1 |       r2++ ->R2
// RE    xx    #E #1 #2 |       r0++      |       r1++      |       r2++
// R1    xx    xx #1 #2 | Copy, r0++ ->UE |       r1++      |       r2++
// R2    xx    xx xx #2 | Copy, r0++ ->UE | Copy, r1++ ->U1 |       r2++

struct MatchDescriptor
{
  uchar errorType;
  uchar r[3];
}; // ~struct MatchDescriptor

} // namespace eland_ms
} // namespace casava



#endif // CASAVA_ELAND_MS_MATCH_DESCRIPTOR_H
