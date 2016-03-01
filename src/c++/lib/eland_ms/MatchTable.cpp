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
 ** \file eland_ms/MatchTable.cpp
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#include "alignment/ELAND_unsquash.h"

#include "alignment/aligner.h"
#include "eland_ms/MatchRequest.hh"

#include "eland_ms/MatchTable.hh"
#include "eland_ms/StateMachine.hh"

#include "eland_ms/ElandDefines.hh"

namespace casava
{
namespace eland_ms
{

// set the limit on the size of the buffer of SeqRequests that should
// be processed
#define MR_REQUEST_SIZE 262144
#define REQUEST_SIZE 262144

// NOISE_TRESHOLD sets the _lower_ limit for the noise ratio:
// noise_ratio = (ad_mismatches-cigar_mismatches)/cigar_mismatches
#define NOISE_THRESHOLD 3.1

// if the Hamming distance between two reads is greater than MIN_HAMMING_DISTANCE,
// then we try to align the reads
#define MIN_HAMMING_DISTANCE 5



// PULLING OUT FRAGMENTS
// method for pulling out the genomic regions of interest
bool unsquashRequests( uint request_cnt,
                       ofstream& out,
                       FragmentFinder& getFragments,
                       const bool& align,
                       StringIndex& files,
                       vector<MatchRequest>& matches,
                       vector<SeqRequest>& frag_requests,
                       vector<char*>& reads,
                       const int& readLength,
                       const int& fragmentLength
                       )
{

  // set up a Q30 quality string
  boost::scoped_array<char> q30_qual_string(new char[readLength+1]);
  for( int q=0;q<readLength;q++ ) {
    q30_qual_string[q] = (char)94;
  }
  q30_qual_string[readLength] = '\0';

  // setting up the Hamming distance object
  Hamming h_dist;


  // setting up the gap settings for ELAND
  ifstream align_score(".align.scores");
  // the actual score values
  ga::alignment::ScoreType match     = 2;
  ga::alignment::ScoreType mismatch  = -1;
  ga::alignment::ScoreType gapopen   = 15;
  ga::alignment::ScoreType gapextend = 3;

  ga::alignment::Aligner aligner(match,mismatch,gapopen,gapextend,(ALIGN_DP_BAND/2),0,0 );
  aligner.init( readLength,fragmentLength,0,0 );
  aligner.allowInserts( 1 );
  aligner.allowDeletions( 1 );

#ifdef READ_PARAMETERS
  if( !align_score ) {
    if( aligner.readAlignScoreFile(align_score,match,mismatch,gapopen,gapextend) == true )
      {
        cerr << "read .align.scores, setting alignment scores to match=" << match
             << ",mismatch=" << mismatch
             << ",gapopen=" << gapopen
             << ",gapextend=" << gapextend << endl;
      }
  } // ~if(!align.score)
#endif



  // print the fragments that we still hold in the bufffer
  vector<char*> tmp_frags;
  vector<char*> frags;
  tmp_frags.resize( request_cnt );
  for (vector<char*>::iterator i(tmp_frags.begin());i!=tmp_frags.end();i++)
    {
      *i=new char[fragmentLength+1];
    } // ~for
  frags.resize( request_cnt );
  for (vector<char*>::iterator i(frags.begin());i!=frags.end();i++)
    {
      *i=new char[fragmentLength+1];
    } // ~for

  vector<char*> frags_cigar;
  frags_cigar.resize( request_cnt );
  for (vector<char*>::iterator i(frags_cigar.begin());i!=frags_cigar.end();i++)
    {
      *i=new char[2*fragmentLength+1];
    } // ~for


  // think of a more elegant way to accomplish this
  vector<int> pos_correction_begin;
  vector<int> pos_correction_end;
  pos_correction_begin.resize( request_cnt );
  pos_correction_end.resize( request_cnt );


  getFragments(frag_requests, reads, frags, files);//,(ALIGN_DP_BAND/2));

  Timer aligner_timer;


  if( align == false )
    {
      for( uint align_cnt=0;align_cnt<request_cnt;align_cnt++ )
        {
          int no_mismatches = 0;
          string alignment_descriptor = aligner.convertToAlignmentDescriptor( reads[ frag_requests[align_cnt].readNum_ ],
                                                                              string(frags[ frag_requests[align_cnt].requestNum_ ]).substr( (ALIGN_DP_BAND/2),readLength),
                                                                              no_mismatches);

          // copy over to frags
          strncpy( tmp_frags[ frag_requests[align_cnt].requestNum_ ],alignment_descriptor.c_str(),alignment_descriptor.size() );
          tmp_frags[ frag_requests[align_cnt].requestNum_ ][alignment_descriptor.size()] = '\0';
        }// for(uint align_cnt...)
    }
  else
    {

      for( uint align_cnt=0;align_cnt<request_cnt;align_cnt++ )
        {
          string cut_fragment = string(frags[ frag_requests[align_cnt].requestNum_ ]).substr( (ALIGN_DP_BAND/2),readLength );
          int no_old_ad_mismatches = 0;
          int no_new_ad_mismatches = 0;


          string alignment_descriptor = aligner.convertToAlignmentDescriptor( reads[ frag_requests[align_cnt].readNum_ ],
                                                                              cut_fragment,
                                                                              no_old_ad_mismatches );

          if( no_old_ad_mismatches < MIN_HAMMING_DISTANCE )
            {
              // do nothing
              strncpy( frags_cigar[ frag_requests[align_cnt].requestNum_ ],alignment_descriptor.c_str(),alignment_descriptor.size() );
              frags_cigar[ frag_requests[align_cnt].requestNum_ ][alignment_descriptor.size()] = '\0';
            }
          else
            {
              // do the full monty alignment

              // perform the alignment
              aligner(q30_qual_string.get(),reads[ frag_requests[align_cnt].readNum_ ],frags[ frag_requests[align_cnt].requestNum_ ],readLength,fragmentLength,(frag_requests[align_cnt].strand_=='F') );



              // compress to CIGAR format
              int offset_begin_indel = 0;
              int offset_end_indel = 0;

              string new_alignment_descriptor = aligner.convertToNewAlignmentDescriptor( aligner.xt_,
                                                                                         aligner.yt_,
                                                                                         no_new_ad_mismatches,
                                                                                         offset_begin_indel,
                                                                                         offset_end_indel );

              // the old alignment_descriptor has been calculated above!


              //              cerr << endl << endl << "no_of_mismatches = " << no_new_ad_mismatches << endl;
              //              cerr << endl <<  "offset begin/offset end = " << offset_begin_indel << "\t" << offset_end_indel << endl;


              bool old_descriptor = false;

              // count the number of mismatches
              int cigar_mismatches = no_new_ad_mismatches;
              int md_mismatches = no_old_ad_mismatches;
              double noise_ratio = 0.0;


              // now decide whether we want to keep the new one or not
              if( new_alignment_descriptor.size() == 0 || (offset_begin_indel>30) || (offset_end_indel>30) )
                {
                  old_descriptor = true;
                } else
                  {
                    noise_ratio = (double)((double)(md_mismatches-cigar_mismatches)/(double)(cigar_mismatches));

                    //cout << noise_ratio << endl;

                    if( noise_ratio < NOISE_THRESHOLD ) {
                      old_descriptor = true;
                    }
                  }

              // if we are fine with the old descriptor, then copy it over
              if( old_descriptor == true ) {
                new_alignment_descriptor = alignment_descriptor;
                pos_correction_begin[ frag_requests[align_cnt].requestNum_ ] = 0;
                pos_correction_end[ frag_requests[align_cnt].requestNum_ ] = 0;
              } else
                {
                  // if we chose the new alignment descriptor, we perhaps
                  // have to adjust the match position
                  pos_correction_begin[ frag_requests[align_cnt].requestNum_ ] = ((ALIGN_DP_BAND/2)-offset_begin_indel);
                  pos_correction_end[ frag_requests[align_cnt].requestNum_ ] = ((ALIGN_DP_BAND/2)-offset_end_indel);

                }



              // copy over to frags
              strncpy( frags_cigar[ frag_requests[align_cnt].requestNum_ ],new_alignment_descriptor.c_str(),new_alignment_descriptor.size() );
              frags_cigar[ frag_requests[align_cnt].requestNum_ ][new_alignment_descriptor.size()] = '\0';

            } // ~else( full_monty_alignment)

        } //~for loop
    } // ~align

  cerr << "alignment done for the moment " << aligner_timer << endl;

  // print the match requests together with the fragment
  int frag_idx = 0;

  if( align == false ) {
    for( uint i=0;i<matches.size();i++ )
      {
        matches[i].print( out,tmp_frags,frag_idx,pos_correction_begin,pos_correction_end );
      }
  } else {
    for( uint i=0;i<matches.size();i++ )
      {
        matches[i].print( out,frags_cigar,frag_idx,pos_correction_begin,pos_correction_end );
      }
  }

  // all fragments should have been used
  if( frag_idx != (int)frags.size() )
    {
      cerr << "frag_idx            = " << frag_idx << endl
           << "frags.size          = " << frags.size() << endl
           << "not all fragments used, this should not happen." << endl;
      exit(1);
    }

  // cleaning up
  for (vector<char*>::iterator i(reads.begin());i!=reads.end();i++)
    {
      delete [] *i;
    } // ~for
  for (vector<char*>::iterator i(frags.begin());i!=frags.end();i++)
    {
      delete [] *i;
    } // ~for
  for (vector<char*>::iterator i(tmp_frags.begin());i!=tmp_frags.end();i++)
    {
      delete [] *i;
    } // ~for
  for (vector<char*>::iterator i(frags_cigar.begin());i!=frags_cigar.end();i++)
    {
      delete [] *i;
    } // ~for

  {
    vector<char*> c_tmp;
    vector<int> i_tmp;
    vector<int> ii_tmp;
    reads.swap(c_tmp);
    frags.swap(c_tmp);
    tmp_frags.swap(c_tmp);
    frags_cigar.swap(c_tmp);
    pos_correction_begin.swap(i_tmp);
    pos_correction_end.swap(ii_tmp);

    reads.clear();
    frags.clear();
    tmp_frags.clear();
    frags_cigar.clear();
    pos_correction_begin.clear();
    pos_correction_end.clear();
  }

  return true;
}




// isInterested__: hash table uses this function to ask
// if it should hash an oligo
bool MatchTable::isInterested__( int oligoNum, int PASS, bool hasNs )
{
  if (PASS==0)
  {
    // if no repeat or QC info set, interested in oligo
    if (size()==0) return true;
    if ((hasQualityFailed__(oligoNum))
        ||(isRepeatMasked__(oligoNum)))
      return false;
    else
      return true;
  } // ~if
  else
  {
    if ((hasQualityFailed__(oligoNum))
        ||(isRepeatMasked__(oligoNum)))
      return false;
    if (!hasNs)
    {
      //      if (this->matchPosition_[oligoNum]==noMatch)
      //   return true;
      if (  ((matchType_[oligoNum].errorType&0x3)==0)
            && (matchType_[oligoNum].r[0]>1))
      {
        //        cout << ((int)matchType_[oligoNum].r[0]) << " " << ((int)matchType_[oligoNum].r[1]) << " " << ((int)matchType_[oligoNum].r[2]) << endl;

        return false;
      } // ~if
      else
      {
        return true;
      } // ~else
    } // ~if
    else
    {
      return true;
    } // ~else

  } // ~else
} // ~MatchTable::wantMatch__( int PASS, bool hasNs )


#if 0
// addMatch: process a match
// oligoNum - number of the oligo found
// matchPos - where it was found in the genome
// matchType - number and position of any errors in the match
// passNumber - which pass of the genome the match was found on
void MatchTable::addMatch( const TableEntryData& entry,
			   const Word thisMask,
			   const MatchPosition matchPos,
			   const FragmentErrorType errorLow,
			   const FragmentErrorType errorHigh )
			   // const int numErrorsz
{
  //  cout << "addMatch: " << " " <<  matchPos-16777216 << " "
  //     << errorLow << " " << errorHigh << endl;


  const uint8_t numErrors
    ( (errorLow>oneError__)+(errorLow>noErrors__)
      +(errorHigh>oneError__)+(errorHigh>noErrors__));

  //   cout << "am " << errorLow << " " << errorHigh << " " << noErrors__ << " " << oneError__ << endl;
  //   cout << "am " << matchPos-16777216 << " " << numErrors << " " << (errorLow&errorPosMask1__) << " "
  // << ((errorLow&errorPosMask2__)>>8) << " "
  // << (errorHigh&errorPosMask1__) << " "
  //  << ((errorHigh&errorPosMask2__)>>8) << endl;



  //  assert(numErrors==numErrorsz);
  assert(numErrors<=2);

  // For oligos with no Ns:
  // On partition 0, interested in 0, 1 and 2 error matches
  // On partition 1, interested in 0 or 1 error matches
  // On partitions 2 to 5, interested in 2 error matches only

  const OligoNumber oligoNum(entry.position
			     &((~isReverseOligo)));
  const bool reverseFlag((entry.position&isReverseOligo)!=0);
  //    oligoNum&=(~isReverseOligo);
    //    cout << 'Y' << endl;

    //    if (this->matchPosition_[oligoNum] != matchPos )
    //    {
#if (NUM_THREADS>1)
  pthread_mutex_lock(&mutex_);
#endif
      if (    ( numErrors < (matchType_[oligoNum].errorType&0x3))
	    || (this->matchPosition_[oligoNum]==noMatch))
      {
	// Then match is better than any match found so far
	if ( numErrors==2 )
	{
	  // store positions of both errors

          //	  matchType_[oligoNum].r[0]
          //   = errorLow&errorMask1;
          //  matchType_[oligoNum].r[1]
          //   = errorHigh&errorMask1;
          // TC 09.08.05 - above lines have a bug - they don't store both
          // match positions if both matches occur in the same fragment

	  matchType_[oligoNum].r[0]
	    = errorLow&errorInfoMask1__;
          if (matchType_[oligoNum].r[0]==0)
          { // then both errors must be in errorHigh
            matchType_[oligoNum].r[0]
              = errorHigh&errorInfoMask1__;
            matchType_[oligoNum].r[1]
              = (errorHigh&errorInfoMask2__)>>errorBits__;
          } // ~if
          else
          { //  assume 2nd error is in upper bits of errorHigh...

            matchType_[oligoNum].r[1]
              = (errorLow&errorInfoMask2__)>>errorBits__;

            // ... but if it's zero ...
            if (matchType_[oligoNum].r[1]==0)
            { // ... 2nd error must instead be in lower bits of errorLow
              matchType_[oligoNum].r[1]
                = (errorHigh&errorInfoMask1__);
            } // ~if
          } // ~else

          assert(matchType_[oligoNum].r[0] !=0);
          assert(matchType_[oligoNum].r[1] !=0);
          assert((matchType_[oligoNum].r[0]&ucharErrorTypeMask) !=0);
          assert((matchType_[oligoNum].r[1]&ucharErrorTypeMask) !=0);


	} // ~if
	else if( numErrors==1 )
	{
	  // store position of single error
	  matchType_[oligoNum].r[0]
	    = (errorLow|errorHigh)&errorInfoMask1__;
          //	    = (errorLow|errorHigh)&errorMask1;
	}
	else if( numErrors==0 )
	{
	  // if only 2 error matches found before, set r[1] to zero
	  // to signify no 1 error matches found so far
	  // else leave r[1] alone! - otherwise we lose all 1 error
	  // matches counted thus far (bug fix over previous version)
	  // TC 05.03.04
	  if ((matchType_[oligoNum].errorType&0x3)==2)
	    matchType_[oligoNum].r[1]=0;
	}
	//	else assert(1==0);

	// this is the first occurrence of this match type
	matchType_[oligoNum].r[numErrors]=1;

	// store its position
	this->matchPosition_[oligoNum] = matchPos;

	Word nInfo(0);
	if (thisMask!=0)
	{
	  nInfo |= 0x1
	    *(((thisMask&entry.suffix.ui)&(thisMask&nonIncorpFirstN))!=0);
	  nInfo |= 0x2
	    *(((thisMask&entry.suffix.ui)&(thisMask&nonIncorpSecondN))!=0);
          nInfo ^= 0x3*reverseFlag;
	} // ~if

	matchType_[oligoNum].errorType
	  = (uchar)(numErrors|(nInfo<<2)|(0x80*reverseFlag));


      } // ~if
      else
      {
	// increment count if we are not already at capacity
	// if match is within 1 of an existing match, ignore it


	//	matchType_[oligoNum].r[numErrors]
	//	  +=( (matchType_[oligoNum].r[numErrors]!=0xFF)
	//	      &&( ((this->matchPosition_[oligoNum]-matchPos)
	//	   *(this->matchPosition_[oligoNum]-matchPos))
	//		  >1));
	// TC 31.10.03
	// Above version can give an erroneous false on very rare occasions
	// (and so miss matches) if mP_[oN]-mP)^2 = -1,0,+1 mod 2^32
	// so instead use mp_[oN]-mP = -1,0,+1 iff mP_[oN]-mP+1 = 0,1,2
	// ... also avoids problem with unsigned quantities
	matchType_[oligoNum].r[numErrors]
	  +=( (matchType_[oligoNum].r[numErrors]!=0xFF)
	      && ( (this->matchPosition_[oligoNum]-matchPos+1) > 2 ));
      } // ~else

  //  cout << this->matchPosition_[oligoNum]-16777216 << " " << (uint) matchType_[oligoNum].errorType << " " << (uint) matchType_[oligoNum].r[0] << " "  << (uint) matchType_[oligoNum].r[1] << " " << (uint) matchType_[oligoNum].r[2] << endl;
#if (NUM_THREADS>1)
      pthread_mutex_unlock(&mutex_);
#endif
} // ~void MatchTable::addMatch( uint oligoNum, MatchPosition oligoPos )
#endif


// adjustMatchPos: make small corrections to the match position to correct
// for various quirks of the matching process (mainly to do with
// incorporation errors and reverse complements)
void MatchTable::adjustMatchPos( const char* pOligo,
		       const char dirChar,
		       char& firstN,
		       char& secondN,
		       MatchPosition& matchPos)
{
  uint headSize, tailSize, numInternalNs;
  countNs( pOligo, OLIGO_LEN_, headSize, tailSize, numInternalNs );

  // Print out '.' as the second N info character to signify there is only
  // one internal N
  if (numInternalNs<=1)
  {
    //    assert(secondN=='D');
    secondN='.';
  } // ~if

  if (numInternalNs==0)
  {
    assert(firstN=='D');
    firstN='.';
  } // ~if

  if (dirChar=='R')
  {
    // Need to shift match position forward by 1 for each non incorporation
    // (compensate for the extra Ns)
    //   cout << matchPos << " -> ";
    matchPos+=(firstN=='I');
    matchPos+=(secondN=='I');
    //  cout << matchPos << endl;
    // All Ns shoved at beginning which pushes match pos out - correct
    matchPos+=(headSize+tailSize);

    //    matchPos-=tailSize; // TC 27.1.6
  } // ~if
  else
  {
    assert(dirChar=='F');
    //    matchPos-=headSize; // TC 27.1.6
  } // ~else
} // ~MatchTable::adjustMatchPos


// For each sequence check if there is at least one match position found
bool MatchTable::getUnmappedReads( vector<bool>& /*unmapped*/ )
{

  return true;
}



#if 0
void MatchTableVerbose::print( OligoSource& oligos,
			       MatchPositionTranslator& getMatchPos,
			       const vector<string>& ,
			       const vector<MatchPosition>& blockStarts,
			       const SuffixScoreTable& ,
			       int )
{
  //  ErrorPositionMapper mapErrors(scoreTable, oligoLength);

  MatchPosition subtractTable[256];
  int chromTable[256];
  int thisChrom(0);
  uint j;
  uint errorPos1, errorPos2;

  //  bool isUniqueExact, isUniqueError;
  //  char buf[lineLength];
  const char* pOligo;
  ushort numErrors;
  char dirChar, firstN = 0, secondN = 0;
  Word nInfo;
  MatchPosition thisBlock, matchPos;

  // extract match info
  MatchPosition extractedMatchPos;
  const char* extractedChromName;
  const char* extractedContigName;


  for (vector<MatchPosition>::const_iterator i(blockStarts.begin());
       i!=blockStarts.end()-1;++i)
  {
    for ( j = *i; j != *(i+1); j+=blockSize)
    {
      subtractTable[j>>blockShift]=*i;
      chromTable[j>>blockShift]=thisChrom;
    } // ~for j
    ++thisChrom;
  } // ~for i

  oligos.rewind();

  //  for (int i(1); i <= numOligos_ ; i++ )
  for (uint i(1); i < this->matchPosition_.size() ; i++ )
  {

      if ( (pOligo=oligos.getNextOligoSelect(true,false)) == NULL )
    {
      cerr << "OligoInfo: unexpectedly ran out of names!" << endl;
      exit (1);
    } // ~if

    fprintf( this->pOut_, "%d\t%s\t%s\t", i, oligos.getLastName(), pOligo );

    if (this->matchPosition_[i]>=blockRepeat)
    {
      if (this->matchPosition_[i]==qualityFailed)
	fprintf( this->pOut_, "QC");
      else if (this->matchPosition_[i]==repeatMasked)
	fprintf( this->pOut_,"RM");
      else
	fprintf( this->pOut_,"RB\t%u", this->matchPosition_[i]-blockRepeat);
    } // ~if
    else if (this->matchPosition_[i]==noMatch)
      fprintf( this->pOut_,"NM\t0\t0\t0");
    else
    {
      numErrors=(this->matchType_[i].errorType&0x3);
      fprintf( this->pOut_,"%c%c\t%u\t%u\t%u",
	     ((this->matchType_[i].r[numErrors]>1)?'R':'U'),
	     '0'+numErrors,
	     ((numErrors==0)?((uint)this->matchType_[i].r[0]):0),
	     ((numErrors<=1)?((uint)this->matchType_[i].r[1]):0),
	     (uint)this->matchType_[i].r[2] );

      assert(this->matchType_[i].r[numErrors]!=0);
      if (this->matchType_[i].r[numErrors]==1)
      {
	// then we have a unique match - print its position

	// get direction
	dirChar = (((this->matchType_[i].errorType&0x80)!=0)?'R':'F');
	// get info on Ns
	nInfo=((Word)(this->matchType_[i].errorType&0x7C))>>2; // NNYYYYYN
	switch (nInfo)
	{
	case (noNonIncorps&0x3):
	  firstN='D';
	  secondN='D';
	  //	  fprintf( this->pOut_,"\tDD");
	  break;
	case (nonIncorpFirstN&0x3):
	  firstN='I';
	  secondN='D';
	  //	  fprintf( this->pOut_,"\tID");
	  break;
	case (nonIncorpSecondN&0x3):
	  firstN='D';
	  secondN='I';
	  //	  fprintf( this->pOut_,"\tDI");
	  break;
	case (nonIncorpBothNs&0x3):
	  firstN='I';
	  secondN='I';
	  //	  fprintf( this->pOut_,"\tII");
	  break;
	default:
	  assert(1==0);
	} // ~switch

	thisBlock = this->matchPosition_[i] >> blockShift;
	matchPos=this->matchPosition_[i]-subtractTable[thisBlock];


        getMatchPos( this->matchPosition_[i],
		     extractedChromName,
		     extractedContigName,
		     extractedMatchPos);

	//	if(extractedMatchPos!=matchPos)
	//	{
	//	  cerr << extractedMatchPos << " " << matchPos << endl;
	//	}


        this->adjustMatchPos( pOligo, dirChar, firstN, secondN, extractedMatchPos );

	fprintf
          ( this->pOut_,"\t%s%s\t%u\t%c\t%c%c",
            extractedChromName,
            extractedContigName,
            extractedMatchPos,
            dirChar,
            firstN,
            secondN );

#ifdef XXX
	this->adjustMatchPos( pOligo, dirChar, firstN, secondN, matchPos );

	fprintf
          ( this->pOut_,"\t%s\t%u\t%c\t%c%c",
            chromNames[chromTable[thisBlock]].c_str(),
            matchPos,//this->matchPosition_[i]-subtractTable[thisBlock],
            dirChar,
            firstN,
            secondN );
#endif

	if ((this->matchType_[i].errorType&0x3)==1)
	{
          errorPos1=(uint)(this->matchType_[i].r[0]&ucharErrorPosMask);
	  fprintf( this->pOut_,"\t%u", errorPos1);
          // fprintf( this->pOut_,"\t%u", (uint)matchType_[i].r[0]);
	} // ~if
	else if ((this->matchType_[i].errorType&0x3)==2)
	{
          errorPos1=(uint)(this->matchType_[i].r[0]&ucharErrorPosMask);
          errorPos2=(uint)(this->matchType_[i].r[1]&ucharErrorPosMask);
	  fprintf( this->pOut_,"\t%u\t%u",
                   min(errorPos1, errorPos2),
                   max(errorPos1, errorPos2) );
          //	  fprintf( this->pOut_,"\t%u\t%u",
          // (uint)min(matchType_[i].r[0],matchType_[i].r[1]),
          // (uint)max(matchType_[i].r[0],matchType_[i].r[1]) );
	} // ~else if

      } // ~if
    } // ~else
    if (0 > fprintf( this->pOut_,"\n"))
    {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "failed to write ELAND output"));
    }
  } // ~for i
} // ~void MatchTableVerbose::print


// parseErrorInfo: parses the information stored in the error char
// - bottom 6 bits are the error position
// - using the top 2 bits and the original read it works
//   out what the erroneous base should have been
void MatchTableAssembly::parseErrorInfo
( const char* oligoBuf,
  const char errorInfo,
  const char dirChar,
  uint& errorPos,
  uchar& shouldBe )
{
  errorPos=(uint)(errorInfo&ucharErrorPosMask);
  //  cout << endl << oligoBuf << endl;
  //  cout << errorPos << endl;
  // cout << (uint)whichBase[oligoBuf[errorPos-1]] << endl;;
  //  cout << (uint)((errorInfo&ucharErrorTypeMask)>>errorPosBits__) << endl;


  //  uint shouldBeIndex
  //  ( whichBase[oligoBuf[errorPos-1]]
  //    ^((uint)((errorInfo&ucharErrorTypeMask)>>errorPosBits__)) );

  uint shouldBeIndex(((dirChar=='F')?(errorPos-1):(OLIGO_LEN_-errorPos)));
  shouldBeIndex=whichBase[(uint)oligoBuf[shouldBeIndex]];
  shouldBeIndex^=((uint)((errorInfo&ucharErrorTypeMask)>>errorPosBits__));
  shouldBeIndex^=( ((uint)0x3)*(dirChar=='R') );

  //  cout << shouldBeIndex << endl;
  //  assert(shouldBeIndex<numDifferentBases);
  //  TBD Next line is a temporary fix. Should never occur - find out why
  if (shouldBeIndex>=numDifferentBases)
    shouldBeIndex=0;

  shouldBe=baseNames[shouldBeIndex];
} // ~MatchTableAssembly::parseErrorInfo



void MatchTableAssembly::print( OligoSource& oligos,
				MatchPositionTranslator& getMatchPos,
				const vector<string>& ,
				const vector<MatchPosition>& blockStarts,
				const SuffixScoreTable& ,
				int )
{

  //  MatchPositionTranslator getMatchPos( chromNames, blockStarts );


  MatchPosition subtractTable[256];
  int chromTable[256];
  int thisChrom(0);
  uint j;
  uint errorPos1, errorPos2;
  uchar shouldBe1, shouldBe2;

  //  bool isUniqueExact, isUniqueError;
  //  char buf[lineLength];
  const char* pOligo;
  ushort numErrors;
  char dirChar, firstN, secondN;
  Word nInfo;
  MatchPosition thisBlock, matchPos;

  // extract match info
  MatchPosition extractedMatchPos;
  const char* extractedChromName;
  const char* extractedContigName;


  for (vector<MatchPosition>::const_iterator i(blockStarts.begin());
       i!=blockStarts.end()-1;++i)
  {
    for ( j = *i; j != *(i+1); j+=blockSize)
    {
      subtractTable[j>>blockShift]=*i;
      chromTable[j>>blockShift]=thisChrom;
    } // ~for j
    ++thisChrom;
  } // ~for i

  oligos.rewind();

  //  for (int i(1); i <= numOligos_ ; i++ )
  for (uint i(1); i < this->matchPosition_.size() ; i++ )
  {

    if ( (pOligo=oligos.getNextOligoSelect(true,false)) == NULL )
    {
      cerr << "OligoInfo: unexpectedly ran out of names!" << endl;
      exit (1);
    } // ~if

    fprintf( this->pOut_, "%s\t%s\t", oligos.getLastName(), pOligo );

    if (this->matchPosition_[i]>=blockRepeat)
    {
      if (this->matchPosition_[i]==qualityFailed)
	fprintf( this->pOut_, "QC");
      else if (this->matchPosition_[i]==repeatMasked)
	fprintf( this->pOut_,"RM");
      else
	fprintf( this->pOut_,"RB\t%u", this->matchPosition_[i]-blockRepeat);
    } // ~if
    else if (this->matchPosition_[i]==noMatch)
      fprintf( this->pOut_,"NM\t0\t0\t0");
    else
    {
      numErrors=(this->matchType_[i].errorType&0x3);
      fprintf( this->pOut_,"%c%c\t%u\t%u\t%u",
	     ((this->matchType_[i].r[numErrors]>1)?'R':'U'),
	     '0'+numErrors,
	     ((numErrors==0)?((uint)this->matchType_[i].r[0]):0),
	     ((numErrors<=1)?((uint)this->matchType_[i].r[1]):0),
	     (uint)this->matchType_[i].r[2] );

     assert(this->matchType_[i].r[numErrors]!=0);
      if (this->matchType_[i].r[numErrors]==1)
      {
	// then we have a unique match - print its position

	// get direction
	dirChar = (((this->matchType_[i].errorType&0x80)!=0)?'R':'F');
	// get info on Ns
	nInfo=((Word)(this->matchType_[i].errorType&0x7C))>>2; // NNYYYYYN
	switch (nInfo)
	{
	case (noNonIncorps&0x3):
	  firstN='D';
	  secondN='D';
	  //	  fprintf( this->pOut_,"\tDD");
	  break;
	case (nonIncorpFirstN&0x3):
	  firstN='I';
	  secondN='D';
	  //	  fprintf( this->pOut_,"\tID");
	  break;
	case (nonIncorpSecondN&0x3):
	  firstN='D';
	  secondN='I';
	  //	  fprintf( this->pOut_,"\tDI");
	  break;
	case (nonIncorpBothNs&0x3):
	  firstN='I';
	  secondN='I';
	  //	  fprintf( this->pOut_,"\tII");
	  break;
	default:
	  assert(1==0);
	} // ~switch

	thisBlock = this->matchPosition_[i] >> blockShift;
	matchPos=this->matchPosition_[i]-subtractTable[thisBlock];

        getMatchPos( this->matchPosition_[i],
		     extractedChromName,
		     extractedContigName,
		     extractedMatchPos);

	//	if(extractedMatchPos!=matchPos)
	//	{
	//  cerr << extractedMatchPos << " " << matchPos << endl;
	//	}

	// TC 4.4.8 - this adjustment now done in Unsquash.cpp
	//	adjustMatchPos( pOligo, dirChar, firstN, secondN,
	//		extractedMatchPos );

	fprintf( this->pOut_,"\t%s%s\t%u\t%c\t%c%c",
		 extractedChromName,
		 extractedContigName,
		 extractedMatchPos,
		 dirChar,
		 firstN,
		 secondN );

#ifdef XXX
	adjustMatchPos( pOligo, dirChar, firstN, secondN, matchPos );

	fprintf( this->pOut_,"\t%s\t%u\t%c\t%c%c",
		 chromNames[chromTable[thisBlock]].c_str(),
		 matchPos,//this->matchPosition_[i]-subtractTable[thisBlock],
		 dirChar,
		 firstN,
		 secondN
	       );
#endif
	if ((this->matchType_[i].errorType&0x3)==1)
	{
          parseErrorInfo
            ( pOligo, this->matchType_[i].r[0], dirChar, errorPos1, shouldBe1 );
          fprintf( this->pOut_, "\t%u%c", errorPos1, shouldBe1 );
          //	  fprintf( this->pOut_,"\t%u", (uint)this->matchType_[i].r[0]);
	} // ~if
	else if ((this->matchType_[i].errorType&0x3)==2)
	{
          parseErrorInfo
            ( pOligo, this->matchType_[i].r[0], dirChar, errorPos1, shouldBe1 );
          parseErrorInfo
            ( pOligo, this->matchType_[i].r[1], dirChar, errorPos2, shouldBe2 );
          //	  fprintf( this->pOut_,"\t%u\t%u",
          //	 (uint)min(this->matchType_[i].r[0],this->matchType_[i].r[1]),
          //	 (uint)max(this->matchType_[i].r[0],this->matchType_[i].r[1]) );
          if (errorPos1<errorPos2)
          {
            fprintf( this->pOut_,"\t%u%c\t%u%c",
                     errorPos1, shouldBe1, errorPos2, shouldBe2);
          } // ~if
          else
          {
            fprintf( this->pOut_,"\t%u%c\t%u%c",
                     errorPos2, shouldBe2, errorPos1, shouldBe1);
          } // ~else

	} // ~else if

      } // ~if
    } // ~else
    if (0 > fprintf( this->pOut_,"\n"))
    {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "failed to write ELAND output"));
    }
  } // ~for i

} // ~void MatchTableAssembly::print
#endif


void MatchTableMulti::initializeTmpFiles
( const char* tmpFilePrefix )
{
  std::string tmpFilePrefixString, tmpFileName;
  if (tmpFilePrefix==NULL)
  {
    if ((pOligoNum_=casava_tmpfile())==NULL)
    {
      BOOST_THROW_EXCEPTION(cc::IoException(errno, "MatchTableMulti could not open oligo num temp file."));
    }
    if ((pMatchType_=casava_tmpfile())==NULL)
    {
      BOOST_THROW_EXCEPTION(cc::IoException(errno, "MatchTableMulti could not open match type temp file."));
    }
  }
  else
  {
    tmpFilePrefixString=tmpFilePrefix;
    tmpFileNameOligoNum_=tmpFilePrefixString + ".num";
    if ((pOligoNum_=fopen(tmpFileNameOligoNum_.c_str(), "w+b"))==NULL)
    {
      BOOST_THROW_EXCEPTION(cc::IoException(errno, "MatchTableMulti could not open oligo num temp file" + tmpFileNameOligoNum_));
    }

    tmpFileNameMatchType_=tmpFilePrefixString + ".type";
    if ((pMatchType_=fopen(tmpFileNameMatchType_.c_str(), "w+b"))==NULL)
    {
      BOOST_THROW_EXCEPTION(cc::IoException(errno, "MatchTableMulti could not open match type temp file" + tmpFileNameMatchType_));
    }

  }

  std::cerr << "Built MatchTableMulti: will store at most "
       <<  maxNumMatchesExact_ << ","
       <<  maxNumMatchesOneError_ << ","
       <<  maxNumMatchesTwoErrors_ << " 0,1,2 error matches per read"
       << endl;
}

MatchTableMulti::~MatchTableMulti()
{
  fclose (pOligoNum_);
  fclose (pMatchType_);
} // MatchTableMulti::~MatchTableMulti



// For each sequence check if there is at least one match position found
bool MatchTableMulti::getUnmappedReads( vector<bool>& unmapped )
{
  unmapped.clear();
  hyperhyper_.clear();
  unmapped.resize(this->matchPosition_.size(),true);
  hyperhyper_.resize(this->matchPosition_.size(),false);

  fseek( pOligoNum_, 0, SEEK_SET);

  uint thisCode;
  uint thisOligo;

  while (1)
  {
    if (1 != fread( &thisCode, sizeof(uint), 1, pOligoNum_))
    {
        const int currentError = errno;
        if (feof(pOligoNum_)!=0)
        {
            break;
        }
        else
        {
            BOOST_THROW_EXCEPTION(cc::IoException(currentError, "failed to read Oligo Num"));
        }
    }

    thisOligo=thisCode&(((uint)~0)>>5);
    //    cout << thisCode << " " << thisOligo << " " << (int)reverseFlag << " "<< thisPos << " " << numErrors << endl;
    assert(thisOligo<this->matchPosition_.size());
    unmapped[thisOligo] = false;

  } // ~while


  // we will map QC reads in the
  for( uint i=1;i<this->matchPosition_.size();i++ )
  {
      if( (this->matchPosition_[i]==qualityFailed) || (this->matchPosition_[i]==repeatMasked) )
      {
          unmapped[i]=false;
      }



      // if sensitive_ == false, we don't take 255:255:255 into account!
      if( (
              (this->matchType_[i].r[0]>maxNumMatchesExact_) &&
           (this->matchType_[i].r[1]>maxNumMatchesOneError_) &&
           (this->matchType_[i].r[2]>maxNumMatchesTwoErrors_)
          )
          &&
          ( this->sensitive_ ||
              ((this->matchType_[i].r[0]<255) &&
              (this->matchType_[i].r[1]<255) &&
              (this->matchType_[i].r[2]<255)
              ))
          )
      {
          unmapped[i]=true;
          this->hyperhyper_[i]=true;
      }

  }


  // count how many occurrences there are
  int count_reads_to_map = 1;
  for( uint i=1;i<unmapped.size();i++ )
  {
      if( unmapped[i]==true ) { count_reads_to_map++; }
  }

  // resize the translator
  this->translator_.resize( count_reads_to_map, 0 );
  //translator_.resize(this->matchPosition_.size(),-1);

  // fill up the translator
  int unmapped_cnt = 1;

  for( uint i=1;i<this->matchPosition_.size();i++ )
    {
        if( unmapped[i]==true )
        {
            //cerr << "mapping " << i << " -> " << unmapped_cnt << endl;
            this->translator_[unmapped_cnt++] = i;
        }
    }

  // calculate the ratio of reads that make it to the second tier
  double secondTierRatio = (double)count_reads_to_map/(double)unmapped.size();
  cerr << "ratio of reads going into the second tier = " << secondTierRatio;


  // store unmapped into unmapped_ which is a member MatchTableMulti
  this->unmapped_ = unmapped;

  return true;
}



// MODIFIED MULTISEED CODE
//#ifdef MODIFIED_MULTISEED_CODE
void MatchTableMulti::
addMatch(MatchIter i, const MatchIter i_end)
{
    static const OligoNumber onMask((~isReverseOligo)&(~seed_bits[3]));

    for(;i!=i_end;++i) {
        
        assert(i->numErrors<=2);

  // For oligos with no Ns:
  // On partition 0, interested in 0, 1 and 2 error matches
  // On partition 1, interested in 0 or 1 error matches
  // On partitions 2 to 5, interested in 2 error matches only


  // splitPrefixMaskLow is not necessary anymore
  // 22/04/09 markus removing splitPrefixMaskLow
  //  const OligoNumber oligoNum(entry.position
  //  			     &((~isReverseOligo)&splitPrefixMaskLow));

        const OligoNumber oligoNum=(i->position&onMask);

        this->matchType_[oligoNum].r[i->numErrors] += (this->matchType_[oligoNum].r[i->numErrors]!=0xFF);

        if (   ((i->numErrors==0)
                &&(this->matchType_[oligoNum].r[0]<=maxNumMatchesExact_))
               || ((i->numErrors==1)
                   &&(this->matchType_[oligoNum].r[0]<=maxNumMatchesExact_)
                   &&(this->matchType_[oligoNum].r[1]<=maxNumMatchesOneError_))
               || ((i->numErrors==2)
                   &&(this->matchType_[oligoNum].r[0]==0)
                   &&(this->matchType_[oligoNum].r[1]<=maxNumMatchesOneError_)
                   &&(this->matchType_[oligoNum].r[2]<=maxNumMatchesTwoErrors_)))
            {
                // TBD: cache then bulk write
                const uint32_t matchCode=((static_cast<uint32_t>(i->numErrors)<<30)|
                           (((i->position&isReverseOligo)!=0)<<29)|
                           (((i->position&(~isReverseOligo))>>29)<<27)|oligoNum);
                ++matchesStored_;

                if( 1 != fwrite (&matchCode,sizeof(uint),1,pOligoNum_) )
                    {
                        cerr << "ERROR: couldn't write match code to temp file: " << strerror(errno) << endl;
                        exit(2);
                    }
                if( 1 != fwrite (&(i->matchPos),sizeof(MatchPosition),1,pMatchType_) )
                    {
                        cerr << "ERROR: couldn't write match pos to temp file: " << strerror(errno) << endl;
                        exit(2);
                    }
            }
    }


#if (NUM_THREADS>1)
      pthread_mutex_unlock(&mutex_);
#endif
} // ~void MatchTableMulti::addMatch( uint oligoNum, MatchPosition oligoPos )
//#endif // MODIFIED_MULTISEED_CODE


void MatchTableMulti::print( OligoSource& oligos,
				MatchPositionTranslator& getMatchPos,
				const vector<string>& ,
				const vector<MatchPosition>& blockStarts,
				const SuffixScoreTable& ,
				int )
{
  cerr << "Info: " << matchesStored_
       << " matches were stored" << endl;
  cerr << "Info: " << ftell(pOligoNum_)
       << " bytes of temp storage used for oligo numbers"
       << endl;
  cerr << "Info: " << ftell(pMatchType_)
       << " bytes of temp storage used for match positions"
       << endl;

  fseek( pOligoNum_, 0, SEEK_SET);
  fseek( pMatchType_, 0, SEEK_SET);
  //  multiPos_.resize(this->matchPosition_.size());
  //  multiType_.resize(this->matchPosition_.size());
  multiMatch_.resize(this->matchPosition_.size());

  //  fclose(pOligoNum_);
  //  fclose(pMatchType_);


  //  if ((pOligoNum_=fopen(tmpFileNameOligoNum_.c_str(), "r"))==NULL)
  //  {
  //   cerr << "Error in MatchTableMulti: could not open file "
  //	 << tmpFileNameOligoNum_ << endl;
  //  exit (1);
  // }

  //  if ((pMatchType_=fopen(tmpFileNameMatchType_.c_str(), "r"))==NULL)
  // {
  //  cerr << "Error in MatchTableMulti: could not open file "
  //	 << tmpFileNameMatchType_ << endl;
  //  exit (1);
  // }



  uint thisOligo;
  uint thisCode;
  MatchPosition thisPos;
  uint numErrors;
  bool reverseFlag;

  while (1)
  {
    if (1 != fread( &thisCode, sizeof(uint), 1, pOligoNum_))
    {
        const int currentError = errno;
        if (feof(pOligoNum_)!=0)
        {
            break;
        }
        else
        {
            BOOST_THROW_EXCEPTION(cc::IoException(currentError, "failed to read Oligo Num"));
        }
    }

    if (1 != fread( &thisPos, sizeof(MatchPosition), 1, pMatchType_))
    {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "failed to read Match Type"));
    }

    numErrors=(thisCode>>30)&0x3;
    reverseFlag=(((thisCode>>29)&0x1)!=0);
    thisOligo=thisCode&(((uint)~0)>>3);
    //    cout << thisCode << " " << thisOligo << " " << (int)reverseFlag << " "<< thisPos << " " << numErrors << endl;
    assert(thisOligo<multiMatch_.size());

  if (   ((numErrors==0)
	  &&(this->matchType_[thisOligo].r[0]<=maxNumMatchesExact_))
      || ((numErrors==1)
	  &&(this->matchType_[thisOligo].r[0]<=maxNumMatchesExact_)
	  &&(this->matchType_[thisOligo].r[1]<=maxNumMatchesOneError_))
      || ((numErrors==2)
	  &&(this->matchType_[thisOligo].r[0]==0)
	  &&(this->matchType_[thisOligo].r[1]<=maxNumMatchesOneError_)
	  &&(this->matchType_[thisOligo].r[2]<=maxNumMatchesTwoErrors_)))
    {
      multiMatch_[thisOligo].push_back( MultiMatch( thisPos, numErrors, 0, reverseFlag) );
    } // ~if
    //    multiMatch_[thisOligo].back().pos_=thisPos;
    //  multiMatch_[thisOligo].back().type_=thisType;

    //    multiPos_[thisOligo].push_back(thisPos);
    //  multiType_[thisOligo].push_back(thisType);
  } // ~while
  //  cout << "READ in " << zz << endl;
  //MatchPosition subtractTable[256];
  int chromTable[256];
  int thisChrom(0);
  //int lastChrom;
  uint j;
  uint nbors0,nbors1,nbors2;

  //  uint errorPos1, errorPos2;
  //  uchar shouldBe1, shouldBe2;

  //  bool isUniqueExact, isUniqueError;
  //  char buf[lineLength];
  const char* pOligo;
  //  ushort numErrors;
  char dirChar;// firstN, secondN;
  // Word nInfo;
  MatchPosition thisBlock /*, matchPos */;

  // extract match info
  MatchPosition extractedMatchPos;
  const char* extractedChromName;
  const char* extractedContigName;
  const char* previousChromName(NULL);
  const char* previousContigName(NULL);

  for (vector<MatchPosition>::const_iterator i(blockStarts.begin());
       i!=blockStarts.end()-1;++i)
  {
    for ( j = *i; j != *(i+1); j+=blockSize)
    {
      //subtractTable[j>>blockShift]=*i;
      chromTable[j>>blockShift]=thisChrom;
    } // ~for j
    ++thisChrom;
  } // ~for i

  oligos.rewind();


  //  for (int i(1); i <= numOligos_ ; i++ )
  for (uint i(1); i < this->matchPosition_.size() ; i++ )
  {

    if ( (pOligo=oligos.getNextOligoSelect(true,false)) == NULL )
    {
      cerr << "OligoInfo: unexpectedly ran out of names!" << endl;
      exit (1);
    } // ~if

    fprintf( this->pOut_, "%s\t%s\t", oligos.getLastName(), pOligo );

    if (this->matchPosition_[i]>=blockRepeat)
    {
      if (this->matchPosition_[i]==qualityFailed)
	fprintf( this->pOut_, "QC");
      else if (this->matchPosition_[i]==repeatMasked)
	fprintf( this->pOut_,"RM");
      else
	fprintf( this->pOut_,"RB\t%u", this->matchPosition_[i]-blockRepeat);
    } // ~if
    else
    {
      numErrors=(this->matchType_[i].errorType&0x3);
      nbors0=((numErrors==0)?((uint)this->matchType_[i].r[0]):0);
      nbors1=((numErrors<=1)?((uint)this->matchType_[i].r[1]):0);
      nbors2=(uint)this->matchType_[i].r[2];

      if ((nbors0==0)&&(nbors1==0)&&(nbors2==0))
      {
	fprintf( this->pOut_, "NM" );
      } // ~if
      else
      {
	fprintf( this->pOut_,"%u:%u:%u",
		 nbors0,nbors1,nbors2 );
      } // ~else
      //	     ((numErrors==0)?((uint)matchType_[i].r[0]):0),
      //	     ((numErrors<=1)?((uint)matchType_[i].r[1]):0),
      //	     (uint)matchType_[i].r[2] );
      sort(multiMatch_[i].begin(),multiMatch_[i].end());

      //lastChrom=-1;
      previousChromName=NULL;

      for (uint j(0); j<multiMatch_[i].size(); ++j)
      {
          dirChar = multiMatch_[i][j].reverse_ ? 'R' : 'F';
          numErrors = multiMatch_[i][j].errors_;
          thisBlock = multiMatch_[i][j].pos_ >> blockShift;
          //matchPos=multiMatch_[i][j].pos_-subtractTable[thisBlock];
          thisChrom=chromTable[thisBlock];

        getMatchPos( multiMatch_[i][j].pos_,
		     extractedChromName,
		     extractedContigName,
		     extractedMatchPos);

	//	if(extractedMatchPos!=matchPos)
	//	{
	//	  cerr << extractedMatchPos << " " << matchPos << endl;
	//	}


	//	if (thisChrom!=lastChrom)
	// NB we are comparing pointers not strings here - OK if strings are
	// always in same location
	if ((extractedChromName!=previousChromName)
	    ||(extractedContigName!=previousContigName))
	{
	  //	  if (lastChrom!=-1)
	  if (previousChromName!=NULL)
	  {
	    fprintf( this->pOut_, "," );
	  }
	  else fprintf( this->pOut_,"\t");
	  //	  fprintf( this->pOut_, "%s:",chromNames[chromTable[thisBlock]].c_str());
	  fprintf( this->pOut_, "%s%s:",extractedChromName, extractedContigName);
	  //lastChrom=thisChrom;
	  previousChromName=extractedChromName;
	  previousContigName=extractedContigName;
	}
	else
	{
	  fprintf( this->pOut_, "," );
	} // ~else
	//	adjustMatchPos( pOligo, dirChar, firstN, secondN, matchPos );

	fprintf( this->pOut_,"%u%c%u",
		 extractedMatchPos,
		 dirChar,
		 numErrors );
	//	fprintf( this->pOut_,"%u%c%u",
	// matchPos,//this->matchPosition_[i]-subtractTable[thisBlock],
	// dirChar,
	// numErrors );
      } // ~for

    } // ~else
    if (0 > fprintf( this->pOut_,"\n"))
    {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "failed to write ELAND output"));
    }
  } // ~for i

} // ~void MatchTableMulti::print



/// Building up multiMatch_, code moved from printSquash
bool MatchTableMulti::buildMatchTable(MatchPositionTranslator& getMatchPos)
{
  cerr << "Info: " << matchesStored_
       << " matches were stored" << endl;
  cerr << "Info: " << ftell(pOligoNum_)
       << " bytes of temp storage used for oligo numbers"
       << endl;
  cerr << "Info: " << ftell(pMatchType_)
       << " bytes of temp storage used for match positions"
       << endl;


  fseek( pOligoNum_, 0, SEEK_SET);
  fseek( pMatchType_, 0, SEEK_SET);

  // firing up the state machinery
  uint table_size = this->matchPosition_.size();
  vector<int> seedOffsets(4,0);

  // if we haven't set the read_length_ accordingly, the assert statements below will fail in singleseed mode
  if( this->read_length_ == 0 )
  {
      this->read_length_ = OLIGO_LEN_;
  }

#ifdef DEBUG
  cerr << "buildMatchTable: hyperhyper_ = " << hyperhyper_.size() << endl;
#endif

  // we don't want to fill the seedOffsets_ array for the
  // MatchTableMulti object from the first pass
  // if size of unmapped_ or hyperhyper_ is == 0, then we are in the multiseed table.
  if( this->hyperhyper_.empty())
  {
      vector<int> so(calculateSeedOffsets(OLIGO_LEN_,this->read_length_));
      using std::swap; swap(seedOffsets, so);
  }

  StateMachine seedsToMatch;
  seedsToMatch.initialize( table_size, seedOffsets );


  uint thisOligo;
  uint thisCode;
  MatchPosition thisPos;
  uint numErrors;
  //bool reverseFlag;
  //uchar thisType;
  //uchar seedNo;

  int oligonums_read = 0;

  int update_cnt = 0;

  vector< vector<bool> > touched( 3, vector<bool>( table_size,false ) );


  while (1)
  {
      oligonums_read++;
      if (1 != fread( &thisCode, sizeof(uint), 1, pOligoNum_))
      {
          const int currentError = errno;
          if (feof(pOligoNum_)!=0)
          {
              break;
          }
          else
          {
              BOOST_THROW_EXCEPTION(cc::IoException(currentError, "failed to read Oligo Num"));
          }
      }

      if (1 != fread( &thisPos, sizeof(MatchPosition), 1, pMatchType_))
      {
          BOOST_THROW_EXCEPTION(cc::IoException(errno, "failed to read Match Type"));
      }

      numErrors=(thisCode>>30)&0x3;
      //reverseFlag=(((thisCode>>29)&0x1)!=0);
      //seedNo=(thisCode>>27)&0x3;
      //thisType=((uchar)numErrors)|(((uchar)reverseFlag)<<7);

      // we now have to shift it by 5, since we also added the seed bits
      //    thisOligo=thisCode&(((uint)~0)>>3);
      thisOligo=thisCode&(((uint)~0)>>5);

      //    cout << thisCode << " " << thisOligo << " " << (int)reverseFlag << " "<< thisPos << " " << numErrors << endl;
      assert(thisOligo<table_size);
      //    assert( thisOligo<matchPosition_.size() );

      // NB:
      // check hyper match; if hyperhyper_.size() > 0, then we are
      // dealing with the object from tier 1. If it was tier 2, then
      // hyperhyper_.size would be 0.
      //
      // If we are in tier 1, we are building the match table now, but
      // we don't want to fill any hypermatches here, because if there
      // were any hits, then we would have filled them up in tier
      // 2. When merging the two match tables, we are setting
      // matchType_ back to 0/0/0, we would insert any hypermatches
      // from tier 1 that were saved to disk.
      bool do_insert = true;
      if( hyperhyper_.size()>0 ) { if( hyperhyper_[thisOligo]==true ) { do_insert=false; } }

      if( do_insert && (   ((numErrors==0)
              &&(this->matchType_[thisOligo].r[0]<=maxNumMatchesExact_))
             || ((numErrors==1)
                 &&(this->matchType_[thisOligo].r[0]<=maxNumMatchesExact_)
                 &&(this->matchType_[thisOligo].r[1]<=maxNumMatchesOneError_))
             || ((numErrors==2)
                 &&(this->matchType_[thisOligo].r[0]==0)
                 &&(this->matchType_[thisOligo].r[1]<=maxNumMatchesOneError_)
                 &&(this->matchType_[thisOligo].r[2]<=maxNumMatchesTwoErrors_))) )
      {
          //      multiMatch_[thisOligo].push_back( MultiMatch( thisPos, thisType ) );
          seedsToMatch.insertSeedHit( thisPos,thisCode,update_cnt,getMatchPos );
          touched[numErrors][thisOligo] = true;
      }
  } // ~while

  // resize the multiMatch_ table
  if( multiMatch_.size() == 0 )
      multiMatch_.resize( table_size );

  // building up the multiMatch table
  for( uint i=1;i<table_size;i++ )
  {
      //      list<MultiMatch> l = seedsToMatch.getHits(i,maxNumMatchesExact_);
      list<MultiMatch> l = seedsToMatch.getHits(i,100, this->hyperhyper_.empty()); // 100 is only an arbitrary high value such that we do not cut any (had the check above)

      if( l.size()>0 )
      {
          multiMatch_[i].insert( multiMatch_[i].end(), l.begin(),l.end() );
          seedsToMatch.clear(i);
      }

      seedsToMatch.matchType_[i].errorType = this->matchType_[i].errorType;

      if( touched[0][i] == false )
      {
          seedsToMatch.matchType_[i].r[0] = this->matchType_[i].r[0];
      }
      if( touched[1][i] == false )
      {
          seedsToMatch.matchType_[i].r[1] = this->matchType_[i].r[1];
      }
      if( touched[2][i] == false )
      {
          seedsToMatch.matchType_[i].r[2] = this->matchType_[i].r[2];
      }



      // in the case that we did a merge with another table, we do not want to overwrite the entries;
      // check if the read that we want to merge is unmapped!
      if( unmapped_.size() > 0 )
      {
          // this branch is true iff the function getUnmappedReads was
          // called, ie this is the MatchTableMulti from the first of
          // the two passes; in this buildMatchTable we are only
          // interested in the reads that *got* mapped in the first
          // pass, because the unmapped reads were dealt with in the
          // mergeTable step! --> therefore we have to check for
          // unmapped_[i]==false
          if( unmapped_[i] == false )
          {
//              matchType_[i] = seedsToMatch.matchType_[i];
              // keep the original matchType_[i] information!
          }
          // in fact, we don't have to do anything in this branch of
          // the code, because for the first MatchTableMulti object we
          // already filled the matchType_ array accordingly
      } else
      {
        this->matchType_[i] = seedsToMatch.matchType_[i];
      }


      if( (this->matchType_[i].r[0]==0)&&
          (this->matchType_[i].r[1]==0)&&
          (this->matchType_[i].r[2]==0) )
      {

          assert(multiMatch_[i].size()==0);
      }

  }

  seedsToMatch.clear();





  return true;
}


void MatchTableMulti::printSquash( OligoSource& oligos,
				   MatchPositionTranslator& getMatchPos,
				   const vector<string>& ,
				   const vector<MatchPosition>& /*blockStarts*/,
				   const SuffixScoreTable& ,
				   int oligoLength,
				   const string& directoryName,
				   const bool& align )
{
  ofstream match_out( this->outputFileName_.c_str() );
  oligos.rewind();
  const char* strlen_oligo;
  // first, deduce the oligo length (nb: the parameter oligoLength contains the length of the seed)
  if ( (strlen_oligo=oligos.getNextOligoSelect(true,false)) == NULL )
    {
      cerr << "printSquash: no results to print as there was no data to align." << endl;
      return; // allow the pipeline to continue
    } // ~if
  int readLength = strlen(strlen_oligo);
  // if we have indels within the read (at most ALIGN_DP_BAND nucleotides)
  int fragmentLength = strlen(strlen_oligo) + ALIGN_DP_BAND;
  // rewind selector
  oligos.rewind();


  // open temporary stream for writing the match locations
  //  hits_out.open("hits.txt");


//  const int reverseStrandStartOffset(readLength-no_of_seeds_*oligoLength);
  const int reverseStrandStartOffset(readLength-oligoLength);
  FragmentFinder getFragments(directoryName,
                              readLength,
                              fragmentLength,
                              reverseStrandStartOffset);


  StringIndex files(directoryName);


  // setting up the gap settings for ELAND
  ifstream align_score(".align.scores");
  // the actual score values

#ifdef READ_PARAMETERS
  if( !align_score ) {
    if( aligner.readAlignScoreFile(align_score,match,mismatch,gapopen,gapextend) == true )
      {
	cerr << "read .align.scores, setting alignment scores to match=" << match
	     << ",mismatch=" << mismatch
	     << ",gapopen=" << gapopen
	     << ",gapextend=" << gapextend << endl;
      }
  } // ~if(!align.score)
#endif


  // code moved from below
  // building up the blockStarts and subtractTable
  // moved code starts here...
  MatchPosition extractedMatchPos;
  const char* extractedChromName;
  const char* extractedContigName;
  const char* previousChromName(NULL);
  const char* previousContigName(NULL);

  uint nbors0,nbors1,nbors2;

  uint numErrors;

  // start building the MatchTableMulti object
  buildMatchTable(getMatchPos);

  const char* pOligo;

  oligos.rewind();

  vector< MatchRequest > matches;
  vector<SeqRequest> frag_requests;
  vector<char*> reads;

  uint request_cnt = 0;
  uint mr_cnt = 0;

  vector<int> allSeedOffsets(calculateSeedOffsets(OLIGO_LEN_,readLength));
  allSeedOffsets.insert(allSeedOffsets.begin(), 0);
  std::cerr << "all seed offsets: " << allSeedOffsets[0] << " "
      << allSeedOffsets[1] << " " << allSeedOffsets[2] << " "
      << allSeedOffsets[3] << " " << allSeedOffsets[4] << std::endl;


  //  for (int i(1); i <= numOligos_ ; i++ )
  for (uint i(1); i < this->matchPosition_.size() ; i++ )
  {

    // ----------------------------------------------------------------------
    // PULLING OUT THE FRAGMENTS:
    // ok, check if the number of fragments that we want to pull out
    // from the genome is greater than REQUEST_SIZE, in this case get
    // the fragments and print everything to pOut_; if not, then go
    // ahead and add new matches to
    if( (request_cnt > REQUEST_SIZE) || (mr_cnt > MR_REQUEST_SIZE) )
      {

	if( unsquashRequests( request_cnt,
			      match_out,
			      getFragments,
			      align,
			      files,
			      matches,
			      frag_requests,
			      reads,
			      readLength,
			      fragmentLength
			      ) == false )
	  {
	    // do nothing for the moment
	  }

	// reset the request counter, clear the vectors
	request_cnt = 0;
	mr_cnt = 0;

	{
	  vector<MatchRequest> mr_tmp;
	  matches.swap(mr_tmp);
	  matches.clear();

	  vector<SeqRequest> sr_tmp;
	  frag_requests.swap(sr_tmp);
	  frag_requests.clear();
	} //~extra brackets for scope
      }



    if ( (pOligo=oligos.getNextOligoSelect(true,false)) == NULL )
      {
	cerr << "OligoInfo: unexpectedly ran out of names!" << endl;
	exit (1);
      } // ~if


    MatchRequest cur_mr;

    cur_mr.header_ = string(oligos.getLastName());
    cur_mr.read_   = string(pOligo);
    cur_mr.matchMode_ = -1;

    reads.push_back( new char[readLength+1] );
    strncpy( reads[ reads.size()-1 ],pOligo,readLength );
    reads[ reads.size()-1 ][readLength] = '\0';

    fprintf( this->pOut_, "%s\t%s\t", oligos.getLastName(), pOligo );

    if (this->matchPosition_[i]>=blockRepeat)
    {
        if (this->matchPosition_[i]==qualityFailed)
        {
            cur_mr.matchMode_ = 0;
            fprintf( this->pOut_, "QC");
        }
        else if (this->matchPosition_[i]==repeatMasked)
        {
            cur_mr.matchMode_ = 1;
            fprintf( this->pOut_,"RM");
        }
        else
        {
            cur_mr.matchMode_ = 2;
            cur_mr.rb_position_ = this->matchPosition_[i]-blockRepeat;
            fprintf( this->pOut_,"RB\t%u", this->matchPosition_[i]-blockRepeat);
        }
    }
    else
    {
        cur_mr.matchMode_ = 3;

        numErrors=(this->matchType_[i].errorType&0x3);

        nbors0=((numErrors==0)?((uint)this->matchType_[i].r[0]):0);
        nbors1=((numErrors<=1)?((uint)this->matchType_[i].r[1]):0);
        nbors2=(uint)this->matchType_[i].r[2];

        if( (nbors0==0)&&(nbors1==0)&&(nbors2==0) )
        {
            assert(multiMatch_[i].size()==0);
        }

        if ((nbors0==0)&&(nbors1==0)&&(nbors2==0))
        {
            cur_mr.nbors0_ = 0;
            cur_mr.nbors1_ = 0;
            cur_mr.nbors2_ = 0;

            fprintf( this->pOut_, "NM" );
        } // ~if
        else
        {
            cur_mr.nbors0_ = nbors0;
            cur_mr.nbors1_ = nbors1;
            cur_mr.nbors2_ = nbors2;


            fprintf( this->pOut_,"%u:%u:%u",
                     nbors0,nbors1,nbors2 );
        } // ~else
        //	     ((numErrors==0)?((uint)this->matchType_[i].r[0]):0),
        //	     ((numErrors<=1)?((uint)this->matchType_[i].r[1]):0),
        //	     (uint)this->matchType_[i].r[2] );
        sort(multiMatch_[i].begin(),multiMatch_[i].end());

        previousChromName=NULL;

 
        for (uint j(0); j<multiMatch_[i].size(); ++j)
        {
            const char dirChar=multiMatch_[i][j].reverse_ ? 'R':'F';

            // offset of the matched seed from the seed 0
            const MatchPosition seedOffset = allSeedOffsets[multiMatch_[i][j].lastSeed_];
            MatchPosition seedPos(multiMatch_[i][j].pos_);

            // determine the number of leading Ns that overlap the first seed
            unsigned int noLeadingNs = 0;
            while( (cur_mr.read_.length() > (noLeadingNs + seedOffset)) && ('N' == cur_mr.read_[noLeadingNs + seedOffset] || 'n' == cur_mr.read_[noLeadingNs + seedOffset]) )
            {
                noLeadingNs++;
            }

            if (multiMatch_[i][j].reverse_)
            {
                seedPos -= seedOffset;
            }
            else
            {
                seedPos += seedOffset;
            }
//            std::cerr << "extracting " << seedPos << " seed " << (int)multiMatch_[i][j].lastSeed_
//                << " seedOffset " << (int)seedOffset
//                << " reverse " << (int)multiMatch_[i][j].reverse_
//                << " errors " << (int)multiMatch_[i][j].errors_
//                << " read " << cur_mr.read_ << std::endl;
            getMatchPos( seedPos, extractedChromName, extractedContigName, extractedMatchPos);
//            std::cerr << "extracted " << seedPos << " " << extractedChromName
//                << " " << extractedContigName << " " <<extractedMatchPos << std::endl;
            // now the extractedMatchPos contains the position of the matched seed in the
            // corresponding contig. Adjust it back to be in respect to the seed 0 position
            if (multiMatch_[i][j].reverse_)
            {
                extractedMatchPos += seedOffset;
            }
            else
            {
                extractedMatchPos -= seedOffset;
            }

            // NB we are comparing pointers not strings here - OK if strings are
            // always in same location
            if ((extractedChromName!=previousChromName) || (extractedContigName!=previousContigName))
            {
                if (previousChromName!=NULL)
                {
                    fprintf( this->pOut_, "," );
                }
                else fprintf( this->pOut_,"\t");
                fprintf( this->pOut_, "%s%s:",extractedChromName, extractedContigName);
                previousChromName=extractedChromName;
                previousContigName=extractedContigName;

                // store the extracted chromosome and contig names in MatchRequest
                cur_mr.chromNames_.push_back( string( string(extractedChromName)+string(extractedContigName)) );

                vector<HitPosition> hp_tmp;
                cur_mr.hits_.push_back( hp_tmp );
            }
            else
            {
                fprintf( this->pOut_, "," );
            }

            // create an entry on the hits_
            HitPosition cur_hit;
            cur_hit.matchPosition_ = extractedMatchPos;

            if( dirChar == 'R' )
            {
                cur_hit.matchPosition_ = cur_hit.matchPosition_  + noLeadingNs - reverseStrandStartOffset;
            }
            else
            {
                cur_hit.matchPosition_ -= noLeadingNs;
            }

            cur_hit.direction_ = dirChar;
            cur_hit.numErrors_ = numErrors;
            cur_mr.hits_[ cur_mr.hits_.size()-1 ].push_back( cur_hit );

            uint chromNum = 0;
            uint contigNum = 0;
            uint pos = 0;


            // retrieve the chromosome number and the contig number
            files.getIndex( cur_mr.chromNames_[cur_mr.chromNames_.size()-1].c_str() ,
                            chromNum,
                            contigNum,
                            pos );

            // create the corresponding SeqRequest
            // NB markus it's frag_requests.size() , because we push the SeqRequest afterwards,
            // but we already pushed the read, so it's reads.size()-1
            frag_requests.push_back( SeqRequest(frag_requests.size(),
                                                reads.size()-1,
                                                chromNum,
                                                contigNum,
                                                (pos+extractedMatchPos),
                                                dirChar,
                                                seedOffset)
                                   );

            fprintf( this->pOut_,"%u%c%u",
                     extractedMatchPos,
                     dirChar,
                     numErrors );
            request_cnt++;
        }

    }

    matches.push_back( cur_mr );
    mr_cnt++;
    if (0 > fprintf( this->pOut_,"\n"))
    {
        BOOST_THROW_EXCEPTION(cc::IoException(errno, "failed to write ELAND output"));
    }

  } // ~for i

  cerr << "REQUEST_CNT = " << request_cnt << endl;


  if( unsquashRequests( request_cnt,
			match_out,
			getFragments,
			align,
			files,
			matches,
			frag_requests,
			reads,
			readLength,
			fragmentLength
			) == false )
    {
      // do nothing for the moment
    }


  match_out.close();
}

// Clear the matchPosition_ and matchType_ vectors
bool MatchTableMulti::clear(void)
{
    // clearing vectors
    vector< vector<MultiMatch> > tmp_mm;
    vector<MatchDescriptor> tmp_md;

    multiMatch_.clear();
    this->matchType_.clear();

    this->multiMatch_.swap(tmp_mm);
    this->matchType_.swap(tmp_md);

    return true;
}


// Given another MatchTableMulti object (most likely the one from a
// second pass through the reads using multiseed mode, merge the
// results from both objects
bool MatchTableMulti::mergeTable( MatchTable* source,MatchPositionTranslator& getMatchPos )
{
    MatchTableMultiSquareSeed* source_multi = static_cast<MatchTableMultiSquareSeed*>(source);

    vector< vector<MultiMatch> > multimatches_second_tier;
    vector<MatchDescriptor> matchdescriptor_second_tier;

    // build up the match table
    source_multi->buildMatchTable(getMatchPos);


    if( source_multi->getMatchInformation(multimatches_second_tier,
                                          matchdescriptor_second_tier) == false )
    {
        return false;
    } else
    {
        // clear the internal vectors, we don't need them anymore here
        source_multi->clear();

        if( multiMatch_.size()==0 )
            multiMatch_.resize( this->matchPosition_.size());

        // merge the information obtained from the multiseed pass into the current object
        for( uint i=1;i<multimatches_second_tier.size();i++ )
        {
            if( matchdescriptor_second_tier[i].r[0] == 0 &&
                matchdescriptor_second_tier[i].r[1] == 0 &&
                matchdescriptor_second_tier[i].r[2] == 0 )
            {
                assert( multimatches_second_tier[i].size() == 0 );
            }


            if( multimatches_second_tier[i].size() > 0 )
            {
                this->multiMatch_[ this->translator_[i] ].insert( this->multiMatch_[this->translator_[i]].end(),
                                                      multimatches_second_tier[i].begin(),
                                                      multimatches_second_tier[i].end() );
            }
            this->matchType_[ this->translator_[i] ].errorType = matchdescriptor_second_tier[i].errorType;
            this->matchType_[ this->translator_[i] ].r[0] = matchdescriptor_second_tier[i].r[0];
            this->matchType_[ this->translator_[i] ].r[1] = matchdescriptor_second_tier[i].r[1];
            this->matchType_[ this->translator_[i] ].r[2] = matchdescriptor_second_tier[i].r[2];

            if( this->matchType_[this->translator_[i]].r[0] == 0 &&
                this->matchType_[this->translator_[i]].r[1] == 0 &&
                this->matchType_[this->translator_[i]].r[2] == 0 )
            {
                assert( this->multiMatch_[this->translator_[i]].size() == 0 );
            }

        }

    }

    // clearing vectors
    vector< vector<MultiMatch> > tmp_mm;
    vector<MatchDescriptor> tmp_md;
    vector<uint> tmp_translator;

    multimatches_second_tier.clear();
    matchdescriptor_second_tier.clear();
    this->translator_.clear();

    multimatches_second_tier.swap(tmp_mm);
    matchdescriptor_second_tier.swap(tmp_md);
    this->translator_.swap(tmp_translator);


    return true;
}


// retrieve the match information from the second pass to merge it
// with the MatchTableMulti object of the singleseed run
bool MatchTableMulti::getMatchInformation( vector< vector<MultiMatch> >& multimatches,vector<MatchDescriptor>& matchdescriptor )
{
    // fill the return objects
    multimatches = multiMatch_;
    matchdescriptor = this->matchType_;
//    matchdescriptor = seedsToMatch.matchType_;


    return true;
}


// ---------------------- extending MatchTableMulti to multiple seeds --------------------------
void MatchTableMultiSquareSeed::
addMatch(MatchIter i, const MatchIter i_end)
{
  static const OligoNumber onMask((~isReverseOligo)&(~seed_bits[3]));

  for(;i!=i_end;++i) {

      assert(i->numErrors<=2);

  // For oligos with no Ns:
  // On partition 0, interested in 0, 1 and 2 error matches
  // On partition 1, interested in 0 or 1 error matches
  // On partitions 2 to 5, interested in 2 error matches only

      const OligoNumber oligoNum(i->position&onMask);

      const uint8_t seedNo((i->position&(~isReverseOligo))>>29);

      ms_matchType_[4*oligoNum+seedNo].r[i->numErrors]
          +=(ms_matchType_[4*oligoNum+seedNo].r[i->numErrors]!=0xFF);


      if( checkNumberOfHits( oligoNum,seedNo,i->numErrors ) ) {
          // TBD: cache then bulk write
          const uint32_t matchCode((i->numErrors<<30)|
                     (((i->position&isReverseOligo)!=0)<<29)|
                     (seedNo<<27)|oligoNum);
          ++this->matchesStored_;

          if( 1 != fwrite (&matchCode,sizeof(uint),1,this->pOligoNum_) )
          {
              cerr << "error writing match code to temp file." << endl;
          }
          if( 1 != fwrite (&(i->matchPos),sizeof(MatchPosition),1,this->pMatchType_) )
          {
              cerr << "error writing match pos to temp file." << endl;
          }
      }
  }

#if (NUM_THREADS>1)
      pthread_mutex_unlock(&mutex_);
#endif
}

//void MatchTableMultiSquareSeed::print( OligoSource& oligos,
//                                 MatchPositionTranslator& getMatchPos,
//                                 const vector<string>& chromNames,
//                                 const vector<MatchPosition>& blockStarts,
//                                 const SuffixScoreTable& scoreTable,
//                                 int oligoLength)
//{
//    this->MatchTableMulti::print(oligos,getMatchPos,chromNames,blockStarts,scoreTable,oligoLength);
//}
//

//void MatchTableMultiSquareSeed::printSquash( OligoSource& oligos,
//                                       MatchPositionTranslator& getMatchPos,
//                                       const vector<string>& chromNames,
//                                       const vector<MatchPosition>& blockStarts,
//                                       const SuffixScoreTable& scoreTable,
//                                       int oligoLength,
//                                       const string& directoryName,
//                                       const bool& align)
//{
//    this->MatchTableMulti::printSquash(oligos,getMatchPos,chromNames,blockStarts,scoreTable,oligoLength,directoryName,align);
//}
//
//bool MatchTableMultiSquareSeed::getUnmappedReads( vector<bool>& unmapped )
//{
//    return this->MatchTableMulti::getUnmappedReads(unmapped);
//}

//bool MatchTableMultiSquareSeed::getMatchInformation( vector< vector<MultiMatch> >& multimatches,vector<MatchDescriptor>& matchdescriptor )
//{
//    return this->MatchTableMulti::getMatchInformation(multimatches,matchdescriptor);
//}
bool MatchTableMultiSquareSeed::mergeTable( MatchTable* source,MatchPositionTranslator& getMatchPos )
{
    MatchTableMultiSquareSeed* source_multiseed = static_cast<MatchTableMultiSquareSeed*>(source);

    vector< vector<MultiMatch> > multimatches_second_tier;
    vector<MatchDescriptor> matchdescriptor_second_tier;

    // build up the match table
    source_multiseed->buildMatchTable(getMatchPos);


    if( source_multiseed->getMatchInformation(multimatches_second_tier,
                                          matchdescriptor_second_tier) == false )
    {
        return false;
    } else
    {
        // clear the internal vectors, we don't need them anymore here
        source_multiseed->clear();

        if( this->multiMatch_.size()==0 )
          this->multiMatch_.resize( this->matchPosition_.size());

        // merge the information obtained from the multiseed pass into the current object
        for( uint i=1;i<multimatches_second_tier.size();i++ )
        {
            if( matchdescriptor_second_tier[i].r[0] == 0 &&
                matchdescriptor_second_tier[i].r[1] == 0 &&
                matchdescriptor_second_tier[i].r[2] == 0 )
            {
                assert( multimatches_second_tier[i].size() == 0 );
            }


            if( multimatches_second_tier[i].size() > 0 )
            {
                this->multiMatch_[ this->translator_[i] ].insert( this->multiMatch_[this->translator_[i]].end(),
                                                      multimatches_second_tier[i].begin(),
                                                      multimatches_second_tier[i].end() );
            }
            this->matchType_[ this->translator_[i] ].errorType = matchdescriptor_second_tier[i].errorType;
            this->matchType_[ this->translator_[i] ].r[0] = matchdescriptor_second_tier[i].r[0];
            this->matchType_[ this->translator_[i] ].r[1] = matchdescriptor_second_tier[i].r[1];
            this->matchType_[ this->translator_[i] ].r[2] = matchdescriptor_second_tier[i].r[2];

            if( this->matchType_[this->translator_[i]].r[0] == 0 &&
                this->matchType_[this->translator_[i]].r[1] == 0 &&
                this->matchType_[this->translator_[i]].r[2] == 0 )
            {
                assert( this->multiMatch_[this->translator_[i]].size() == 0 );
            }

        }

    }

    // clearing vectors
    vector< vector<MultiMatch> > tmp_mm;
    vector<MatchDescriptor> tmp_md;
    vector<uint> tmp_translator;

    multimatches_second_tier.clear();
    matchdescriptor_second_tier.clear();
    this->translator_.clear();

    multimatches_second_tier.swap(tmp_mm);
    matchdescriptor_second_tier.swap(tmp_md);
    this->translator_.swap(tmp_translator);


    return true;
}

bool MatchTableMultiSquareSeed::buildMatchTable( MatchPositionTranslator& getMatchPos )
{
  cerr << "Info: " << this->matchesStored_
       << " matches were stored" << endl;
  cerr << "Info: " << ftell(this->pOligoNum_)
       << " bytes of temp storage used for oligo numbers"
       << endl;
  cerr << "Info: " << ftell(this->pMatchType_)
       << " bytes of temp storage used for match positions"
       << endl;

  if (!this->matchesStored_)
  {
      return false;
  }

  fseek( this->pOligoNum_, 0, SEEK_SET);
  fseek( this->pMatchType_, 0, SEEK_SET);

  // firing up the state machinery
  uint table_size = this->matchPosition_.size();
  vector<int> seedOffsets(4,0);

  // we don't want to fill the seedOffsets_ array for the
  // MatchTableMulti object from the first pass
  // if unmapped_ or hyperhyper_ is == 0, then we are in the multiseed table.
  if( this->hyperhyper_.empty())
  {
      vector<int> so(calculateSeedOffsets(OLIGO_LEN_,this->read_length_));
      using std::swap; swap(seedOffsets, so);
  }

  StateMachine seedsToMatch;
  seedsToMatch.initialize( table_size, seedOffsets );


  uint thisOligo;
  uint thisCode;
  MatchPosition thisPos;
  uint8_t numErrors;
  //bool reverseFlag;
  //uchar thisType;
  uint8_t seedNo;

  int oligonums_read = 0;

  int update_cnt = 0;

  vector< vector< vector<bool> > > touched( table_size,vector< vector<bool> >( this->no_of_seeds_, vector<bool>(3,false ) ) );


  while (1)
  {
      oligonums_read++;
      if (1 != fread( &thisCode, sizeof(uint), 1, this->pOligoNum_))
      {
          const int currentError = errno;
          if (feof(this->pOligoNum_)!=0)
          {
              break;
          }
          else
          {
              BOOST_THROW_EXCEPTION(cc::IoException(currentError, "failed to read Oligo Num"));
          }
      }

      if (1 != fread( &thisPos, sizeof(MatchPosition), 1, this->pMatchType_))
      {
          BOOST_THROW_EXCEPTION(cc::IoException(errno, "failed to read Match Type"));
      }

      numErrors=(thisCode>>30)&0x3;
      //reverseFlag=(((thisCode>>29)&0x1)!=0);
      seedNo=(thisCode>>27)&0x3;
      //thisType=((uchar)numErrors)|(((uchar)reverseFlag)<<7);

      // we now have to shift it by 5, since we also added the seed bits
      //    thisOligo=thisCode&(((uint)~0)>>3);
      thisOligo=thisCode&(((uint)~0)>>5);

//      cout << thisCode << " " << thisOligo << " " << (int)reverseFlag << " "<< thisPos << " " << numErrors << endl;
      assert(thisOligo<table_size);
      //    assert( thisOligo<this->matchPosition_.size() );

      // NB:
      // check hyper match; if hyperhyper_.size() > 0, then we are
      // dealing with the object form tier 1. If it was tier 2, then
      // hyperhyper_.size would be 0.
      //
      // If we are in tier 1, we are building the match table now, but
      // we don't want to fill any hypermatches here, because if there
      // were any hits, then we would have filled them up in tier
      // 2. When merging the two match tables, we are setting
      // this->matchType_ back to 0/0/0, we would insert any hypermatches
      // from tier 1 that were saved to disk.
      bool do_insert = true;
      if( this->hyperhyper_.size()>0 ) { if( this->hyperhyper_[thisOligo]==true ) { do_insert=false; } }

      if( do_insert && (checkNumberOfHits(thisOligo,seedNo,numErrors)) )
      {
          //      this->multiMatch_[thisOligo].push_back( MultiMatch( thisPos, thisType ) );
          seedsToMatch.insertSeedHit( thisPos,thisCode,update_cnt,getMatchPos );
          touched[thisOligo][seedNo][numErrors] = true;
      }
  } // ~while

  // resize the this->multiMatch_ table
  if( this->multiMatch_.size() == 0 )
      this->multiMatch_.resize( table_size );

  // building up the multiMatch table
  for( uint i=1;i<table_size;i++ )
  {
      //      list<MultiMatch> l = seedsToMatch.getHits(i,maxNumMatchesExact_);
      list<MultiMatch> l = seedsToMatch.getHits(i,10000, this->hyperhyper_.empty()); // 10000 is only an arbitrary high value such that we do not cut any (had the check above)

      // TODO: MB 09/11/12
      //
      // parse the hits, if all the hits have the same number of seed
      // hits, then check if the number of hits < default_number; if
      // some hits have more seed matches than the others, only take
      // the hits having showing the most seed hits
      if( l.size()>0 )
      {
          this->multiMatch_[i].insert( this->multiMatch_[i].end(), l.begin(),l.end() );
          seedsToMatch.clear(i);
      }

      seedsToMatch.matchType_[i].errorType = ms_matchType_[i].errorType;

      // iterate over touched, because now we have it various
      vector<bool> cumulative_touched(3,false);

      for( short j=0;j<3;j++ )
          for( short k=0;k<this->no_of_seeds_;k++ )
          {
              // cumulative touch for 0/1/2 error matches
              cumulative_touched[j] = ( cumulative_touched[j] || touched[i][k][j] );
          }

      if( cumulative_touched[0] == false )
      {
          const int cumulativeSum = (ms_matchType_[4*i].r[0]+ms_matchType_[4*i+1].r[0]+ms_matchType_[4*i+2].r[0]+ms_matchType_[4*i+3].r[0]);

          seedsToMatch.matchType_[i].r[0] = (cumulativeSum<=255)?cumulativeSum:255;
      }
      if( cumulative_touched[1] == false )
      {
          const int cumulativeSum = (ms_matchType_[4*i].r[1]+ms_matchType_[4*i+1].r[1]+ms_matchType_[4*i+2].r[1]+ms_matchType_[4*i+3].r[1]);

          seedsToMatch.matchType_[i].r[1] = (cumulativeSum<=255)?cumulativeSum:255;
      }
      if( cumulative_touched[2] == false )
      {
          const int cumulativeSum = (ms_matchType_[4*i].r[2]+ms_matchType_[4*i+1].r[2]+ms_matchType_[4*i+2].r[2]+ms_matchType_[4*i+3].r[2]);

          seedsToMatch.matchType_[i].r[2] = (cumulativeSum<=255)?cumulativeSum:255;
      }



      // in the case that we did a merge with another table, we do not want to overwrite the entries;
      // check if the read that we want to merge is unmapped!
      if( this->unmapped_.size() > 0 )
      {
          // this branch is true iff the function getUnmappedReads was
          // called, ie this is the MatchTableMulti from the first of
          // the two passes; in this buildMatchTable we are only
          // interested in the reads that *got* mapped in the first
          // pass, because the unmapped reads were dealt with in the
          // mergeTable step! --> therefore we have to check for
          // unmapped_[i]==false
      } else
      {
          this->matchType_[i] = seedsToMatch.matchType_[i];
      }


      if( (this->matchType_[i].r[0]==0)&&
          (this->matchType_[i].r[1]==0)&&
          (this->matchType_[i].r[2]==0) )
      {

          assert(this->multiMatch_[i].size()==0);
      }

  }

  seedsToMatch.clear();





  return true;
}


} //namespace eland_ms
} //namespace casava
// end of ELAND_outer.cpp
