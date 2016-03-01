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
 ** \file eland_ms/MatchRequest.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 ** Contains interface functionality from ELAND - specifically, stuff to do
 ** with scoring and storage of alignments
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_MATCH_REQUEST_H
#define CASAVA_ELAND_MS_MATCH_REQUEST_H

namespace casava
{
namespace eland_ms
{
// 09/04/09 markus
// define MatchRequest and Hit here
//
// Hit:
// - uint extractedMatchPos
// - char dirChar
// - short numErrors
//
// MatchRequest:
// - string header
// - string read
// - char matchPosition
// - short nbors0, short nbors1, short nbors2
// - vector< vector< Hit > >

// markus todo
// - add constructor
struct HitPosition {
  // has to be long and not MatchPosition, because we can end up having negative alignment positions!
  long matchPosition_;
  char direction_;
  int numErrors_;
};

// markus todo
// - add constructor
struct MatchRequest {
  string header_;
  string read_;
  // 0 = QC , 1 = RM , 2 = RB (need to store another value) , 3 = full alignment
  short matchMode_;
  uint rb_position_;
  short nbors0_,nbors1_,nbors2_;
  vector< string > chromNames_; // extracted chromosome names
  vector< vector< HitPosition > > hits_; // the actual hits

  // print the information to out
    void print( ostream& out,vector<char*>& frags,int& frag_idx,const vector<int>& pos_correction_begin,const vector<int>& pos_correction_end )
    {
        out << header_ << "\t"
            << read_ << "\t";

        switch( matchMode_ ) {
        case 0:
            out << "QC\t-";
            break;
        case 1:
            out << "RM\t-";
            break;
        case 2:
            out << "RB\t" << rb_position_;
            break;
        case 3:
            if ((nbors0_==0)&&(nbors1_==0)&&(nbors2_==0))
            {
                out << "NM\t-";
            } // ~if
            else
            {
                out << nbors0_ << ":"
                    << nbors1_ << ":"
                    << nbors2_;



                // if we have matches to list, then add a tab
                if( chromNames_.size() > 0 )
                {
                    out << "\t";
                    // ok, now print the single hits
                    for( uint i=0;i<chromNames_.size();i++ )
                    {
                        // add proper ',' placement
                        if( i > 0 )
                        {
                            out << ",";
                        }
                        out << chromNames_[i] << ":";

                        for( uint j=0;j<(hits_[i].size()-1);j++ )
                        {
                            if( (pos_correction_begin[frag_idx] != 0) || (pos_correction_end[frag_idx] != 0 ) )
                            {
                                hits_[i][j].matchPosition_ -= ( (hits_[i][j].direction_=='R')?pos_correction_end[frag_idx]:pos_correction_begin[frag_idx] );
                            }

                            out << hits_[i][j].matchPosition_ << hits_[i][j].direction_ << frags[frag_idx++] << ",";
                            //		      frags.erase( frags.begin() );

                        }

                        if( (pos_correction_begin[frag_idx] != 0) || (pos_correction_end[frag_idx] != 0 ) )
                        {
                            hits_[i][ hits_[i].size()-1 ].matchPosition_ -= ( (hits_[i][hits_[i].size()-1].direction_=='R')?pos_correction_end[frag_idx]:pos_correction_begin[frag_idx] );
                        }


                        out << hits_[i][ hits_[i].size()-1 ].matchPosition_
                            << hits_[i][ hits_[i].size()-1 ].direction_
                            << frags[frag_idx++];
                        //		  frags.erase( frags.begin() );
                    }
                }
                else
                {
                    out << "\t-";
                }
            }
            break;
        default:
            cerr << "switch reached default, should not happen.";
            exit(1);
        }


        out << endl;
    }

};


} // namespace eland_ms
} // namespace casava



#endif // CASAVA_ELAND_MS_MATCH_REQUEST_H
