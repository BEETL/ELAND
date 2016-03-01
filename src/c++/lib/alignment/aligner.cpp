// Copyright (c) 2009 Illumina
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).

// PROJECT: ELAND
// MODULE: aligner.cpp
// AUTHOR: M. Bauer

// Description: perform a gapped alignment with reads that have not been mapped onto the genome


#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

#include "alignment/aligner.h"

typedef unsigned char uchar;
typedef unsigned int uint;

#define QUALSCALE 0.0333
#define EPSILON 0.0000001
#define MIN_REPEAT_THRESHOLD 0.01


namespace ga { namespace alignment {


    int Aligner::operator()
        ( const char* qvals_, const char* x, const char* y, const int x_size, const int y_size,const bool& flag )
    {
        // the matrices have to be initialized first
        if( init_ == false ) {
            cerr << "DP matrices not initialized...exiting." << endl;
            return 0;
        }


        //  const ScoreType w_match(3),w_mismatch(1),w_open(2),w_extend(2);
        ScoreType max = -1;
        ScoreType whichMatrix = -1;
        ScoreType thisMax = -1;
        ScoreType thisMatrix = -1;
        ScoreType nextMatrix = -1;

        xt_="";
        yt_="";
        at_.clear();

        int ii = -1;
        int jj = -1;

        for (int i(0);i<=x_size;i++)
        {
            //E_[i][0]=-w_open_-(i*w_extend_);
            E_[i][0]=0;
            F_[i][0]=-10000;
            G_[i][0]=-10000;
        }

        for (int j(0);j<=y_size;j++)
        {
            E_[0][j]=-10000;
            //F_[0][j]=-w_open_-(j*w_extend_);
            F_[0][j]=0;
            G_[0][j]=-10000;
        }



        for (int i(1);i<=x_size;i++)
        {
#ifdef DEBUG
            cout << i << ": ";
#endif
            //      int startIdx = ( (i-width_)< 1 )? 1:(i-width_);
            //      int endIdx = ( (i+width_)> y_size )? y_size:(i+width_);

//	int startIdx = ( (i-width_)< 1 )? 1:(i-width_);
//	int endIdx = ( (i+width_)> y_size )? y_size:(i+width_);

            // we are shifting the DP band to the right 
            int startIdx = i;
            int endIdx = ( (i+2*width_)> y_size )? y_size:(i+2*width_);


            for (int j(startIdx);j<=endIdx;j++)
            {
#ifdef DEBUG
                cout << "i/j" << i << "/" << j << endl;
#endif
                // update G_
                max3(G_[i][j],TG_[i][j],
                     G_[i-1][j-1],E_[i-1][j-1],F_[i-1][j-1],flag);

                //      G_[i][j]=V_[i-1][j-1];
                if (x[i-1]==y[j-1])
                {
                    //	G_[i][j]+=w_match_;

                    // changed to incorporate quality scores
                    G_[i][j]+=(w_match_*(qvals_[i-1]-64)*QUALSCALE);
                }
                else
                {
                    //	G_[i][j]-=w_mismatch_;
                    //	G_[i][j]+=w_mismatch_;

                    // incorporate quality values
                    G_[i][j]+=(w_mismatch_*(qvals_[i-1]-64)*QUALSCALE);
                }
                //      max =G_[i][j];
                //      dir=0;

                // update E_
                max3(E_[i][j],TE_[i][j],
                     G_[i][j-1]-w_open_,
                     E_[i][j-1]-w_extend_,
                     F_[i][j-1]-w_open_,
                     flag );

                //      E_[i][j]-=w_extend_;

                // update F_
                max3(F_[i][j],TF_[i][j],
                     G_[i-1][j]-w_open_,
                     E_[i-1][j]-w_open_,
                     F_[i-1][j]-w_extend_,
                     flag );

                //      F_[i][j]-=w_extend_;

#ifdef DEBUG
                cout << i << " " << j << " => ";
                cout << G_[i][j] << " " << E_[i][j] << " "
                     << F_[i][j] << " " << TG_[i][j] << " "
                     << TE_[i][j] << " " << TF_[i][j] << endl;
#endif
            } // ~for j
#ifdef DEBUG
            cout << endl;
#endif
        } // ~for i

        max=-10000;
        for (int i(0);i<=x_size;i++)
        {
            max3( thisMax, thisMatrix,
                  G_[i][y_size], E_[i][y_size], F_[i][y_size],flag );
            if (thisMax>max)
            {
                max=thisMax;
                ii=i;
                jj=y_size;
                whichMatrix=thisMatrix;
            } // ~if

        } // ~for i

        for (int j(0);j<=y_size;j++)
        {
            max3( thisMax, thisMatrix,
                  G_[x_size][j], E_[x_size][j], F_[x_size][j],flag );
            if (thisMax>max)
            {
                max=thisMax;
                ii=x_size;
                jj=j;
                whichMatrix=thisMatrix;
            } // ~if

        }

        score_=max;
#ifdef DEBUG
        cout << "start cell = " << ii << " " << jj << " " << whichMatrix << endl;
#endif
        x_start_=ii;
        y_start_=jj;

        // ================== filling up the end gaps ==================
        // ================== start ==================
        int x_fillEnd = x_size;
        int y_fillEnd = y_size;
        while( y_fillEnd > y_start_ ) {
            xt_+='-';
            yt_+=y[y_fillEnd-1];
            at_+=' ';
            y_fillEnd--;
        }
        while( x_fillEnd > x_start_ ) {
            xt_+=x[x_fillEnd-1];
            yt_+='-';
            at_+=' ';
            x_fillEnd--;
        }
        // ================== end ==================


        // check if this is a repeat candidate!

#ifdef ACTIVATE_REPEAT_DETECTION
        // in the case of negative scores, we have to make sure that we take the absolute value!
        ScoreType repeatThreshold = (max-( abs(max)*MIN_REPEAT_THRESHOLD));

#ifdef DEBUG
        cerr << "settting the repeat threshold = " << repeatThreshold << endl;
#endif

        bool repeat_found = false;
        for (int i(x_size);i<=x_size;i++)
        {
            int startIdx = ( (i-width_)< 1 )? 1:(i-width_);
            int endIdx = ( (i+width_)> y_size )? y_size:(i+width_);
          
            for (int j(startIdx);j<=endIdx;j++)
            {
                if( ( ( repeatThreshold < G_[i][j] ) ||
                      ( repeatThreshold < E_[i][j] ) ||
                      ( repeatThreshold < F_[i][j] )
                        )
                    && 
                    ( abs(j-jj)>(x_size/2) ) 
                    )
                {
                    if( repeat_found == false )
                    {
                        repeat_found = true;

#ifdef DEBUG
                        cerr << "G comparison    = " << ( repeatThreshold < G_[i][j] ) << "\t" << repeatThreshold << "\t" << G_[i][j] << endl
                             << "E comparison    = " << ( repeatThreshold < E_[i][j] ) << endl
                             << "F comparison    = " << ( repeatThreshold < F_[i][j] ) << endl
                             << "coordinate diff = " << ( abs(j-jj)>(x_size/2) ) << endl;
                      
                        cerr << "POTENTIAL REPEAT CANDIDATE" << endl;
                        cerr << "ii/jj/max = " << ii << "\t" << jj << "\t" << max << endl
                             << "i/j/G/E/F = " << i << "\t" << j << "\t" << G_[i][j] << "\t" << E_[i][j] << "\t" << F_[i][j] << endl;
#endif
                        break;
                    }
                }
            }
        }

#endif // END_OF_REPEAT_DETECTION


        DPMatrix* pT_ =&TE_;


        while ((ii>0)&&(jj>0))
        {
            if (whichMatrix==0)
            {
                pT_=&TG_;
            }
            else if (whichMatrix==1)
            {
                pT_=&TE_;
            }
            else
            {
                assert(whichMatrix==2);
                pT_=&TF_;
            }
            //  cout << TG_[6][7] << " " << TE_[6][7] << " " << TF_[6][7] << endl;

#ifdef DEBUG
            cout << ii << " " << jj << " " << whichMatrix << " ";
            cout << " "<< G_[ii][jj]<< ":"<< E_[ii][jj] << ":" << F_[ii][jj];
            cout << "/";
            cout << " "<< TG_[ii][jj]<< ":"<< TE_[ii][jj] << ":" << TF_[ii][jj];
            cout << " - " << (*pT_)[ii][jj] << endl;
#endif
            nextMatrix=(*pT_)[ii][jj];
            if (whichMatrix==0)
            {
                xt_+=x[ii-1];
                yt_+=y[jj-1];
                at_+=((x[ii-1]==y[jj-1])?'|':'.');
                ii--; jj--;
                // whichMatrix=0;
            } // ~if
            else if (whichMatrix==1)
            {
                xt_+='-';
                yt_+=y[jj-1];
                at_+=' ';
                jj--;
                //    whichMatrix=1;
            } // ~else if
            else
            {
                assert(whichMatrix==2);
                xt_+=x[ii-1];
                yt_+='-';
                at_+=' ';
                ii--;
                //      whichMatrix=2;
            } // ~else
            whichMatrix=nextMatrix;
        } // ~while
#ifdef DEBUG
        cout << "end cell = " << ii << " " << jj << endl;
#endif

        while( jj > 0 ) {
            xt_+='-';
            yt_+=y[jj-1];
            at_+=' ';
            jj--;
        }
        while( ii > 0 ) {
            xt_+=x[ii-1];
            yt_+='-';
            at_+=' ';
            ii--;
        } // ~else

#ifdef DEBUG
        cout << "end cell = " << ii << " " << jj << endl;
#endif





        x_end_=ii;
        y_end_=jj;

        reverse(xt_.begin(),xt_.end());
        reverse(at_.begin(),at_.end());
        reverse(yt_.begin(),yt_.end());
        //cout << xt_ <<endl;
        // cout << at_ <<endl;
        //  cout << yt_ <<endl;


        // everything went well
        return 1;

    } // ~ScoreType Aligner::operator()( const string& x, const string& y)





    // readAlignScoreFile: fn specifies the filename of the scoring file
    // the scoring files essentially contains 4 values for match/mismatch
    // and gap open and gap extend, lines beginning with "#" are ignored
    bool Aligner::readAlignScoreFile( ifstream& in, ScoreType& m, ScoreType& mm, ScoreType& go, ScoreType& ge )
    {
        m = 0;
        mm = 0;
        go = 0;
        ge = 0;

        string curLine;
        int nonCommentLinesRead = 0;
        while( !in.eof() ) {
            getline( in,curLine );
            if( curLine != "" ) {
                if( curLine[0] != '#' ) {

                    switch( nonCommentLinesRead ) {
                    case 0:
                        m = atoi( curLine.c_str() );
                        break;
                    case 1:
                        mm = atoi( curLine.c_str() );
                        break;
                    case 2:
                        go = atoi( curLine.c_str() );
                        break;
                    case 3:
                        ge = atoi( curLine.c_str() );
                        break;
                    default:
                        cerr << "too many lines to read in the alignment scoring file\n" << endl;
                        return false;
                    }
                    nonCommentLinesRead++;
                }
            }
        } // while

        if( m == 0 && mm == 0 && go == 0 && ge == 0 )
            return false;
        else
        {
            w_match_ = m;
            w_mismatch_ = mm;
            w_open_ = go;
            w_extend_ = ge;
        }

        return true;
    }

    string Aligner::convertToCIGAR( const string& alignedA,const string& alignedB )
    {
        if( alignedA.size() != alignedB.size() )
        {
            return "";
        }

        int match = 0;
        ostringstream out;

        for( uint i=0;i<alignedA.size();i++ ) {
            if( alignedA[i]==alignedB[i] ) {
                match++;
            } else if( alignedA[i] == '-' ) {
                // deletion
                if( match != 0 )
                {
                    out << match;
                }
                match = 0;

                out << "d";
            } else if( alignedB[i] == '-' ) {
                // insertion
                if( match != 0 )
                {
                    out << match;
                }
                match = 0;

                out << char(tolower(alignedA[i]));
            } else {
                // mismatch
                if( match != 0 )
                {
                    out << match;
                }
                match = 0;

                out << alignedB[i];
            }

        }

        if( match != 0 )
        {
            out << match;
        }

        return out.str();

    }


    // method for computing the standard match descriptor
    string Aligner::convertToAlignmentDescriptor( const string& alignedA,const string& alignedB,int& no_mismatches )
    {
        if( alignedA.size() != alignedB.size() )
        {
            return "";
        }

        no_mismatches = 0;

        int match = 0;
        ostringstream out;

        for( uint i=0;i<alignedA.size();i++ ) {
            if( alignedA[i]==alignedB[i] ) {
                match++;
            } else if( alignedA[i] == '-' || alignedB[i] == '-' ) {
                // we do not allow indels in the standard alignment descriptor
                return "";
            } else {
                // mismatch
                if( match != 0 )
                {
                    out << match;
                }
                match = 0;

                out << alignedB[i];
                no_mismatches++;
            }

        }

        if( match != 0 )
        {
            out << match;
        }

        return out.str();

    }


    string Aligner::convertToNewAlignmentDescriptor( const string& alignedA,const string& alignedB,int& no_mismatches,int& begin_offset,int& end_offset )
    {
        if( alignedA.size() != alignedB.size() )
        {
            return "";
        }

        if( (alignedA[0]!='-') && (alignedB[0]=='-') )
        {
            // the alignment does not seem right, something is going on
            return "";
        }

        if( (alignedA[alignedA.size()-1]!='-') && (alignedB[alignedB.size()-1]=='-') )
        {
            // the alignment does not seem right, something is going on
            return "";
        }


        no_mismatches = 0;

        // initializing
        begin_offset = 0;
        end_offset = 0;

        // determine the beginning of the read
        int start_idx = 0;
        int end_idx = alignedA.size()-1;

        while( (alignedA[start_idx]=='-') && (alignedB[start_idx]!='-') ) { start_idx++; }
        while( (alignedA[end_idx]=='-') && (alignedB[end_idx]!='-') ) { end_idx--; end_offset++; }

        // take indels in the first seed into account
        begin_offset = start_idx;


        int match = 0;
        int deleted_bases = 0;
        bool escape_mode = false;
        ostringstream out;

        for( uint i=start_idx;static_cast<int>(i)<=end_idx;i++ ) {
            if( alignedA[i]==alignedB[i] ) {
                if( escape_mode == true ) {
                    // print a "$" at the end
                    if( deleted_bases > 0 ) {
                        out << deleted_bases;
                        deleted_bases = 0;
                    }
                    out << "$";
                    escape_mode = false;
                }
                match++;

            } else if( alignedB[i] == '-' ) {
                // deletion
                if( match != 0 )
                {
                    out << match;
                }
                match = 0;

                if( escape_mode == false ) {
                    out << "^";
                    escape_mode = true;
                }

                //	out << "d";
                deleted_bases++;
            } else if( alignedA[i] == '-' ) {
                // insertion
                if( match != 0 )
                {
                    out << match;
                }
                match = 0;

                if( escape_mode == false ) {
                    out << "^";
                    escape_mode = true;
                }
                if( deleted_bases > 0 ) {
                    out << deleted_bases;
                    deleted_bases = 0;
                }

                //	out << char(tolower(alignedA[i]));
                out << alignedB[i];
            } else {
                // mismatch
                if( escape_mode == true ) {
                    // print a "$" at the end
                    if( deleted_bases > 0 ) {
                        out << deleted_bases;
                        deleted_bases = 0;
                    }
                    out << "$";
                    escape_mode = false;
                }

                if( match != 0 )
                {
                    out << match;
                }
                match = 0;

                out << alignedB[i];
                no_mismatches++;
            }

        }

        if( match != 0 )
        {
            out << match;
        }
        if( deleted_bases > 0 ) {
            out << deleted_bases;
        }
        if( escape_mode == true ) {
            out << "$";
        }

        return out.str();

    }


    // put your sanity checks here
    bool Aligner::checkAlignmentSanity( const string& read, const string& ref )
    {
        if( (read.size() != ref.size()) || (read.size()==0) || (ref.size()==0) )
            return false;

        // we expect some gaps in front of the read (ALIGN_DP_BAND/2 bases
        // are pulled out to the left and right of the read)
        if( read[0]!='-' && ref[0]=='-')
            return false;


        return true;
    }



} } // namespace ga
