// Copyright (c) 2008 Illumina
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).

// PROJECT: ELAND (Efficient Large-Scale Alignment of Nucleotide Databases)
// MODULE: alignExtendedFrags.cpp
// AUTHOR: M. Bauer

// Description: perform a gapped alignment with reads that have not been mapped onto the genome


#ifndef GA_ALIGNER_H
#define GA_AlIGNER_H

using namespace std;

typedef unsigned char uchar;

#define QUALSCALE 0.0333
#define EPSILON 0.0000001
#define UINT_INIT 1048576 // 2^20

namespace ga { namespace alignment {

    typedef double ScoreType;

    struct DPMatrix : public   vector<vector<ScoreType> >
    {
        DPMatrix( void ) {}
        DPMatrix( int x, int y)
        {
            resize(x,y);
        } // ~ctor

        void resize( int x, int y)
        {
            vector<vector<ScoreType> >::resize(x);
            for (int j(0);j<x;j++)
                (*this)[j].resize(y);
        } // ~resize;

    };

    struct Aligner
    {
        public:
Aligner( ScoreType w_match,
         ScoreType w_mismatch,
         ScoreType w_open,
         ScoreType w_extend,
         int w,
         int insertsize,
         int stdDeviation ) :
    w_match_(w_match),
            w_mismatch_(w_mismatch),
            w_open_(w_open),
            w_extend_(w_extend),
            width_(w),
            init_(false),
            expected_insertsize_(insertsize),
            expected_stdDeviation_(stdDeviation)
        {
            qualScaling_ = QUALSCALE;
            allowInserts_ = 1;
            allowDeletions_ = 1;
        }

        int operator()( const char* qvals_,const char* x, const char* y,
                              const int x_size, const int y_size,const bool& flag );

        int operator()( const char* qvals_,const string& x, const string& y,const bool& flag)
        {
            return (*this)(qvals_,x.c_str(),y.c_str(),x.size(),y.size(),flag);
        } // ~operator()


        bool readAlignScoreFile( ifstream& in, ScoreType& m, ScoreType& mm, ScoreType& go, ScoreType& ge );


        //protected:
        ScoreType w_match_;
        ScoreType w_mismatch_;
        ScoreType w_open_;
        ScoreType w_extend_;
        int width_;
        double qualScaling_; // scaling factor for the quality string values
        short allowInserts_;
        short allowDeletions_;
        bool init_;

        const int expected_insertsize_;
        const int expected_stdDeviation_;
        int expectedGapHit_;


        DPMatrix E_,F_,G_,TE_,TF_,TG_;

        string xt_;
        string yt_;
        string at_;
        int x_start_;
        int y_start_;
        int x_end_;
        int y_end_;

        ScoreType score_;


        void init( const int x_size, const int y_size, const int expected_gap, const int )
        {
            E_ = DPMatrix(x_size+1,y_size+1);
            F_ = DPMatrix(x_size+1,y_size+1);
            G_ = DPMatrix(x_size+1,y_size+1);
            TE_ = DPMatrix(x_size+1,y_size+1);
            TF_ = DPMatrix(x_size+1,y_size+1);
            TG_ = DPMatrix(x_size+1,y_size+1);

            expectedGapHit_ = expected_gap;

            init_ = true;
        } // init

        void max3( ScoreType& max, ScoreType& which,
                   const ScoreType& v0,
                   const ScoreType& v1,
                   const ScoreType& v2,
                   const bool& forward)
        {
            if( forward == true )
            {
                max=v0; which=0;
                if ( (v1-v0) > EPSILON )  { max=v1; which=1; }
                if ( (v2-max) > EPSILON ) { max=v2; which=2; }
            }
            else
            {
                max=v2; which=2;
                if ( (v1-v2) > EPSILON )  { max=v1; which=1; }
                if ( (v0-max) > EPSILON ) { max=v0; which=0; }
            }

        } // ~Aligner::max3



        // GET/SET METHODS -------------------------------------------------------------------------
        void setQualityScaling( const double& scaling ) { qualScaling_ = scaling; }
        double getQualityScaling( void ) { return qualScaling_; }

/*     void setWidth( const int& val ) { width_ = val; } */
/*     int getWidth( void ) { return width_; } */

        void allowInserts( const bool& val ) { allowInserts_ = val; }
        bool insertsAllowed( void ) { return( allowInserts_ == 1); }

        void allowDeletions( const bool& val ) { allowDeletions_ = val; }
        bool deletionsAllowed( void ) { return( allowDeletions_ == 1); }

        string convertToCIGAR( const string& alignedA,const string& alignedB );
        string convertToAlignmentDescriptor( const string& alignedA,const string& alignedB,int& no_mismatches );
        string convertToNewAlignmentDescriptor( const string& alignedA,const string& alignedB,int& no_mismatches,int& begin_offset,int& end_offset );

        bool checkAlignmentSanity( const string& read, const string& ref );


    }; // ~class Aligner


} } // namespace ga


#endif
