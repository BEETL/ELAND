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
 ** \file orphanAligner.cpp
 **
 ** \brief Phage align.
 **
 ** Phage align.
 **
 ** \author Markus Bauer
 **/

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <locale>

#include <boost/shared_ptr.hpp>

#include "alignment/ELAND_unsquash.h"
#include "alignment/aligner.h"
#include "common/ExtendedFileReader.h"


#define FRAGMENT_LENGTH 450
#define UPPER_BOUND_OCC 1
#define MAX_ERROR_RATE 0.1

#define REQUEST_SIZE 262144


// singleton request
struct SingletonRequest
{
    SingletonRequest(
        int readNum,
        int left_or_right,
        uint fileIndex,
        uint contigNum,
        uint filePos,
        char strand,
        string orphan ) :
        readNum_(readNum),
        left_or_right_(left_or_right),
        fileIndex_(fileIndex),
        contigNum_(contigNum),
        filePos_(filePos),
        strand_(strand),
        orphan_(orphan)
    {}
    int readNum_;
    int left_or_right_;
    uint fileIndex_;
    uint contigNum_;
    uint filePos_;
    char strand_;
    string orphan_;
}; // ~struct SeqRequest

// singleton request
struct SingletonAlignment
{
    SingletonAlignment(
        int readNum,
        int left_or_right,
        uint fileIndex,
        uint contigNum,
        int filePos,
        char strand,
        string orphan,
        int aligned_position,
        string matchDesc ) :
        readNum_(readNum),
        left_or_right_(left_or_right),
        fileIndex_(fileIndex),
        contigNum_(contigNum),
        filePos_(filePos),
        strand_(strand),
        orphan_(orphan),
        aligned_position_(aligned_position),
        matchDesc_(matchDesc)
    {}
    int readNum_;
    int left_or_right_;
    uint fileIndex_;
    uint contigNum_;
    int filePos_;
    char strand_;
    string orphan_;
    int aligned_position_;
    string matchDesc_;
}; // ~struct SeqRequest


bool parseMatch( const string& match,string &chr,char& strand,uint& pos );
void parseXYZ( const char* xyz,int& oneError,int& twoError );
string reverseComplement( const string& s );

int global_upper_bound_occ = UPPER_BOUND_OCC;
int global_fragment_length = FRAGMENT_LENGTH;



class OrphanAligner {
public:
    OrphanAligner( const string& squashed_genome,
                   const ga::alignment::ScoreType match,
                   const ga::alignment::ScoreType mismatch,
                   const ga::alignment::ScoreType gapopen,
                   const ga::alignment::ScoreType gapextend,
                   const int& read_length,
                   const int& maxNumberMismatches,
                   const int& expected_insertsize,
                   const int& expDeviation
        ) : squashed_genome_(squashed_genome),
            align_forward_(match,mismatch,gapopen,gapextend,global_fragment_length,expected_insertsize,expDeviation),
            align_reverse_(match,mismatch,gapopen,gapextend,global_fragment_length,expected_insertsize,expDeviation),
            maxNumberMismatches_(maxNumberMismatches)
    {
        // if the read length was greater then the expected insert size, we'd have overlapping reads
//        assert(expected_insertsize>read_length);

        // calculate where to jump, F -> jump over the read + 50 bases more
        jump_forward_ = read_length; // for the moment, do not jump forward any additional bases

        // TODO: estimate the actual fragment size from the expected insert size and the read lengths;
        // right now it's only hardcoded-values


        // 300bp quality string --> string of Q30 qualities scales the match/mismatch contributions to 1
        qual_string_ = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";

        int expected_gap_forward = -1;
        int expected_gap_reverse = -1;
        // set the gaps favourably for the forward alignment
        expected_gap_forward = global_fragment_length - expected_insertsize;

        // set the gaps favourably for the reverse alignment
        expected_gap_reverse = expected_insertsize - ( jump_forward_ );

        cerr << "expeced_gap_forward/expected_gap_reverse : "
             << expected_gap_forward << "\t" << expected_gap_reverse << endl;


        align_forward_.init( read_length,global_fragment_length,expected_gap_forward,expDeviation );
        align_forward_.allowInserts( 1 );
        align_forward_.allowDeletions( 1 );

        align_reverse_.init( read_length,global_fragment_length,expected_gap_reverse,expDeviation );
        align_reverse_.allowInserts( 1 );
        align_reverse_.allowDeletions( 1 );

#ifdef DEBUG
        cerr << "setting jump_forward = " << jump_forward_ << endl;
#endif

    }


    ~OrphanAligner() {};

    vector<SingletonAlignment> pullOutFragments( StringIndex& files,vector<SingletonRequest>& req,uint& orphans_rescued );


private:
    boost::shared_ptr<SquashFile> pSquash_;
    string squashed_genome_;
    string qual_string_;
    ga::alignment::Aligner align_forward_;
    ga::alignment::Aligner align_reverse_;
    uint jump_forward_;
    uint jump_backward_;
    // threshold for orphans to be rescued. or not.
    const int maxNumberMismatches_;

};





bool lessThanSingleton(const SingletonRequest& a, const SingletonRequest& b)
{
  //  return ((a.fileIndex_<b.fileIndex_)
  //	  ||((a.fileIndex_==b.fileIndex_)&&(a.filePos_<b.filePos_)));
  return ((a.fileIndex_<b.fileIndex_)
	  ||((a.fileIndex_==b.fileIndex_)&&(a.contigNum_<b.contigNum_))
	  ||((a.fileIndex_==b.fileIndex_)
	     &&(a.contigNum_==b.contigNum_)
	     &&(a.filePos_<b.filePos_)));
} // ~bool lessThanRequest;

bool lessThanSingletonAlignment(const SingletonAlignment& a, const SingletonAlignment& b)
{
    return ( (a.readNum_<b.readNum_)
             ||( (a.readNum_==b.readNum_)&&(a.fileIndex_<b.fileIndex_) )
        );
}



int main( int numArgs, const char** args)
{
    if (numArgs != 6 )
    {
        cerr << "usage : " << args[0] << " <extended1> <extended2> <genome> <output_suffix> <upper_bound>" << endl
             << "(input file has to be in the export file format, genome has to point to a squashed genome directory)" << endl;
        exit (1);
    }


    string squashed_genome = string(args[3]);
    StringIndex files(squashed_genome);
    int cur_request = 0;
    vector<SingletonRequest> req;


    ExtendedFileReaderActual i_export_left(args[1]);
    ExtendedFileReader left(&i_export_left);
    ExtendedFileReaderActual i_export_right(args[2]);
    ExtendedFileReader right(&i_export_right);

    ofstream out_left_extended( string(string(args[1])+ string(args[4]) ).c_str() );
    ofstream out_right_extended( string(string(args[2])+ string(args[4]) ).c_str() );

    // check that the eland_extended files are not empty
    if (!left.getNextEntry())
    {
        cerr << "WARNING: input file " << args[1] << " is empty" << endl;
        return 0;
    }
    left.rewind();

    cerr << "estimating insert size...";
    vector<int> insertSizes;
    // come up with statistics for the insert size distribution
    while( left.getNextEntry() == true )
    {
        right.getNextEntry();

        int left_matchcounter   = left.getMatchCounter();
        string left_matches   = string(left.getMatches());

        int right_matchcounter   = right.getMatchCounter();
        string right_matches   = string(right.getMatches());

        int currentInsertSize = -1;

        // we are only looking at read pairs where either read can be placed uniquely
        if( (left_matchcounter == 1) && (right_matchcounter==1) )
        {
            string left_match_chr = "";
            uint left_match_position = 0;
            char left_strand = 'F';
            string right_match_chr = "";
            uint right_match_position = 0;
            char right_strand = 'F';

            parseMatch( left_matches,left_match_chr,left_strand,left_match_position );
            parseMatch( right_matches,right_match_chr,right_strand,right_match_position );

            // we have to be on the very same chromosome
            if( left_match_chr != right_match_chr )
            {
                continue;
            }

            // we are only looking at read pairs that have proper orientation
            if( left_strand=='F' && right_strand=='R' )
            {
                // left_match_position < right_match_position
                if( left_match_position < right_match_position )
                {
                    currentInsertSize = right_match_position - left_match_position;
                }
            }
            else if( left_strand=='R' && right_strand=='F' )
            {
                // right_match_position < left_match_position
                if( right_match_position < left_match_position )
                {
                    currentInsertSize = left_match_position - right_match_position;
                }
            }

        }



        if( currentInsertSize != -1 )
        {
            insertSizes.push_back( currentInsertSize );
//            cout << currentInsertSize << endl;
        }
    }
    cerr << "done." << endl;

    // insertSizes now hold all the information
    // get the median
    sort( insertSizes.begin(),insertSizes.end() );
    int medianDistribution = 0;
    // calculate variance + standard deviation
    double varianceDistribution = 0.0;
    double d_stdDeviationDistribution = 0.0;
    int i_stdDeviationDistribution = 0;
    int meanDistribution = 0;


    if( insertSizes.size() > 0 )
    {
        medianDistribution = insertSizes[ insertSizes.size()/2 ];

        // get the mean
        int sumOfAll = 0;
        for( unsigned int i=0;i<insertSizes.size();i++ )
        {
            sumOfAll += insertSizes[i];
        }
        meanDistribution = sumOfAll/insertSizes.size();

        cerr << "median of distribution : " << medianDistribution << endl;

        // calculate variance + standard deviation
        for( unsigned int i=0;i<insertSizes.size();i++ )
        {
            varianceDistribution += (double)((insertSizes[i]-meanDistribution)*(insertSizes[i]-meanDistribution))/(double)(insertSizes.size()-1);
        }
        d_stdDeviationDistribution = sqrt(varianceDistribution);
        i_stdDeviationDistribution = static_cast<int>(d_stdDeviationDistribution);

        cerr << "variance of distribution : " << varianceDistribution << endl;
        cerr << "standard deviation of distribution : " << d_stdDeviationDistribution << endl;
    }

    // rewind everything
    left.rewind();
    right.rewind();

    // setting the number of occurrences, upper bound
    global_upper_bound_occ=atoi(args[5]);


    ga::alignment::ScoreType match = 6;
    ga::alignment::ScoreType mismatch = -1;
    ga::alignment::ScoreType gapopen = 15;
    ga::alignment::ScoreType gapextend = 3;

    // determine read length
    left.getNextEntry();
    int left_read_length = string(left.getRead()).size();
    left.rewind();
    right.getNextEntry();
    int right_read_length = string(right.getRead()).size();
    right.rewind();

    int read_length = (left_read_length>right_read_length)?left_read_length:right_read_length;



    // compute the number of mismatches that we allow for an orphan to be rescued
    const int maxNumberMismatches = (int)(read_length*MAX_ERROR_RATE);
    cerr << "Setting the orphan rescue threshold to " << maxNumberMismatches << endl;


    OrphanAligner oa( squashed_genome,
                      match,
                      mismatch,
                      gapopen,
                      gapextend,
                      read_length,
                      maxNumberMismatches,
                      medianDistribution,
                      i_stdDeviationDistribution );
    vector<SingletonAlignment> alns;

    uint cur_line_idx = 0;
    uint orphans_rescued = 0;
    uint total_no_of_candidates = 0;


    // main loop
    while( left.getNextEntry()==true ) {
        right.getNextEntry();

        int left_matchcounter   = left.getMatchCounter();
        string left_matches   = string(left.getMatches());

        int right_matchcounter   = right.getMatchCounter();
        string right_matches   = string(right.getMatches());

        string match_chr = "";
        uint i_match_position = 0;
        char strand = 'F';
        string orphan_read = string(left.getRead());

        short candidate = 0;


        /// split up the input at the ',' marks
        vector<string> v_chrom;
        vector<char> v_strand;
        vector<uint> v_match_position;
        vector<string> split_matches;


        // matchcounter returns the minimal number of the X:Y:Z read:
        // this is 0 for NM and QC
        // if min_numer > 10 then we do not report any hits
        if( (left_matchcounter > 0) && (left_matchcounter <= global_upper_bound_occ) && right_matchcounter == 255 )
        {
            size_t pos_in_string = 0;
            size_t next_comma = -1;

#ifdef DEBUG
            cerr << "left_matches: " << left_matches << "\t" << left.getXYZ() << "\t" << right.getXYZ() << endl;
#endif

            next_comma = left_matches.find(",",pos_in_string);
            while( next_comma != string::npos )
            {
                string cur_substr = left_matches.substr( pos_in_string,(next_comma-pos_in_string) );

#ifdef DEBUG
                cerr << "SPLIT : " << cur_substr << endl;
#endif

                split_matches.push_back( cur_substr );

                pos_in_string = next_comma + 1;
                next_comma = left_matches.find(",",pos_in_string);
            }
#ifdef DEBUG
            cerr << "SPLIT : " << left_matches.substr( pos_in_string,(next_comma-pos_in_string) ) << endl;
#endif

            split_matches.push_back( left_matches.substr( pos_in_string,1024 )  );


            for( vector<string>::iterator it=split_matches.begin();it!=split_matches.end();it++ )
            {
                // take left.getMatches() and parse the hit
                // chrX.fa:55351604F102C5
                parseMatch( *it,match_chr,strand,i_match_position );
                orphan_read = string(right.getRead());
                candidate = 2;

#ifdef DEBUG
                cerr << match_chr << "/" << strand << "/" << i_match_position << endl;
#endif

                // push them back
                if( i_match_position > 0 )
                {
                    v_chrom.push_back( match_chr );
                    v_strand.push_back( strand );
                    v_match_position.push_back( i_match_position );
                }
            }

#ifdef DEBUG
            cerr << endl << endl;
#endif


        }

        if( left_matchcounter == 255 && (right_matchcounter > 0) && (right_matchcounter <= global_upper_bound_occ) )
        {
            size_t pos_in_string = 0;
            size_t next_comma = -1;

#ifdef DEBUG
            cerr << "right_matches: " << right_matches << "\t" << left.getXYZ() << "\t" << right.getXYZ() << endl;
#endif

            next_comma = right_matches.find(",",pos_in_string);
            while( next_comma != string::npos )
            {
                string cur_substr = right_matches.substr( pos_in_string,(next_comma-pos_in_string) );

#ifdef DEBUG
                cerr << "SPLIT : " << cur_substr << endl;
#endif

                split_matches.push_back( cur_substr );

                pos_in_string = next_comma + 1;
                next_comma = right_matches.find(",",pos_in_string);
            }

#ifdef DEBUG
            cerr << "SPLIT : " << right_matches.substr( pos_in_string,(next_comma-pos_in_string) ) << endl;
#endif
            split_matches.push_back( right_matches.substr( pos_in_string,1024 )  );


            for( vector<string>::iterator it=split_matches.begin();it!=split_matches.end();it++ )
            {
                // take right.getMatches() and parse the hit
                // chrX.fa:55351604F102C5
                parseMatch( *it,match_chr,strand,i_match_position );
                orphan_read = string(left.getRead());
                candidate = 1;

#ifdef DEBUG
                cerr << match_chr << "/" << strand << "/" << i_match_position << endl;
#endif

                if( i_match_position > 0 )
                {
                    v_chrom.push_back( match_chr );
                    v_strand.push_back( strand );
                    v_match_position.push_back( i_match_position );
                }

            }

#ifdef DEBUG
            cerr << endl << endl;
#endif

        }

        // we need the following to make sure that the assert statement does not kick off
        // when we have a hit at position 0 or < 0
        if( v_match_position.size() > 0 )
        {
            for( unsigned int ii=0;ii<v_match_position.size();ii++ )
            {
                i_match_position = v_match_position[ii];
                match_chr = v_chrom[ii];
                strand = v_strand[ii];

                // sanity checks
                assert( i_match_position != 0 );
                assert( match_chr != "" );

                uint chromNum = 0;
                uint contigNum = 0;


                files.getIndex( match_chr.c_str(),
                                chromNum,
                                contigNum,
                                i_match_position
                    );

//                cerr << "putting together: " << match_chr << "\t" << chromNum << "\t" << contigNum << "\t" << i_match_position << "\t" << strand << "\t" << orphan_read << endl;

                SingletonRequest singleton(
                    cur_line_idx,
                    candidate,
                    chromNum,
                    contigNum,
                    i_match_position,
                    strand,
                    orphan_read );

                req.push_back(singleton);
                cur_request++;


            }


            if( cur_request > REQUEST_SIZE )
            {
                cerr << "pulling out fragments..." << endl;
                total_no_of_candidates += cur_request;
                vector<SingletonAlignment> orphan_alignments  = oa.pullOutFragments( files,req,orphans_rescued );
                alns.insert( alns.end(),orphan_alignments.begin(),orphan_alignments.end() );


                {
                    vector<SingletonRequest> req_tmp;
                    req.swap( req_tmp );
                    req.clear();
                }
                cur_request = 0;


            } // else

        } // candidate

        cur_line_idx++;
    } // while


    if( req.size() > 0 )
    {
        vector<SingletonAlignment> orphan_alignments  = oa.pullOutFragments( files,req,orphans_rescued );
        alns.insert( alns.end(),orphan_alignments.begin(),orphan_alignments.end() );

        total_no_of_candidates += cur_request;


        req.clear();
        cur_request = 0;
    } // if( req.size() ... )

    cerr << "total number of candidates = " << total_no_of_candidates << endl;
    cerr << "orphans rescued            = " << orphans_rescued << endl;
    cerr << "list size                  = " << alns.size() << endl;

    cerr << "sorting singleton alignments...";
    sort( alns.begin(),alns.end(),lessThanSingletonAlignment );
    cerr << "done." << endl;


    // =================================== OUTPUT ===================================
    cerr << "writing output...";
    // print the new alignments to the extended files
    cerr << "we got " << alns.size() << " orphans to write." << endl;

    left.rewind();
    right.rewind();
    int line_cnt_extended = 0;
    uint cur_idx_alns = 0;
    while( left.getNextEntry()==true ) {
        right.getNextEntry();

        if( alns.size() > 0 )
        {
            if( alns[cur_idx_alns].readNum_ == line_cnt_extended )
            {
                short left_or_right = alns[cur_idx_alns].left_or_right_;
                int cnt_matches = 0;

                // build up the new match string, extend it to a contig level,
                // we only have to extend the loop below
                ostringstream new_match;
                uint cur_chrom = UINT_INIT; // initializes to 2^XX
                uint cur_contig = UINT_INIT; // initializes to 2^XX
                int cur_offset = 0;

                while( (alns[cur_idx_alns].readNum_ == line_cnt_extended) && (cur_idx_alns<alns.size()) )
                {
                    if( cur_chrom != UINT_INIT )
                    {
                        new_match << ",";
                    }

                    bool changed_contig = false;
                    if( alns[cur_idx_alns].fileIndex_ != cur_chrom )
                    {
                        cur_chrom = alns[cur_idx_alns].fileIndex_;
                        changed_contig = true;
                    }
                    if( alns[cur_idx_alns].contigNum_ != cur_contig )
                    {
                        cur_contig = alns[cur_idx_alns].contigNum_;
                        changed_contig = true;
                    }

                    if( changed_contig == true )
                    {
                        string contigName = files.getContigName( cur_chrom,cur_contig,cur_offset );
                        new_match << files.names_[ cur_chrom ];
                        if( contigName != "" )
                        {
                            new_match << "/" << contigName;
                        }
                        new_match << ":";
                    }
                    new_match << (alns[cur_idx_alns].aligned_position_-cur_offset) << alns[cur_idx_alns].strand_ << alns[cur_idx_alns].matchDesc_;

                    // don't forget to increment
                    cnt_matches++;
                    cur_idx_alns++;
                }


                if( left_or_right == 1 )
                {
                    // parse left.getXYZ()
                    int oneE = 0;
                    int twoE = 0;
                    parseXYZ( left.getXYZ(),oneE,twoE );
                    ostringstream string_xyz;
//                    string_xyz << "1:" << oneE << ":" << twoE;
                    string_xyz << cnt_matches << ":" << oneE << ":" << twoE;


                    // write new hit in read1 file
                    out_left_extended << left.getMachine() << "\t" << left.getRead() << "\t" << string_xyz.str() << "\t" << new_match.str() << endl;
                    // write standard for read2
                    out_right_extended << right.getMachine() << "\t" << right.getRead() << "\t" << right.getXYZ() << "\t" << right.getMatches();
                }
                else
                {
                    // parse right.getXYZ()
                    int oneE = 0;
                    int twoE = 0;
                    parseXYZ( right.getXYZ(),oneE,twoE );
                    ostringstream string_xyz;
//                    string_xyz << "1:" << oneE << ":" << twoE;
                    string_xyz << cnt_matches << ":" << oneE << ":" << twoE;


                    // write standard hit in read1 file
                    out_left_extended << left.getMachine() << "\t" << left.getRead() << "\t" << left.getXYZ() << "\t" << left.getMatches();
                    // write new hit in read2 file
                    out_right_extended << right.getMachine() << "\t" << right.getRead() << "\t" << string_xyz.str() << "\t" << new_match.str() << endl;
                }

                // we already increment the number above
                // cur_idx_alns++;
            }
            else
            {
                out_left_extended << left.getMachine() << "\t" << left.getRead() << "\t" << left.getXYZ() << "\t" << left.getMatches();
                out_right_extended << right.getMachine() << "\t" << right.getRead() << "\t" << right.getXYZ() << "\t" << right.getMatches();
            }
        }
        else
        {
            out_left_extended << left.getMachine() << "\t" << left.getRead() << "\t" << left.getXYZ() << "\t" << left.getMatches();
            out_right_extended << right.getMachine() << "\t" << right.getRead() << "\t" << right.getXYZ() << "\t" << right.getMatches();
        }


        line_cnt_extended++;
    }


    cerr << "done." << endl << "mission accomplished." << endl;




} // main



bool parseMatch( const string& match,string &chr,char& strand,uint& pos )
{
    const string const_chr_del = ":";
    const string const_comma_del = ",";
    const string const_strand_del = "RF";

    size_t chr_del = -1;
    size_t strand_del = -1;
    //size_t comma_del = -1; // denotes the pos
    size_t pos_in_string = 0;
    size_t pos_after_chr_del = 0;

//    cerr << "parsing " << match << endl;

    // determine chromosome
    chr_del = match.find(const_chr_del,pos_in_string);
    if( chr_del != string::npos ) { pos_after_chr_del = chr_del+1; }
    strand_del = match.find_first_of(const_strand_del,pos_after_chr_del);
    // comma_del = match.find(const_comma_del,pos_in_string);

    if( strand_del == string::npos )
    {
        return false;
    }

    if( chr_del != string::npos )
    {
        // update the current chromosome, also update pos_in_string, because now we are parsing the string after ":"
        chr = match.substr(0,chr_del);
        pos_in_string = chr_del+1;
    }
    strand = match.substr(strand_del,1)[0];

    int pos_tmp = (int)atoi( match.substr(pos_in_string, (strand_del - pos_in_string) ).c_str() );
    if( pos_tmp >= 0 ) { pos = pos_tmp; }
    else { pos = 0; }



    return true;


}


vector<SingletonAlignment> OrphanAligner::pullOutFragments( StringIndex& files,vector<SingletonRequest>& req,uint& orphans_rescued )
{
    vector<SingletonAlignment> res;

    if( req.size() == 0 )
    {
        return res;
    }

    sort( req.begin(),req.end(),lessThanSingleton );

    uint cur_file_idx = UINT_INIT;

    for( uint j=0;j<req.size();j++ )
    {
#ifdef DEBUG
        cerr << "CUR_REQUEST:\t" << req[j].fileIndex_ << "\t" << req[j].filePos_ << "\t" << req[j].strand_ << endl;
#endif

        if( cur_file_idx != req[j].fileIndex_ )
        {
            pSquash_=boost::shared_ptr<SquashFile>(
                new SquashFile(squashed_genome_, files.names_[req[j].fileIndex_], files));
            cur_file_idx = req[j].fileIndex_;
        }

        uint adapted_pos = (req[j].filePos_-1);

        if( req[j].strand_=='R' )
        {

            // check that we are not at the beginning of a chromosome,
            // creating an integer underflow
	  //            if( ((int)adapted_pos - global_fragment_length)<0 ) // OLD CHECK DOES NOT WORK FOR unsigned ints!
	  if( adapted_pos < (uint)global_fragment_length )
            {
                continue;
            }

            adapted_pos -= global_fragment_length;

        }
        else
        {
            adapted_pos += jump_forward_;
        }

#ifdef DEBUG
        cerr << "adapted position = " << adapted_pos << endl;
#endif


        string result = "";
        string non_reversed_result = "";
        pSquash_->goToPos
            ( req[j].contigNum_,adapted_pos );

        if( req[j].strand_=='R' )
        {
            for (int jj1(0);jj1<global_fragment_length;jj1++)
            {
                result+=pSquash_->getNextBase();
            }
        }
        else
        {
            // pull out the reverse complement
            for (int jj1(0);jj1<global_fragment_length;jj1++)
            {
                char cur_base = pSquash_->getNextBase();
                result += reverseCharASCII[(uint)(cur_base) ];
                non_reversed_result+=cur_base;
            }
            reverse( result.begin(),result.end() );
        }

#ifdef DEBUG
        cerr << "FRAGMENT PULLED OUT: " << result << endl;
        cerr << req[j].strand_ << endl
             << req[j].orphan_ << endl
             << result << endl << endl;
#endif

        ga::alignment::Aligner* curAligner = &align_forward_;
        if( req[j].strand_ == 'F' )
        {
            // we need to take the align_forward_object
            curAligner = &align_reverse_;
        }

        int retCode = (*curAligner)( qual_string_.c_str(),req[j].orphan_.c_str(),result.c_str(),req[j].orphan_.size(),result.size(),(req[j].strand_=='R') );

        assert(retCode>0); // assume everything went fine

        // return code of 2 is a repeat candidate
        if( retCode == 2 )
        {
#ifdef DEBUG
            cerr << "------------------------------------------------------------" << endl;
#endif
        }


#ifdef DEBUG
        int i=0;
        int reverse_i = curAligner->xt_.size()-1;
        while( i<curAligner->xt_.size() && (curAligner->xt_[i]=='-') ) {i++;}
        while( reverse_i>0 && (curAligner->xt_[reverse_i]=='-') ) {reverse_i--;}
        assert(reverse_i>i);
        int no_matches = 0;
        int sizeOverlap = reverse_i-i;
        for( int overlap_i=i;overlap_i<=reverse_i;++overlap_i )
        {
            if( curAligner->xt_[overlap_i]==curAligner->yt_[overlap_i] ) { no_matches++; }
        }

        if( ((double)(no_matches)/(double)(sizeOverlap))>0.95 )
        {
            cerr << req[j].strand_ << "\t" << ((double)(no_matches)/(double)(sizeOverlap)) << endl
                 << curAligner->xt_ << endl
                 << curAligner->yt_ << endl << endl;
        }

#endif

#ifdef SANITY_CHECK
        if( req[j].strand_ == 'F' )
        {
            // check if we get the same alignment if we reverse-complement the read and align it against the non-reversed reference
            string orphan_revComplement = reverseComplement( req[j].orphan_ );

            align_( qual_string_.c_str(),orphan_revComplement.c_str(),non_reversed_result.c_str(),orphan_revComplement.size(),non_reversed_result.size() );

            cerr << endl << "NON REVERSED" << endl;
            cerr << curAligner->xt_ << endl
                 << curAligner->yt_ << endl;

        } else
        {
            cerr << endl;
        }
#endif

#ifdef CUT_OUT_THE_ALIGNED_FRAGMENT
        int i=0;
        while( i<curAligner->xt_.size() && (curAligner->xt_[i]=='-') ) {i++;}

        if( (i-20+160)<curAligner->xt_.size() && (i>20) )
        {
            cerr << curAligner->xt_.substr(i-20,160) << endl
                 << curAligner->yt_.substr(i-20,160) << endl;

            int cur_score = 0;
            while(i<curAligner->xt_.size() && curAligner->xt_[i]!='-' )
            {
                if( curAligner->xt_[i]==curAligner->yt_[i] ) { cur_score+=2; } else { cur_score-=1; }
                i++;
            }
        }
#endif
        int mismatches = 0;
        int offset_begin = 0;
        int offset_end = 0;
        string new_ad = curAligner->convertToNewAlignmentDescriptor(curAligner->xt_,curAligner->yt_,mismatches,offset_begin,offset_end );

#ifdef DEBUG
        cerr << "DEBUG:" << endl
             << curAligner->xt_ << endl
             << curAligner->yt_ << endl
             << "# mismatches/ad: " << mismatches << "\t" << new_ad << endl;
#endif


        if( new_ad != "" && (mismatches<maxNumberMismatches_) )
        {
            // compute the match position
            int orphan_match_position = adapted_pos + 1; // offset_begin;
            if( req[j].strand_ == 'F' )
            {
                orphan_match_position += offset_end;
            }
            else
            {
                orphan_match_position += offset_begin;
            }

            // pull together SingletonAlignment
            res.push_back( SingletonAlignment(
                               req[j].readNum_,
                               req[j].left_or_right_,
                               req[j].fileIndex_,
                               req[j].contigNum_,
                               req[j].filePos_,
                               (req[j].strand_=='F')?'R':'F',
                               req[j].orphan_,
                               orphan_match_position,
                               new_ad )
                );
            orphans_rescued++;

#ifdef DEBUG
            cerr << "ORPHAN RESCUED" << endl;
#endif
        }
        else
        {
#ifdef DEBUG
            cerr << "ORPHAN NOT RESCUED" << endl;
#endif
        }

#ifdef DEBUG
        cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;
#endif

    } // for

    return res;
}


void parseXYZ( const char* xyz,int& oneError,int& twoError )
{
    // now we also have to consider the case in which the second reads
    // did not map multiple times, but was unmapped!
    if( xyz[0] == 'N' || xyz[0] == 'Q' || xyz[0] == 'R' ) //expected non-digit cases are NM, QC or RM
    {
        oneError = 0;
        twoError = 0;
    }
    else
    {
        xyz=strchr(xyz,':');
        BOOST_ASSERT(xyz && "Cannot parse one-error match count out of neighbourhood record");

        oneError = atoi(++xyz);

        xyz=strchr(xyz,':');
        BOOST_ASSERT(xyz && "Cannot parse two-error match count out of neighbourhood record");

        twoError = atoi(++xyz);
    }
}


string reverseComplement( const string& s )
{
  string result = "";

  for( int i=s.size()-1;i>=0;i-- ) {
    if( s[i] == 'A' ) {
      result += 'T';
    } else if( s[i] == 'T' ) {
      result += 'A';
    } else if( s[i] == 'G' ) {
      result += 'C';
    } else if( s[i] == 'C' ) {
      result += 'G';
    } else {
      result += s[i];
    }
  }

  return result;
}
