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
 ** \file eland_ms/MatchTable.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 ** Contains interface functionality from ELAND - specifically, stuff to do
 ** with scoring and storage of alignments
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_MATCH_TABLE_H
#define CASAVA_ELAND_MS_MATCH_TABLE_H

#include "ElandConstants.hh"
#include "MatchDescriptor.hh"
#include "MatchPositionTranslator.hh"
#include "SuffixScoreTable.hh"
#include "MultiMatch.hh"
#include "TableEntry.hh"
#include "ElandDefines.hh"
#include "common/StreamUtil.h"

namespace casava
{
namespace eland_ms
{


// stores arguments to MatchTable.addMatch(), a set of MatchStores is
// held by MatchCache
//
// position - number of the oligo found
// matchPos - where it was found in the genome
// numErrors - number of errors in the oligo suffix
//
struct MatchStore {
    inline
    void
    set(const OligoNumber p,
        const MatchPosition m,
        const uint8_t n) {
        position=p;
        matchPos=m;
        numErrors=n;
    }
    
    OligoNumber position;
    MatchPosition matchPos;
    uint8_t numErrors;
};


typedef std::vector<MatchStore> MatchCacheStore;
typedef MatchCacheStore::const_iterator MatchIter;



// struct MatchTable: this stores all the match information accumulated
// for each oligo during the search
class MatchTable
{
public:
    //  static MatchTable<OLIGO_LEN>* getMatchTable( const char* outputFileName );

  MatchTable( const int OLIGO_LEN, const char* outputFileName,const bool& write_multi=false,const short& no_of_seeds=1 ) :
    OLIGO_LEN_(OLIGO_LEN),
    no_of_seeds_(no_of_seeds),
    write_multi_(write_multi)
  {
      // initializing member variables
      read_length_ = 0;
      sensitive_ = false;

#if (NUM_THREADS>1)
    pthread_mutex_init(&mutex_, NULL);
#endif
    if( write_multi == true )
    {
        pOut_ = fopen( string(string(outputFileName)+".multi").c_str(), "w");
    } else
    {
        // open /dev/null to avoid creating the file
        pOut_ = fopen( "/dev/null", "w");
    }


    if (pOut_==NULL)
    {
      cerr << "Error in MatchTable: could not open file "
           << outputFileName << endl;
      exit (1);
    } else
      {
	outputFileName_ = outputFileName;
      }// ~if
  } // ~ctor
  virtual ~MatchTable()
  {
    fclose (pOut_);
  } // ~dtor

  // addMatch: process a set of matches
  // MatchIter iterates over MatchStore objects.
  //
  virtual void addMatch(MatchIter i, const MatchIter i_end) = 0;

  // print: print the search results to standard output at the end of a run
  virtual void print( OligoSource& oligos,
		      MatchPositionTranslator& getMatchPos,
		      const vector<string>& chromNames,
		      const vector<MatchPosition>& blockStarts,
		      const SuffixScoreTable& scoreTable,
		      int oligoLength )=0;

  virtual void printSquash( OligoSource& oligos,
			    MatchPositionTranslator& getMatchPos,
			    const vector<string>& chromNames,
			    const vector<MatchPosition>& blockStarts,
			    const SuffixScoreTable& scoreTable,
			    int oligoLength,
			    const string& directoryName,
			    const bool& align )=0;

  virtual bool getUnmappedReads( vector<bool>& unmapped ) =0;

  virtual bool mergeTable( MatchTable* source,MatchPositionTranslator& getMatchPos )=0;

    virtual bool buildMatchTable(MatchPositionTranslator& getMatchPos)=0;


  virtual void resize( int n )
  {
    matchPosition_.resize(n);
    matchType_.resize(n);
  } // ~void resize( int n )

  virtual void setQualityFailed__( int i )
  {
    if (matchPosition_.size()!=0)
    {
      assert((int)matchPosition_.size()>i);
      matchPosition_[i]=qualityFailed;
    }
  } // ~setQuality__

  virtual void setRepeatMasked__( int i )
  {
    matchPosition_[i]=repeatMasked;
  } // ~setRepeatMask__

  virtual bool isRepeatMasked__( int i ) const
  {
    return (matchPosition_[i]==repeatMasked);
  } // ~isRepeatMasked__

  virtual bool hasQualityFailed__( int i ) const
  {
    return (matchPosition_[i]==qualityFailed);
  } // ~hasQualityFailed__

  //  virtual bool hasQualityFailed__( int i ) const
  //  {
  //   return (matchPosition_[i]==qualityFailed);
  //  } // ~isRepeat__

  virtual bool isInterested__( int oligoNum, int PASS, bool hasNs );


  virtual void setSameAs__( OligoNumber newOligo, OligoNumber existingOligo )
  {
    if (  (matchPosition_[newOligo]<blockRepeat)
          && (matchPosition_[existingOligo]!=blockRepeat+newOligo) )
    {
      matchPosition_[newOligo] = blockRepeat+existingOligo;
    } // ~if
  } // ~setSameAs__
  size_t size( void ) const
  {
    assert(matchPosition_.size()==matchType_.size());
    return matchPosition_.size();
  } // ~size_t size( void ) const


  // adjustMatchPos: make small corrections to the match position to correct
  // for various quirks of the matching process (mainly to do with
  // incorporation errors and reverse complements)
  void adjustMatchPos( const char* pOligo,
		       const char dirChar,
		       char& firstN,
		       char& secondN,
		       MatchPosition& matchPos);

  // set the number of seeds
  void setNoOfSeeds( const short& no_of_seeds ){no_of_seeds_ = no_of_seeds;}


  void setReadLength( const short& read_length ){read_length_ = read_length;}


  void setSensitivity( const bool &sensitive ){sensitive_ = sensitive;}


protected:
  int OLIGO_LEN_;

  vector<MatchPosition> matchPosition_;
  vector<MatchDescriptor> matchType_;
  vector<uint> translator_;

  string outputFileName_;

  FILE* pOut_;

  short no_of_seeds_;
  short read_length_;
  bool sensitive_;

  bool write_multi_;

#if (NUM_THREADS>1)
  pthread_mutex_t mutex_;
#endif

}; // ~struct MatchTable



#if 0
// MatchTableVerbose: produce verbose match output
class MatchTableVerbose : public MatchTable
{
public:
MatchTableVerbose( const int OLIGO_LEN, const char* outputFileName ) :
    MatchTable( OLIGO_LEN, outputFileName ) {} // ~ctor

  virtual void print( OligoSource& oligos,
		      MatchPositionTranslator& getMatchPos,
		      const vector<string>& chromNames,
		      const vector<MatchPosition>& blockStarts,
		      const SuffixScoreTable& scoreTable,
		      int oligoLength);

  virtual void printSquash( OligoSource& oligos,
			    MatchPositionTranslator& getMatchPos,
			    const vector<string>& chromNames,
			    const vector<MatchPosition>& blockStarts,
			    const SuffixScoreTable& scoreTable,
			    int oligoLength,
			    const string& directoryName,
			    const bool& align){}

  virtual bool getUnmappedReads( vector<bool>& unmapped ){return true;}

   virtual bool mergeTable( MatchTable* source,MatchPositionTranslator& getMatchPos ){return true;}

  virtual bool buildMatchTable( MatchPositionTranslator& getMatchPos ){return true;}


}; // ~class MatchTableVerbose


// MatchTableAssembly : print output in a format compatible with the
// assembly module
class MatchTableAssembly : public MatchTable
{
public:
  MatchTableAssembly( const int OLIGO_LEN, const char* outputFileName ) :
    MatchTable( OLIGO_LEN, outputFileName ) {} // ~ctor

  virtual void print( OligoSource& oligos,
		      MatchPositionTranslator& getMatchPos,
		      const vector<string>& chromNames,
		      const vector<MatchPosition>& blockStarts,
		      const SuffixScoreTable& scoreTable,
		      int oligoLength);

  virtual void printSquash( OligoSource& oligos,
			    MatchPositionTranslator& getMatchPos,
			    const vector<string>& chromNames,
			    const vector<MatchPosition>& blockStarts,
			    const SuffixScoreTable& scoreTable,
			    int oligoLength,
			    const string& directoryName,
			    const bool& align ){}

  virtual bool getUnmappedReads( vector<bool>& unmapped ){return true;}

    virtual bool mergeTable( MatchTable* source,MatchPositionTranslator& getMatchPos ){return true;}

    virtual bool buildMatchTable( MatchPositionTranslator& getMatchPos ){return true;}



  // parseErrorInfo: parses the information stored in the error char
  // - bottom 6 bits are the error position
  // - using the top 2 bits and the original read it works
  //   out what the erroneous base should have been
  void parseErrorInfo( const char* oligoBuf,
                       const char errorInfo,
                       const char dirChar,
                       uint& errorPos,
                       uchar& shouldBe );


}; // ~class MatchTableAssembly
#endif

// MatchTableMulti : save multiple matches for each oligo in a temporary file
class MatchTableMulti : public MatchTable
{
public:
  enum
  {
    maxNumMatchesDefault=10
    //  maxNumMatchesOneError=10,
    //  maxNumMatchesTwoErrors=10
  };

  //  MatchTableMulti( const char* outputFileName,
  //		   const char* tmpFilePrefix=NULL );

  // Specify separate max numbers of 0,1,2-error matches
  MatchTableMulti( const int OLIGO_LEN,
                   const char* outputFileName,
                   const bool& write_multi,
		   const int maxNumMatchesExact,
		   const int maxNumMatchesOneError,
		   const int maxNumMatchesTwoErrors,
		   const char* tmpFilePrefix=NULL) :
    MatchTable( OLIGO_LEN,outputFileName,write_multi ),
    maxNumMatchesExact_(maxNumMatchesExact),
    maxNumMatchesOneError_(maxNumMatchesOneError),
    maxNumMatchesTwoErrors_(maxNumMatchesTwoErrors),
    matchesStored_(0)
  {
      initializeTmpFiles( tmpFilePrefix );
  } // ~ctor


  // Use default number of matches per read
  MatchTableMulti (const int OLIGO_LEN,
                   const char* outputFileName,
                   const bool& write_multi,
		   const char* tmpFilePrefix=NULL) :
    MatchTable( OLIGO_LEN, outputFileName,write_multi ),
    maxNumMatchesExact_(maxNumMatchesDefault),
    maxNumMatchesOneError_(maxNumMatchesDefault),
    maxNumMatchesTwoErrors_(maxNumMatchesDefault),
    matchesStored_(0)
  {
      initializeTmpFiles( tmpFilePrefix );
  } // ~ctor


  // Store same number of matches for each type of match
  MatchTableMulti (const int OLIGO_LEN,
                   const char* outputFileName,
                   const bool& write_multi,
		   const int maxNumMatches,
		   const char* tmpFilePrefix=NULL) :
    MatchTable( OLIGO_LEN,outputFileName,write_multi ),
    maxNumMatchesExact_(maxNumMatches),
    maxNumMatchesOneError_(maxNumMatches),
    maxNumMatchesTwoErrors_(maxNumMatches),
    matchesStored_(0)
  {
      initializeTmpFiles( tmpFilePrefix );
  } // ~ctor


  ~MatchTableMulti();

  // addMatch: process a set of matches
  // MatchIter iterates over MatchStore objects.
  //
  virtual void addMatch(MatchIter i, const MatchIter i_end);

  virtual void print( OligoSource& oligos,
		      MatchPositionTranslator& getMatchPos,
		      const vector<string>& chromNames,
		      const vector<MatchPosition>& blockStarts,
		      const SuffixScoreTable& scoreTable,
		      int oligoLength);


  virtual void printSquash( OligoSource& oligos,
			    MatchPositionTranslator& getMatchPos,
			    const vector<string>& chromNames,
			    const vector<MatchPosition>& blockStarts,
			    const SuffixScoreTable& scoreTable,
			    int oligoLength,
			    const string& directoryName,
			    const bool& align);

  virtual bool getUnmappedReads( vector<bool>& unmapped );

    bool getMatchInformation( vector< vector<MultiMatch> >& multimatches,vector<MatchDescriptor>& matchdescriptor );
    bool mergeTable( MatchTable* source,MatchPositionTranslator& getMatchPos );
  virtual  bool buildMatchTable( MatchPositionTranslator& getMatchPos );
    bool clear(void);


  const int maxNumMatchesExact_;
  const int maxNumMatchesOneError_;
  const int maxNumMatchesTwoErrors_;

  uint matchesStored_;
  FILE* pOligoNum_;
  FILE* pMatchType_;
  string tmpFileNameOligoNum_;
  string tmpFileNameMatchType_;
  vector<vector<MultiMatch> > multiMatch_;

  vector<bool> hyperhyper_;
  vector<bool> unmapped_;

  //  vector<vector<MatchPosition> > multiPos_;
  //  vector<vector<uchar> > multiType_;
private:
    // initialize does some setup that is common to all constructora
    void initializeTmpFiles( const char* tmpFilePrefix=NULL);

}; // ~class MatchTableMulti




// ------------------------ BEGIN -----------------------
// ---------- MatchTableMultiSquareSeed -----------------
// this match table handles multiple hits and multi seeds
// --> MultiSquared
class MatchTableMultiSquareSeed : public MatchTableMulti
{
public:
  // Specify separate max numbers of 0,1,2-error matches
  MatchTableMultiSquareSeed( const int OLIGO_LEN,
                       const char* outputFileName,
                       const bool& write_multi,
                       const int maxNumMatchesExact,
                       const int maxNumMatchesOneError,
                       const int maxNumMatchesTwoErrors,
                       const char* tmpFilePrefix=NULL) :
      MatchTableMulti( OLIGO_LEN,outputFileName,write_multi,maxNumMatchesExact,maxNumMatchesOneError,maxNumMatchesTwoErrors,tmpFilePrefix )
    {
    } // ~ctor


  // Use default number of matches per read
  MatchTableMultiSquareSeed (const int OLIGO_LEN,
                             const char* outputFileName,
                       const bool& write_multi,
		   const char* tmpFilePrefix=NULL) :
    MatchTableMulti( OLIGO_LEN, outputFileName,write_multi,MatchTableMulti::maxNumMatchesDefault,tmpFilePrefix )
  {
  } // ~ctor


  // Store same number of matches for each type of match
  MatchTableMultiSquareSeed (const int OLIGO_LEN,
                             const char* outputFileName,
                       const bool& write_multi,
		   const int maxNumMatches,
		   const char* tmpFilePrefix=NULL) :
    MatchTableMulti( OLIGO_LEN, outputFileName,write_multi,maxNumMatches,tmpFilePrefix)
  {
  } // ~ctor


    ~MatchTableMultiSquareSeed() {};

  virtual void resize( int n )
  {
    MatchTableMulti::matchPosition_.resize(n);
    MatchTableMulti::matchType_.resize(n);
    ms_matchType_.resize(4*n);

  }

  // addMatch: process a set of matches
  // MatchIter iterates over MatchStore objects.
  //
  virtual void addMatch(MatchIter i, const MatchIter i_end);

//  virtual void print( OligoSource& oligos,
//		      MatchPositionTranslator& getMatchPos,
//		      const vector<string>& chromNames,
//		      const vector<MatchPosition>& blockStarts,
//		      const SuffixScoreTable& scoreTable,
//		      int oligoLength);
//
//
//  virtual void printSquash( OligoSource& oligos,
//			    MatchPositionTranslator& getMatchPos,
//			    const vector<string>& chromNames,
//			    const vector<MatchPosition>& blockStarts,
//			    const SuffixScoreTable& scoreTable,
//			    int oligoLength,
//			    const string& directoryName,
//			    const bool& align);
//
//  virtual bool getUnmappedReads( vector<bool>& unmapped );

//    bool getMatchInformation( vector< vector<MultiMatch> >& multimatches,vector<MatchDescriptor>& matchdescriptor );
    inline
    bool checkNumberOfHits( const OligoNumber oligoNum, const uint8_t seedID, const uint8_t numErrors ) {
        return (((numErrors==0)
                 &&(ms_matchType_[4*oligoNum+seedID].r[0]<=this->maxNumMatchesExact_))
                || ((numErrors==1)
                    &&(ms_matchType_[4*oligoNum+seedID].r[0]<=this->maxNumMatchesExact_)
                    &&(ms_matchType_[4*oligoNum+seedID].r[1]<=this->maxNumMatchesOneError_))
                || ((numErrors==2)
                    &&(ms_matchType_[4*oligoNum+seedID].r[0]==0)
                    &&(ms_matchType_[4*oligoNum+seedID].r[1]<=this->maxNumMatchesOneError_)
                    &&(ms_matchType_[4*oligoNum+seedID].r[2]<=this->maxNumMatchesTwoErrors_)));
    }

    bool mergeTable( MatchTable* source,MatchPositionTranslator& getMatchPos );
    bool buildMatchTable( MatchPositionTranslator& getMatchPos );

    vector< MatchDescriptor > ms_matchType_;


  //  vector<vector<MatchPosition> > multiPos_;
  //  vector<vector<uchar> > multiType_;


}; // ~class MatchTableMultiSquareSeed


// ------------- END ----------------------
// --------- MatchTableMultiSquareSeed ----------

#if 0
MatchTable* MatchTable::getMatchTable( const char* outputFileName )
{
  const char suffixName[] =".vmf";
  const char* p(strstr(outputFileName, suffixName));
  if ((p==NULL)
      ||(p!=outputFileName
         +strlen(outputFileName)
         -strlen(suffixName)))
  {
    cerr << "Will output assembly format match details to file "
         << outputFileName << endl;
    return new MatchTableAssembly(outputFileName);
  } // ~if
  else
  {
    cerr << "Will output verbose format match details to file "
         << outputFileName << endl;
    return new MatchTableVerbose(outputFileName);
  } // ~else
}
#endif


inline
vector <int> calculateSeedOffsets(const int OLIGO_LEN, int readLength)
{
    BOOST_ASSERT( OLIGO_LEN<=readLength && "OLIGO_LEN must not be bigger than the read length");
    const int noUncoveredBases = (readLength - OLIGO_LEN);
    // take off half of OLIGO_LEN, then split the remaining read into four seeds
    const int baseLine = (noUncoveredBases/4)<(OLIGO_LEN/2)?(noUncoveredBases/4):(OLIGO_LEN/2); // place the first seed at OLIGO_LENGTH/2
    int moveAhead = (readLength - baseLine - OLIGO_LEN)/3; // let's use all the available four seeds
    if( moveAhead < 0 ) { moveAhead = 0; } // in the case that we don't

#ifdef DEBUG
  std::cerr << "calculateSeedOffsets: readLength = " << readLength
       << " noUncoveredBases = " << noUncoveredBases
       << " baseLine = " << baseLine
       << " moveAhead = " << moveAhead << std::endl;
#endif

    vector<int> ret;
    for( uint i = 0; i < 4; ++i )
    {
        ret.push_back(baseLine+i*moveAhead);
        BOOST_ASSERT( (ret[i]+OLIGO_LEN)<=readLength );
    }
    return ret;
}

} // namespace eland_ms
} // namespace casava



#endif // CASAVA_ELAND_MS_MATCH_TABLE_H
