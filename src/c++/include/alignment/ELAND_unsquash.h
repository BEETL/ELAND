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
// MODULE: unsquashGenome.cpp
// AUTHOR: A. J. Cox

// Description: grab sequence fragments from a directory of sequence files
// that have been 'squashed' by SquashGenome.cpp into the 2-bits-per-base
// format used by ELAND.
// For convenience, text files in ELAND output format are used as input.

//#define DEBUG

#ifndef INCLUDED_ELAND_UNSQUASH_H
#define INCLUDED_ELAND_UNSQUASH_H

#include <map>
#include <fstream>
#include <sys/resource.h>
#include <sys/time.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>


#include "alignment/GlobalUtilities.hh"


// bufSize: maximum size of an input line
const uint bufSize(2048);
// batchSize: this many requests will be stored then sorted and fulfilled
// in order of file name and position
const uint batchSize(262144);

// maxReadLengthELAND: longest read length allowable by ELAND
const uint maxReadLengthELAND(32);

// LineBuffer: read in ASCII text line-by-line, performing
// buffer resizing as necessary
// Simple usage:
//  LineBuffer lb(stdin);
//  lb();
//  while (lb())
//  {
//    printf ("%s", *lb);
//  }

// NB1: if last line of file does not end in a carriage return then neither
// will the output line
// NB2: does not check for Windows-format text files
// TBD move to GlobalUtilities.h if it proves useful elsewhere
class LineBuffer
{
public:
  LineBuffer( FILE* , int initialBufSize=2048) :
    bufSize_(initialBufSize)
  {
    buf_=new char[initialBufSize];
  } // ~ctor
  ~LineBuffer()
  {
    delete [] buf_;
  } // ~dtor
  char* operator*( void )
  {
    return buf_;
  } // ~operator*
  bool operator()( void )
  {
    buf_[bufSize_-1]='A'; // this can be anything apart from '\0'

    if (fgets(buf_,bufSize_,stdin)==NULL)
    {
      fprintf(stderr,"end of file reached\n");
      return false;
    } // ~if
    else if ((buf_[bufSize_-1]!='\0')||(buf_[bufSize_-2]=='\n'))
    {
      return true; // OK, fits in buffer
    } // ~else if
    else
    {
      do
      {
        int bufStart=bufSize_-1;
        bufSize_*=2;
        fprintf(stderr, "resizing to %d\n", bufSize_);
        char* tmpBuf=new char[bufSize_];
        strcpy(tmpBuf,buf_);
        delete [] buf_;
        buf_=tmpBuf;
        buf_[bufSize_-1]='A'; // this can be anything apart from '\0'
        if (fgets(buf_+bufStart,bufSize_-bufStart,stdin)==NULL)
        {
          return true; // if here then already got some data, so return that
        } // ~if
      } // ~do
      while ((buf_[bufSize_-1]=='\0')&&(buf_[bufSize_-2]!='\n'));
    } // ~else
    return true;
  } // ~operator*

private:
  char* buf_;
  int bufSize_;

}; // ~class LineBuffer


struct ContigIndex
{

  ContigIndex( const string& squashDirName, const char* name ) { open(squashDirName + "/" + name + ".idx"); }
  ContigIndex(const string& indexFilePath) { open(indexFilePath); }

  void open(const string& indexFilePath)
  {
    cerr << "Making index " << indexFilePath << endl;
    FILE* pIndex=fopen(indexFilePath.c_str(),"r");
    if (pIndex==NULL)
    {
      cerr << "Could not find index file " << indexFilePath << endl;
      exit (1);
    } // ~if
    const uint bufSize(2048);
    char buf[bufSize];
    uint offset;

    int numItems;
    //    LineBuffer lb(pIndex);
    while( (numItems=fscanf( pIndex, "%u\t%s\n", &offset, buf))==2)
    {
      offsets_.push_back(offset);
	names_.push_back((string)(buf+1));

#ifdef DEBUG
	//cerr << "Parsed offset=" << offsets_.back() << ", name="
	//     << names_.back() << endl;
#endif
    } // ~while

    if (numItems!=EOF)
    {
	cerr << "Problem parsing line in index file :"
	     << " " << buf << endl;
	exit (1);
    } // ~if

    fclose(pIndex);

    // populate map at end, avoids any problems with out of date pointers
    for (uint i(0);i<offsets_.size();++i)
    {
      //      index_[names_[i].c_str()]=offsets_[i];
      index_[names_[i].c_str()]=i;
    } // ~for
  } // ~ctor

  // Given a contig name, returns the contig's number (contigs are ordered as
  // presented in the originating file, numbering starts from zero) and
  // adjusts the position so it is relative to the start of the file not
  // the start of the contig
  void adjustPos (const char* contigName, uint& contigNum, uint& chromPos )
  {
      IndexT::iterator
      i(index_.find(contigName));
    if (i!=index_.end())
    {
      contigNum=i->second;
      chromPos+=offsets_[contigNum];
      //      chromPos+=i->second;
    } // ~if
    else
    {
      cerr << "Could not find entry for " << contigName << endl;
      exit (1);
    } // ~else
  } // ~adjustPos

    string getContigName( uint& contigNum, int& offset )
    {
        IndexT::iterator it;
        for( it=index_.begin();it!=index_.end();it++ )
        {
            if( (uint)it->second==contigNum )
            {
                offset = offsets_[contigNum];
                return it->first;
            }
        }
        offset = 0;
        return "";
    }


  struct LessThanString
  {
    bool operator()( const char* a, const char* b ) const
    {
      return (strcmp(a, b)<0);
    } // ~bool lessThanGenomeFrag
  };
  typedef map<const char*, uint, LessThanString > IndexT;

  vector<string> names_;
  vector<uint> offsets_;
  IndexT index_;

  //  int getIndex(const char* name)
  //  {
  //   if (index_.find(name)==index_.end())
  //  { // add name to index if not already seen
  //   names_.push_back(string(name));
  //  index_[names_.back().c_str()]=names_.size()-1;
  //  } // ~if
  //   return index_[name];
  //  } // ~get
}; // ~struct ContigIndex


// LessThanString: function object to enable STL algorithms to compare
// C-style strings. Was a member class of StringIndex

struct LessThanString
{
  bool operator()( const char* a, const char* b ) const
  {
    return (strcmp(a, b)<0);
  } // ~bool operator()
};



// StringIndex: assigns an index number to each different string seen.
// This index number can be used to retrieve the corresponding name.
// map<std::string,int> will
// do pretty much the same job but I want to avoid having to instantiate a
// std::string object for each line of input parsed. Std::strings are just
// used as convenient stores for const char* strings.
// TC 1.8.7 - modified to distinguish between 'chromName' and
// 'chromName/contigName' - in the latter case, adjust the match position
// from position in contig to position in file
struct StringIndex
{
  StringIndex( const string& squashDirName ) :
    squashDirName_(squashDirName)
  {
    cerr << "Building string index for " << squashDirName << endl;
    // Need to make sure the vector is not reallocated as all the pointers
    // then become invalid. Assuming each entry refers to a different file
    // is the most conservative assumption.
    names_.reserve(batchSize);
  } // ~ctor

  ~StringIndex()
  {
    for ( vector<ContigIndex*>::iterator i(contig_.begin());
	  i!=contig_.end(); ++i )
      delete *i;
  } // ~dtor



  //  int getIndex(const char* name)
  void getIndex(const char* name,
		uint& chromNum,
		uint& contigNum,
		uint& chromPos)
  {
    //    cerr << "getIndex " << name << endl;
    const char* contigName(NULL);
    char* nameBuf(NULL);
    if ((contigName=strchr(name,'/'))!=NULL)
    {
      nameBuf=new char[strlen(name)+1];
      strcpy(nameBuf,name);
      char * tmpContigName(nameBuf+(contigName-name));
      *tmpContigName='\0';
      contigName=tmpContigName+1;
      name=(const char*) nameBuf;
    } // ~if

    if (index_.find(name)==index_.end())
    { // add name to index if not already seen
      names_.push_back(string(name));
      index_[names_.back().c_str()]=names_.size()-1;
      contig_.push_back(NULL);
    } // ~if
    chromNum=index_[name];
    if (contigName!=NULL)
    { // specified chrom/contig, adjust contig pos to file pos
      //    cerr << "found contig name " << contigName << endl;
      if (contig_[chromNum]==NULL)
      {
	contig_[chromNum]= new ContigIndex
	  ( squashDirName_, name );
      } // ~if
      contig_[chromNum]->adjustPos( contigName, contigNum, chromPos );
    } // ~if
    else
    {
      contigNum=0;
    }
    if (nameBuf!=NULL) delete [] nameBuf;

  } // ~getIndex


    string getContigName(uint& chromNum,
                         uint& contigNum,
                         int& offset
        )
    {
        if (contig_[chromNum]!=NULL)
        {
            return string( contig_[chromNum]->getContigName( contigNum,offset ) );
        } // ~if

        return "";
    } // getContigName()


  const string& squashDirName_;
  vector<string> names_;
  vector<ContigIndex*> contig_;
  map<const char*, int, LessThanString > index_;

}; // ~struct StringIndex



struct SeqRequest
{
  SeqRequest( int requestNum,
	      int readNum,
	      int fileIndex,
	      int contigNum,
	      uint filePos,
	      char strand,
	      unsigned int seedOffset) :
    requestNum_(requestNum),
    readNum_(readNum),
    fileIndex_(fileIndex),
    contigNum_(contigNum),
    filePos_(filePos),
    strand_(strand),
    seedOffset_(seedOffset)
  {}
  int requestNum_;
  int readNum_;
  int fileIndex_;
  int contigNum_;
  uint filePos_;
  char strand_;
  unsigned int seedOffset_;
}; // ~struct SeqRequest

inline bool lessThanRequest (const SeqRequest& a, const SeqRequest& b)
{
  //  return ((a.fileIndex_<b.fileIndex_)
  //	  ||((a.fileIndex_==b.fileIndex_)&&(a.filePos_<b.filePos_)));
  return ((a.fileIndex_<b.fileIndex_)
	  ||((a.fileIndex_==b.fileIndex_)&&(a.contigNum_<b.contigNum_))
	  ||((a.fileIndex_==b.fileIndex_)
	     &&(a.contigNum_==b.contigNum_)
	     &&(a.filePos_<b.filePos_)));
} // ~bool lessThanRequest;

struct ContigValidRegion
{
  ContigValidRegion() :
    contigNum_(0),
    start_(0),
    finish_(0) {}
  ContigValidRegion
  ( int contigNum, long long start, long long finish ) :
    contigNum_(contigNum),
    start_(start),
    finish_(finish) {}
  int contigNum_;
  long long start_;
  long long finish_;
};



inline bool lessThanFinishContigValidRegion( const ContigValidRegion& a,
		const ContigValidRegion& b )
{
  return(a.finish_<b.finish_);
} // ~bool lessThanGenomeFrag

inline bool lessThanContigNumContigValidRegion( const ContigValidRegion& a,
		const ContigValidRegion& b )
{
  return(a.contigNum_<b.contigNum_);
} // ~bool lessThanGenomeFrag


// SquashFile: provide interface to a 'squashed' chromosome file. Assumes:
// First line might be a comment starting with '>', rest of lines are sequence.
// All sequence lines contain same number of bases, except possibly the last.
// NB Numbering of bases in file starts at zero.
class SquashFile
{
public:
  //  SquashFile( string name, const StringIndex& files );
  ~SquashFile() { munmap( pMap, fileSize ); }

  uint getNextBase( long long i)
  {
    return ( (pWord[i>>4])>>(2*((i&0xF)^0xF)) ) & 0x3;
  } // ~getNextBase

  char operator[]( long long pos )
  {
    //if (pos<0) return 'N';
    //    assert(pos<=numChars);
    //  --pos;
    //  pos+=(pos/lineLength);
    // return pStart[pos];
    return(baseNames[getNextBase(pos)]);
  } // ~operator[]

  void goToPos( int contig, long long pos )
  {
#ifdef DEBUG
    cerr << "goToPos: " << contig << " " << pos << endl;
#endif

    if (pos_.contigNum_!=contig)
    { // Get the set of regions (=start/finish pairs) for desired contigs
      pos_.contigNum_=contig;
      contigRegions_=equal_range
	(valid_.begin(),valid_.end(),pos_,lessThanContigNumContigValidRegion);
#ifdef DEBUG
      for (vector<ContigValidRegion>::iterator i(contigRegions_.first);
	   i!=contigRegions_.second; ++i)
      {
	cerr << "contig region: "
	     << i->contigNum_ << " "
	     << i->start_ << " "
	     << i->finish_ << endl;
      }
#endif
    } // ~if

    pos_.finish_=pos;
    // Get the first region in the contig witrh a finish position that exceeds
    // pos
    thisRegion_=lower_bound
      (contigRegions_.first,contigRegions_.second,
       pos_,lessThanFinishContigValidRegion);
  }

  char getNextBase(void)
  {
    char base;
    if (thisRegion_==contigRegions_.second)
    { // have gone off the end of the region
      base= 'N';
    } // ~if
    else if (pos_.finish_<thisRegion_->start_)
    { // not yet at start of region
      base='N';
    } // ~else if
    else
    {
      base=(*this)[pos_.finish_];
    } // ~else
    pos_.finish_++;
    // do we need to change regions?
    if ((thisRegion_!=contigRegions_.second)
	&&(pos_.finish_>thisRegion_->finish_))
    {
      thisRegion_++;
    } // ~if
    return base;
  } // ~operator[]


  int getNumChars( void ) const { return numChars; }
private:
  void *pMap;
  const char* pStart;
  const char* pEnd;
  const Word* pWord;
  size_t fileSize;
  int lineLength;
  int numChars;
  vector<ContigValidRegion> valid_;
  pair<vector<ContigValidRegion>::iterator,
       vector<ContigValidRegion>::iterator > contigRegions_;
  vector<ContigValidRegion>::iterator thisRegion_;
  ContigValidRegion pos_;

public:
//SquashFile::SquashFile( string name, const StringIndex& files )
SquashFile( string squashDirName,
			string chromFileName,
			const StringIndex& files ) :
  pos_(-1,-1,-1)
{
  // populate list of valid regions
  string name(squashDirName+(string)"/"+chromFileName);



  FileReader fileReader(name.c_str());
  for (const ValidRegion* i(fileReader.getFirstValid());
       i!=fileReader.getLastValid(); ++i)
  {
    valid_.push_back(ContigValidRegion(0, i->start, i->finish));
#ifdef DEBUG
    cerr << "valid region: " << valid_.back().start_ << " "
	 << valid_.back().finish_ << endl;
#endif
  } // ~for

  // try to relate them to contigs
  map<const char*, int, LessThanString>::const_iterator i
    (files.index_.find(chromFileName.c_str()));

  if (i==files.index_.end())
  {
    cerr << "Could not find entry in index for " << name << endl;
  }
  else
  {
    ContigIndex* pContigIndex(files.contig_[i->second]);
    if (pContigIndex==NULL)
    {
      cerr << "No contig index found for " << name << endl;
    }
    else
    {
      thisRegion_=valid_.begin();
      int contigNum(0);
      // can assume at least two offsets, else squash would have been done in
      // single-contig mode (ie no .idx file at all)
      vector<uint>::iterator i(pContigIndex->offsets_.begin());
      while (++i!=pContigIndex->offsets_.end())
      {
	while((thisRegion_!=valid_.end())&&(thisRegion_->finish_<=*i))
	{
	  thisRegion_->contigNum_=contigNum;
	  thisRegion_++;
	} // ~while
	contigNum++;
      } // ~while
      for (;thisRegion_!=valid_.end();++thisRegion_)
      {
	thisRegion_->contigNum_=contigNum;
      }
#ifdef DEBUG
      for (thisRegion_=valid_.begin();thisRegion_!=valid_.end();thisRegion_++)
      {
	cerr << "Region: " << thisRegion_->contigNum_ << " " << thisRegion_->start_ << " " << thisRegion_->finish_ << endl;
      }
#endif
    } // ~else
  } // ~else


  // memory map sequence file
  name+=(string)".2bpb";

  int fd;
  fd = open(name.c_str(), O_RDONLY, S_IRUSR );
  if (fd==-1)
  {
    cerr << "Error: could not open file "
         << name << endl;
    exit (1);
  } // ~if
  fileSize = lseek(fd, 0, SEEK_END);

  pMap=mmap(0, fileSize, PROT_READ, MAP_SHARED, fd, 0);
  if (pMap==MAP_FAILED)
  {
    cerr << "Error: could not memory map file "
         << name << ": " << strerror(errno) << endl;
    exit (1);
  } // ~if
  close (fd);

  pStart=(char*)pMap;
  pEnd=pStart+fileSize;

  pWord = (const Word*) pStart;

} // ~ctor

}; // ~class SquashFile


struct NInfo
{
  static const uint initSize_ = (uint)-1;
  NInfo( void ) :
  headSize_(initSize_) /*, tailSize_(initSize_)*/ {}
  uint headSize_;
  //uint tailSize_;
}; // struct NInfo

struct FragmentFinder
{

  FragmentFinder ( const string& squashDirName,
		   const int readLength,
		   const int fragmentLength,
		   const int reverseStrandStartOffset ) :
  readLength_(readLength),
  fragmentLength_(fragmentLength),
  bases_ahead_((fragmentLength-readLength)/2),
  squashDirName_(squashDirName),
  reverseStrandStartOffset_(reverseStrandStartOffset)
  { assert(fragmentLength>=readLength); } // ~c'tor

private:
  const int readLength_;
  const int fragmentLength_;
  const int bases_ahead_;
  const string squashDirName_;
  const int reverseStrandStartOffset_;
  //  vector<int> numNs_;
  vector<NInfo> numNs_;

// default (= base class) behaviour of fetchFragment is to pull
// genomic fragments from the squash files, orienting them in the
// same direction as the read
void fetchFragment
(SquashFile& squash, const SeqRequest& i, const vector<char*>& reads, char* buf)
{
      buf[fragmentLength_]=0;
      if (NInfo::initSize_ == numNs_[i.readNum_].headSize_)
      {
	numNs_[i.readNum_].headSize_ = countHeadNs(reads[i.readNum_],readLength_);

        // Only count the leading Ns that overlap the seed
        if (i.seedOffset_ > 0)
        {
          numNs_[i.readNum_].headSize_ = (numNs_[i.readNum_].headSize_ <= i.seedOffset_ ? 0 : numNs_[i.readNum_].headSize_ - i.seedOffset_);
        }
      } // ~if

      if (i.strand_=='F')
      {
	squash.goToPos
	  (i.contigNum_,static_cast<long long>(i.filePos_)-numNs_[i.readNum_].headSize_-1-bases_ahead_);
	for (int j(0);j<fragmentLength_;j++)
	{
	  buf[j]=squash.getNextBase();
	  //	  buf[j]=(*squash_)[i.filePos_+j-numNs_[i.readNum_].headSize_-1];
	} // ~for j
      } // ~if
      else
      {

	squash.goToPos(i.contigNum_,
			  static_cast<long long>(i.filePos_)
			  +numNs_[i.readNum_].headSize_
			  -reverseStrandStartOffset_-1-bases_ahead_);

	for (int j(0);j<fragmentLength_;j++)
	{

	  buf[fragmentLength_-1-j]
	    =reverseCharASCII[(uint)(squash.getNextBase()) ];
	  //	  buf[readLength_-1-j]
	  //=reverseCharASCII[(uint)(*squash_)
	  //	      [i.filePos_+j
	  //	       +numNs_[i.readNum_].headSize_
	  //	       -reverseStrandStartOffset_-1]];
	} // ~for
      } // ~else
      //      cout << reads[i.readNum_] << " " << i.filePos_
      // << i.strand_ << " " << buf << endl;
} // ~FragmentFinder::fetchFragment


public:
void operator()
( vector<SeqRequest>& requests,
  const vector<char*>& reads,
  vector<char*>& frags,
  const StringIndex& files)
{
  if (requests.empty()) return;

  std::auto_ptr<SquashFile> pSquash;
  int lastFile=requests[0].fileIndex_+1; // ensure false on first go
  sort( requests.begin(), requests.end(), lessThanRequest );

  numNs_.clear();
  numNs_.resize( reads.size());

  typedef std::vector<SeqRequest>::iterator sriter;
  sriter i(requests.begin()),i_end(requests.end());
  for (;i!=i_end;++i){
    if (i->fileIndex_!=lastFile)
    {
	//	squashFileName=squashDirName_
	//  +(string)"/"
	//  +files.names_[i->fileIndex_]
	//  +(string)".2bpb";

	//	squashFileName=squashDirName_
	//	  +(string)"/"
	//  +files.names_[i->fileIndex_];

      pSquash.reset(new SquashFile(squashDirName_,
                                   files.names_[i->fileIndex_],
                                   files));
    } // ~if
    lastFile=i->fileIndex_;

#ifdef DEBUG
      cerr << "fetchFragment: DUMP SeqRequest = requestNum/readNum/fileIndex/contigNum/filePos/strand = "
	   << i->requestNum_ << "/"
	   << i->readNum_ << "/"
	   << i->fileIndex_ << "/"
	   << i->contigNum_ << "/"
	   << i->filePos_ << "/"
	   << i->strand_ << endl;
#endif

    fetchFragment(*(pSquash.get()), *i, reads, frags[i->requestNum_]);
  } // ~for i
} // ~FragmentFinder::operator()
}; // ~struct FragmentFinder

#if 0
// struct FragmentFinderAlign : public FragmentFinder
// Instead of printing the genomic fragment, replace bases that match
// the read with '-' characters
struct FragmentFinderAlign : public FragmentFinder
{

  FragmentFinderAlign ( const string& squashDirName,
		   const int readLength,
                   const int fragmentLength,
		   const int reverseStrandStartOffset ) :
  FragmentFinder( squashDirName,
		  readLength,
                  fragmentLength,
		  reverseStrandStartOffset )
  {} // ~c'tor

protected:
virtual void fetchFragment
( SquashFile& squash, const SeqRequest& i, const vector<char*>& reads, char* buf)
{
  FragmentFinder::fetchFragment( squash, i, reads, buf);

#ifdef DEBUG // debug flag markus
  cerr << endl << endl;
  cerr << "XXXXXXXXXXXXX" << endl;
  cerr << "i.readNum_ = " << i.readNum_ << "\t"
       << "readLength_ = " << readLength_ << "\t" << "reads = " << reads[i.readNum_] << "\tbuf = " << buf << endl;
#endif


  for (int j(0);j<readLength_;j++)
  {
    if (toupper(reads[i.readNum_][j])==buf[j+bases_ahead_])
      buf[j]='-';
  } // ~for j

} // ~FragmentFinderAlign::fetchFragment

}; // struct ~FragmentFinderAlign

// struct FragmentFinderAlignCompact : public FragmentFinderAlign
// Replace runs of matching bases ('-' characters) with an integer,
// eg '---A-----C' becomes '3A5C'
// Perl 1 liner to convert back:
// cat s_3_frag.txt | perl -lane '{ while (/(\d+)/) { $a="-" x $1; s/$1/$a/ } print }'
struct FragmentFinderAlignCompact : public FragmentFinderAlign
{

  FragmentFinderAlignCompact ( const string& squashDirName,
		   const int readLength,
		   const int fragmentLength,
		   const int reverseStrandStartOffset ) :
  FragmentFinderAlign( squashDirName,
		  readLength,
		  fragmentLength,
		  reverseStrandStartOffset )
  {} // ~c'tor
  char buf_[bufSize];

protected:
virtual void fetchFragment
( SquashFile& squash, const SeqRequest& i, const vector<char*>& reads, char* buf)
{
  FragmentFinderAlign::fetchFragment( squash, i, reads, buf_);
  int exactRun(0), adjustedLength(strlen(buf_));
  char* p(buf);
  char c;
  //  cout << adjustedLength << endl;
  for (int j(0);j<adjustedLength;j++)
  {
    if (buf_[j]=='-')
    {
      exactRun++;
    } // ~if
    else
    {
      c=buf_[j];
      if (exactRun!=0)
      {
	sprintf(p,"%u%c",exactRun,c);
	exactRun=0;
      } // ~if
      else
      {
	sprintf(p,"%c",c);
      } // ~else
      p+=strlen(p);
    } // ~else
  } // ~for j
  if (exactRun!=0) sprintf(p,"%u",exactRun);
} // ~FragmentFinderAlign::fetchFragment


}; // struct ~FragmentFinderAlign
#endif


struct ResultsParser
{
  virtual ~ResultsParser() {}
  virtual void operator()( const char* buf,
			   StringIndex& files,
			   vector<char*>& reads,
			   vector<SeqRequest>& requests ) const=0;
}; // ~struct ResultsParser

//struct ResultsParserStandard : public ResultsParser
//{
//  virtual void operator()( const char* buf,
//			   StringIndex& files,
//			   vector<char*>& reads,
//			   vector<SeqRequest>& requests ) const
//  {
//    char chromName[bufSize];
//    char readBuf[bufSize];
//    uint chromNum;
//    uint contigNum;
//    uint chromPos;
//    uint readLength;
//    char strand;
//    if (sscanf(buf, "%*s\t%s\t%*s\t%*s\t%*s\t%*s\t%s\t%u\t%c",
//	       readBuf, chromName, &chromPos, &strand )==4)
//    {
//      files.getIndex(chromName, chromNum, contigNum, chromPos );
//
//      requests.push_back
//	(SeqRequest(requests.size(),
//		    reads.size(),
//		    chromNum,
//		    contigNum,
//		    chromPos,
//		    strand));
//      readLength=strlen(readBuf);
//      reads.push_back( new char[readLength+1] );
//      strcpy(reads.back(),readBuf);
//      // cout << reads.back() << endl;
//    } // ~if
//    return;
//  } // ~operator()
//
//}; // ~struct ResultsParserStandard : public ResultsParser
//
//struct ResultsParserMulti : public ResultsParser
//{
//  virtual void operator()( const char* buf,
//			   StringIndex& files,
//			   vector<char*>& reads,
//			   vector<SeqRequest>& requests ) const
//  {
//    char chromName[bufSize];
//    const char *p1, *p2;
//    uint chromNum;
//    uint contigNum;
//    uint pos;
//    char strand;
//
//    p1=strchr(buf,'\t'); // scan past read ID
//    assert(p1!=NULL);
//    p1++;
//    p2=strchr(p1,'\t'); // scan past read
//    assert(p2!=NULL);
//    uint readLength(p2-p1);
//    p2++; // now first character of match type field
//    if (isalpha(*p2)) return; // NM,QC,RM...always a letter if no match found
//
//    reads.push_back( new char[readLength+1] );
//    strncpy(reads.back(),p1,readLength);
//    reads.back()[readLength]='\0'; // TBD needed? maybe done by strncpy?
//
//    p1=strchr(p2,'\t'); // scan past match type
//    //    cout << buf << endl;
//    if (p1==NULL) return; // too many repeats for matches to be stored
//
//    chromName[0]='\0';
//    while (1)
//    {
//      p1++;
//      p2=strpbrk(p1,":,");
//      if ((p2==NULL)||(*p2==','))
//      {
//	assert (sscanf(p1,"%u%c%*c", &pos, &strand)==2);
//	assert ((strand=='F')||(strand=='R'));
//	assert (chromName[0]!='\0');
//	files.getIndex(chromName, chromNum, contigNum, pos );
//	requests.push_back
//	  (SeqRequest(requests.size(),
//		      reads.size()-1,
//		      chromNum,
//		      contigNum,
//		      pos,
//		      strand));
//	if (p2==NULL) return;
//      } // ~if
//      else
//      { // must be a colon, therefore it's a chrom name
//	strncpy(chromName,p1,p2-p1);
//	chromName[p2-p1]='\0';
//      } // ~else
//      p1=p2;
//    } // ~while
//
//  }   // ~operator()
//
//}; // ~struct ResultsParserMulti : public ResultsParser
//

inline void printFragments( const vector<char*>& frags,
		     const vector<int>& requestsPerLine )
{
      vector<char*>::const_iterator j(frags.begin());
      for ( vector<int>::const_iterator i(requestsPerLine.begin());
	    i!=requestsPerLine.end(); ++i)
      {
	for (int k(0); k<*i; k++)
	{
	  if (k!=0) printf(":");
	  assert(j!=frags.end());
	  printf("%s",*j++);
	} // ~for k
	printf("\n");
      } // ~for i

} // ~void printFragments


#endif
