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
 ** \file eland_ms/ContigNameFinder.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 ** Contains interface functionality from ELAND - specifically, stuff to do
 ** with scoring and storage of alignments
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_CONTIG_NAME_FINDER_H
#define CASAVA_ELAND_MS_CONTIG_NAME_FINDER_H

namespace casava
{
namespace eland_ms
{
// class ContigNameFinder
// If required, this converts a match position into a contig name + position
// in contig, else returns a null string and leaves position unchanged

class ContigNameFinder
{
 public:
 static ContigNameFinder* getContigNameFinder
   ( const string& directoryName, const string& chromName );
 virtual void operator() ( const char*& pContigName, uint& outputPos )=0;
 virtual ~ContigNameFinder() {}
};

class ContigNameFinderNull : public ContigNameFinder
{
 public:
  ContigNameFinderNull( void ): nullContigName('\0') {}
  virtual void operator() ( const char*& pContigName, uint& )
  {
    pContigName=&nullContigName;
  } // ~operator()
 private:
  const char nullContigName;
};

class ContigNameFinderIndex : public ContigNameFinder
{
 public:
  ContigNameFinderIndex( FILE* pIndex )
  {
    const uint bufSize(2048);
    char buf[bufSize];
    uint offset;

    int numItems;
    while( (numItems=fscanf( pIndex, "%u\t%s\n", &offset, buf))==2)
    {
      offsets_.push_back(offset);
      names_.push_back('/'+(string)(buf+1));
      //cerr << "Parsed offset=" << offsets_.back() << ", name="
      //	   << names_.back() << endl;
    } // ~while
    if (numItems!=EOF)
    {
	cerr << "Problem parsing line in index file :"
	     << " " <<buf << endl;
	exit (1);
    } // ~if
    fclose(pIndex);
  } // ~ctor

  virtual void operator() ( const char*& pContigName, uint& outputPos )
  {
    vector<uint>::iterator i
      (lower_bound(offsets_.begin(),offsets_.end(),outputPos));
    --i;
    outputPos-=*i;
    pContigName=names_[i-offsets_.begin()].c_str();
    // TBD
  } // ~operator()
 private:
  vector<uint> offsets_;
  vector<string> names_;

};

} // namespace eland_ms
} // namespace casava



#endif // CASAVA_ELAND_MS_CONTIG_NAME_FINDER_H
