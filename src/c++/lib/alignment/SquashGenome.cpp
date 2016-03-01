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
 ** \file SquashGenome.cpp
 **
 ** \brief Implementation of squash genome apis.
 **
 ** \author Tony Cox
 **/

#include <dirent.h>
#include <set>
#include <string>

#include "alignment/GlobalUtilities.hh"

#include "alignment/ELAND_unsquash.h"

// SquashGenome.cpp
// Takes a set of DNA sequence files in ASCII format and squashes them
// into a 2 bits per base format. Produces two files for each sequence
// file.
// originalFileName.2bpb - compressed sequence data
// originalFileName.vld - gives the start and end positions
// of each valid region in the sequence file (i.e. regions containing
// only AGCT - can't handle '-'s, 'N's etc. using 2bpb format)
// Can also unsquash a single compressed file.

namespace casava
{
namespace alignment
{

enum
{
  seqLineLength = 60
};


const std::string getContigNameForbiddenCharacters()
{
    const static std::string invalidChars("?()[]/\\=+<>:;\"',*^&");
    return invalidChars;
}

// class IndexMaker
// creates an index to the entries of a multiple-entry fasta file
class IndexBuilder
{
    int logLevel_;
  public:
  IndexBuilder( const string& idxFileName, bool validateNames, int logLevel ) :
      logLevel_(logLevel), idxFile_(idxFileName.c_str()), validateChromNames_(validateNames)
  {
    // open output files
    if (idxFile_.fail())
    {
      if (logLevel > 0) cerr << "ERROR: Error in squash: could not open file " << idxFileName << endl;
      exit (1);
    } // ~if
    else
    {
      if (logLevel > 1) cerr << "Opened file " << idxFileName << endl;
    } // ~else
  } // ~ctor

  ~IndexBuilder()
  {
    idxFile_.close();
  } // ~dtor

  void addEntry( const string& strEntry, off_t entryPos)
  {
    string thisName(strEntry.substr(0, strEntry.find_first_of(" \t\n\r")));

    if (thisName.empty())
    {
       if (logLevel_ > 0) cerr << "ERROR: empty entry name in fasta file. "
          "Fasta headers must have a non-whitespace character after the '>'\n";
       exit (1);
    }
    else if (validateChromNames_ &&
            std::string::npos != thisName.find_first_of(getContigNameForbiddenCharacters()))
    {
        if (logLevel_ > 0) cerr << "ERROR: invalid entry name in fasta file ("
        << thisName << ") - cannot contain the following characters: " << getContigNameForbiddenCharacters()
        << endl;
        exit (1);
    }
    else if (names_.count(thisName)!=0)
    {
       if (logLevel_ > 0) cerr << "ERROR: duplicate entry name in fasta file ("
	   << thisName << ") - each entry must have distinct name"
	   << endl;
      exit (1);
    } // ~else if

    if (!(idxFile_ << entryPos << "\t>" << thisName << endl))
    {
        if (logLevel_ > 0) cerr << "ERROR: unable to write idx file\n";
        exit (2);
    }
    names_.insert(thisName);
  } // ~addEntry

  unsigned int getNumEntries( void ) const
  {
    return names_.size();
  } // ~getNumEntries( void ) const

private:
  std::ofstream idxFile_;
  std::set<std::string> names_;
  const bool validateChromNames_;
};

Word getNextBase( Word* pWord, unsigned int i)
{
  return ( (pWord[i>>4])>>(2*((i&0xF)^0xF)) ) & 0x3;
}

void outputSizesToXML( const char* dirName, DIR* pDir , int logLevel)
{
    static const std::string suffixName(".vld");
    const fs::path vldPath(dirName);

    std::string vldFileName, vldFileFullPath;

    FILE* pVldFile;
    dirent* dirEntry;

    std::cout << "<sequenceSizes>\n";

    while( (dirEntry = readdir(pDir)) )
    {
        vldFileName = dirEntry->d_name;
        if (std::equal(suffixName.rbegin(), suffixName.rend(), vldFileName.rbegin()))
        {
            fs::path vldFileFullPath = vldPath / vldFileName;

            // memory map validity file
            if ((pVldFile = fopen(vldFileFullPath.string().c_str(), "r"))==NULL)
            {
                if (logLevel > 0) cerr << "ERROR: Error in outputSizesToXML: could not open file "
                << vldFileFullPath << endl;
                exit (1);
            } // ~if

            char headerChar(0);
            while (1 == fread( &headerChar, sizeof(headerChar), 1, pVldFile) && '\n' != headerChar) ;

            std::string fileName(vldFileName.substr(0, vldFileName.length() - suffixName.length()));
            ContigIndex index(dirName, fileName.c_str());
            std::vector<unsigned int>::size_type contig(1);
            ValidRegion lastValidRegion, validRegion;

            while (1 == fread( &validRegion, sizeof(validRegion), 1, pVldFile))
            {
                if (contig < index.offsets_.size() && validRegion.finish >= index.offsets_[contig])
                {
                    std::cout << "\t<chromosome fileName=\""
                        << fileName << "\" contigName=\""
                        << index.names_[contig-1] << "\" totalBases=\""
                        << lastValidRegion.finish - index.offsets_[contig - 1] + 1 << "\"/>\n";
                    ++contig;
                }
                lastValidRegion = validRegion;
                if (logLevel > 2) cerr << lastValidRegion.start << "\t" << lastValidRegion.finish << "\t>" << index.names_[contig-1] << "\n";
            }
            std::cout << "\t<chromosome fileName=\""
                << fileName << "\" contigName=\""
                << index.names_[contig-1] << "\" totalBases=\""
                << lastValidRegion.finish - index.offsets_[contig - 1] + 1 << "\"/>\n";
        } // ~if
    } // ~while

  std::cout << "</sequenceSizes>\n";

} // ~outputSizesToXML


void unsquash( const char* squashName , int logLevel)
{
  string seqFileName(squashName), vldFileName(squashName);
  seqFileName += ".2bpb";
  vldFileName += ".vld";

  void *pMapSeq;
  void *pMapVld;
  int fd;
  int seqFileSize, vldFileSize;

  Word* pWord;
  //Word* pLastWord;

  ValidRegion* pValid;
  ValidRegion* pLastValid;

  char* pChar;


  // memory map sequence file
  fd = open(seqFileName.c_str(), O_RDONLY, S_IRUSR );
  if (fd==-1)
  {
    if (logLevel > 0) cerr << "ERROR: Error in unsquash: could not open file " << seqFileName << endl;
    exit (1);
  } // ~if
  seqFileSize = lseek(fd, 0, SEEK_END);
  if (logLevel > 1) cerr << "squash: opened file " << seqFileName << " of size "
       << seqFileSize << " bytes." << endl;

  if (0 == seqFileSize)
  {
    close(fd);
    return;
  }
  pMapSeq=mmap(0, seqFileSize, PROT_READ, MAP_SHARED, fd, 0);
  if (pMapSeq==MAP_FAILED)
  {
    if (logLevel > 0) cerr << "ERROR: Error in unsquash: could not memory map file "
	 << seqFileName << ": " << strerror(errno) << endl;
    exit (1);
  } // ~if
  close (fd);

  // memory map validity file
  fd = open(vldFileName.c_str(), O_RDONLY, S_IRUSR );
  if (fd==-1)
  {
    if (logLevel > 0) cerr << "ERROR: Error in unsquash: could not open file " << vldFileName << endl;
    exit (1);
  } // ~if
  vldFileSize = lseek(fd, 0, SEEK_END);
  if (logLevel > 1) cerr << "squash: opened file " << vldFileName << " of size "
       << vldFileSize << " bytes." << endl;

  pMapVld=mmap(0, vldFileSize, PROT_READ, MAP_SHARED, fd, 0);
  if (pMapVld==MAP_FAILED)
  {
    if (logLevel > 0) cerr << "ERROR: Error in unsquash: could not memory map file "
	 << vldFileName << ": " << strerror(errno) << endl;
    exit (1);
  } // ~if
  close (fd);

  // set up pointers for Word

  pChar= (char*) pMapSeq;
  pWord= (Word*) pMapSeq;
  //pLastWord = (Word*)(pChar+seqFileSize);

  // set up last pointer for valid region

  pChar= (char*) pMapVld;
  pLastValid = (ValidRegion*)(pChar+vldFileSize);

  // print out the comment line ...
  while (*(pChar++) != '\n') ;

  // set up
  pValid = (ValidRegion*) pChar;
  assert( pValid!=pLastValid );

  ContigIndex index((std::string(squashName) + ".idx").c_str());

  off_t i(0);
  off_t printedPos(0);
  int contig(0);
  while (pValid != pLastValid)
  {
      if (logLevel > 2) cerr << pValid->start << "\t" << pValid->finish << "\t>" << index.names_[contig] << "\n";
      if (pValid->finish >= index.offsets_[contig])
      {
          cout << (printedPos ? "\n>" : ">") << index.names_[contig++];
          printedPos = 0;
      }

      for (; i < pValid->start; ++i)
      {
          printedPos++ % seqLineLength || cout.put('\n');
          assert( !getNextBase(pWord, i));
          cout.put('N');
      }

      for (; i < (pValid->finish + 1); ++i)
      {
          printedPos++ % seqLineLength || cout.put('\n');
          cout.put(baseNames[getNextBase(pWord, i)]);
      }
      ++pValid;
  }

} // ~unsquash

void squash( const char* directoryName, const char* fileName, bool validateNames, bool allowManyContigs, int logLevel)
{
  // make file names
  if (logLevel > 1) cerr << "Full file path: " << fileName << endl;
  const char* pFileName( strrchr(fileName,'/') );
  if (pFileName==NULL) pFileName=fileName;
  else ++pFileName;
  if (logLevel > 1) cerr << "Extracted file name:" << pFileName << endl;

  string seqFileName(directoryName),
    vldFileName(directoryName),
    idxFileName(directoryName);

  seqFileName += '/';
  vldFileName += '/';
  idxFileName += '/';

  seqFileName += (string)pFileName;
  vldFileName += (string)pFileName;
  idxFileName += (string)pFileName;

  seqFileName += (string)".2bpb";
  vldFileName += (string)".vld";
  idxFileName += (string)".idx";


  if (logLevel > 1) cerr << seqFileName << " " << vldFileName << endl;

  // open output files
  ofstream seqFile(seqFileName.c_str());
  if (seqFile.fail())
  {
    if (logLevel > 0) cerr << "ERROR: Error in squash: could not open file " << seqFileName << endl;
    exit (1);
  } // ~if

  ofstream vldFile(vldFileName.c_str());
  if (vldFile.fail())
  {
    if (logLevel > 0) cerr << "ERROR: Error in squash: could not open file " << seqFileName << endl;
    exit (1);
  } // ~if

  // memory map input file
  std::ifstream reference( fileName );

  if (reference.fail())
  {
    if (logLevel > 0) cerr << "ERROR: Error in squash: could not open file " << fileName << endl;
    exit (1);
  } // ~if

  char ch(0);
  std::string lastHeader;
  if (!reference.get(ch) || '>' != ch || !std::getline(reference, lastHeader))
  {
    if (logLevel > 0) cerr << "ERROR: Error in squash: could not read fasta header " << fileName << endl;
    exit (1);
  } // ~if
  vldFile << lastHeader << '\n';


  off_t thisPos(0), totalValid(0);
  vector<Word> words;
  vector<ValidRegion> valids;
  bool inValidRegion(false);
  Word thisBase;

  IndexBuilder indexBuilder(idxFileName, validateNames, logLevel);
  indexBuilder.addEntry(lastHeader, thisPos );


  while(reference >> ch)
  {
    if ('\n' == ch) continue;

    if ('>' == ch)
    {
      if (!allowManyContigs)
      {
          if (logLevel > 0) cerr << "ERROR: Error in squash: multiple contigs are not allowed in: " << fileName << "\n";
          exit (1);
      }

      if (inValidRegion==true)
      {
        valids.back().finish=thisPos-1;
        totalValid+=valids.back().finish-valids.back().start+1;
        inValidRegion=false;
      } // #if
      else
      {
          // when the chromosome ends with non ACGT characters, the length does not include those.
          // add empty region to compensate for it
          valids.push_back(ValidRegion(thisPos, thisPos - 1));
      }
      if (logLevel > 2) cerr << valids.back().start << "\t" << valids.back().finish << "\t>" << lastHeader << '\n';
      if (!std::getline(reference, lastHeader))
      {
        if (logLevel > 0) cerr << "ERROR: Error in squash: empty '>' line found in fasta file " << fileName << endl;
        exit (1);
      }
      indexBuilder.addEntry( lastHeader, thisPos );
    }
    else
    {
        if ((thisPos%maxBasesPerWord)==0)
        {
          words.push_back(0);
        } // ~if

        thisBase = whichBase[(unsigned int)ch];

        if (thisBase!=nv)
        {
          if (inValidRegion==false)
          {
        valids.push_back(ValidRegion());
        valids.back().start=thisPos;
        inValidRegion=true;
          } // ~if
        } // ~if
        else
        {
          thisBase=0; // codes for A
          if (inValidRegion==true)
          {
        valids.back().finish=thisPos-1;
        if (logLevel > 2) cerr << valids.back().start << "\t" << valids.back().finish << "\t>" << lastHeader << '\n';
        totalValid+=valids.back().finish-valids.back().start+1;
        inValidRegion=false;
          } // #if
        } // ~else
        words.back() <<= numBitsPerBase;
        words.back() |= thisBase;
        ++thisPos;
    }
  }

  if (inValidRegion==true)
  {
    valids.back().finish=thisPos-1;
    totalValid+=valids.back().finish-valids.back().start+1;
  }
  else
  {
      // when the chromosome ends with non ACGT characters, the length does not include those.
      // add empty region to compensate for it
      valids.push_back(ValidRegion(thisPos, thisPos - 1));
  }

  if (logLevel > 2) cerr << valids.back().start << "\t" << valids.back().finish << "\t>" << lastHeader << "\n";
  words.back() <<= 2*(16-thisPos%maxBasesPerWord);

  // write sequence file
  if (!seqFile.write( (const char*)&words[0], words.size()*sizeof(Word))) {
      if (logLevel > 0) cerr << "ERROR: could not store squashed sequence in" << seqFileName << "\n";
      exit (2);
  }

  // write validFile
  if (!vldFile.write( (const char*)&valids[0], valids.size()*sizeof(ValidRegion))) {
      if (logLevel > 0) cerr << "ERROR: could not region data in" << vldFileName << "\n";
      exit (2);
  }

  if (logLevel > 0) cerr << "INFO: finished file " << fileName << endl
       << thisPos << " bases" << endl
       << totalValid << " valid bases ("
       << (100.*totalValid/thisPos) << "\045)" << endl
       << valids.size() << " valid regions" << endl
       << indexBuilder.getNumEntries() << " entries" << endl;
}

} //namespace alignment
} //namespace casava
