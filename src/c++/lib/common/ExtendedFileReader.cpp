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
 **/

#include "common/ExtendedFileReader.h"
#include <cstdlib>
#include <cstring>
#include <iostream>

typedef unsigned int uint;

// *** member functions definitions for ExtendedFileReader ***
// These just call their counterparts in  ExtendedFileReaderImp

// Returns pointer to ASCII sequence of next oligo, or null if at end
bool ExtendedFileReader::getNextEntry( void )
{ return p_->getNextEntry(); }


// Rewind - next oligo read will be first in list
void ExtendedFileReader::rewind ( void )
{ return p_->rewind(); }


const char* ExtendedFileReader::getMachine( void ) const
{ return p_->getMachine(); }


const char* ExtendedFileReader::getRead ( void ) const
{ return p_->getRead(); }


int ExtendedFileReader::getMatchCounter ( void ) const
{ return p_->getMatchCounter(); }


const char* ExtendedFileReader::getXYZ ( void ) const
{ return p_->getXYZ(); }


const char* ExtendedFileReader::getMatches ( void ) const
{ return p_->getMatches(); }



ExtendedFileReaderActual::ExtendedFileReaderActual( const char* exportFileName ) :
  pFile_(fopen(exportFileName, "r")),entry_(NumberOfEntries)
{
  if (pFile_==NULL)
  {
    cerr << "Error in ExtendedFileReaderActual: could not open file "
	 << exportFileName << endl;
    exit (1);
  } // ~if
} // ~ctor

ExtendedFileReaderActual::~ExtendedFileReaderActual()
{
  fclose(pFile_);
} // ~dtor



bool ExtendedFileReaderActual::getNextEntry( void )
{
  const char separator('\t');

  if (fgets(buf_,120000, pFile_)==NULL) return false;
  char* p(buf_);
  p--;

  for (uint i(0); i<entry_.size();i++)
  {
    if (p==NULL)
    {
      cerr << "Error - could not parse line: "
	   << buf_ << endl;
      return false;
    } // ~if
    p++;
    entry_[i]=p;
    p=strchr(p,separator);
    if (p!=NULL) *p='\0';
    //    printf ("%d %s\n", i, entry_[i]);
  } // ~for
  // strchr should return NULL on last iter, could check
  return true;
} // ~ExtendedFileReaderActual::getNextEntry( void )

void ExtendedFileReaderActual::rewind ( void )
{
  ::rewind(pFile_);
} // ~ExportFileImp::rewind ( void )

// *** member functions definitions for ExtendedFileReaderActual ***
// These parse the relevant fields in buf_ into the relevant type


const char* ExtendedFileReaderActual::getMachine( void ) const
{ return entry_[Machine]; }

const char* ExtendedFileReaderActual::getRead ( void ) const
{ return entry_[Read]; }


const char* ExtendedFileReaderActual::getXYZ( void ) const
{
  return entry_[MatchCounter];
}


int ExtendedFileReaderActual::getMatchCounter ( void ) const
{

    const char* p(entry_[Matches]);

    // check for NM
    if( entry_[MatchCounter][0]=='N' )
    {
        // in this case also return 255 (pretend it's a repetitve sequence), because
        // then the orphanAligner will kick
        return 255;
    }


    // check for QC and RM
    if( !isdigit(entry_[MatchCounter][0]) )
    {
        return 0;
    }
    else
    {
        if( entry_[Matches][0]=='-' )
            return 255;
    }


    // loop through the matches and count how many hits we have
    int no_hits = 1;

//    cerr << "parsing " << p << " -> yielding : ";

    p=strchr(p,',');
    while (p!=NULL)
    {
        no_hits++;
        p=strchr(p+1,',');
    }
//    cerr << no_hits << endl;

    return no_hits;
}

const char* ExtendedFileReaderActual::getMatches ( void ) const
{ return entry_[Matches]; }

