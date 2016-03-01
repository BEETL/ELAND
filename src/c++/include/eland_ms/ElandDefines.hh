/**
 ** copyright (c) 2007-2010 Illumina, Inc.
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
 ** \file eland_ms/ElandDefines.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_ELAND_DEFINES_H
#define CASAVA_ELAND_MS_ELAND_DEFINES_H

namespace casava
{
namespace eland_ms
{

// uncomment this one to print each oligo being checked
//#define DEBUG_SCAN

#define MAX_HASH_BITS 25

#define NUM_THREADS 1

//#define DONT_SEARCH_REVERSE_STRAND

// if next line is uncommented, only do the first pass out of three
// This means all single error matches are found but only some
// double error matches
//#define ONE_ERROR_PER_OLIGO

// if next line is uncommented, interpret N characters as deletions
// as well as missing bases when doing alignments - deprecated, see documents
//#ifdef ALLOW_N_TO_BE_DELETION

} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_ELAND_DEFINES_H
