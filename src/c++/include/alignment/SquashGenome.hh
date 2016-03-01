/**
 ** Copyright (c) 2007-2009 Illumina, Inc.
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
 ** \file SquashGenome.hh
 **
 ** \brief declarations for squash genome apis
 **
 ** \author Roman Petrovski
 **/

#ifndef CASAVA_ALIGNMENT_SQUASH_GENOME_HH
#define CASAVA_ALIGNMENT_SQUASH_GENOME_HH

namespace casava
{
namespace alignment
{

void outputSizesToXML( const char* dirName, DIR* pDir, int logLevel);
void unsquash( const char* squashName, int logLevel );
void squash( const char* directoryName, const char* fileName, bool validateNames, bool allowManyContigs, int logLevel );
const std::string getContigNameForbiddenCharacters();

} // namespace alignment
} // namespace casava

#endif // #ifndef CASAVA_ALIGNMENT_SQUASH_GENOME_HH
