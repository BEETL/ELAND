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
 ** @file ElandExtendedReader.h
 **
 ** @brief This class is responsible for parsing ELAND extended files.
 **
 ** @author Michael Stromberg
 **/

#pragma once

#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <cstdlib>
#include <iostream>
#include <stdint.h>
#include <string>
#include "common/CasavaRead.hh"
#include "common/LineReader.hh"
#include "kagu/KaguDataTypes.h"

namespace casava {
namespace common {

class ElandExtendedReader : public LineReader {
public:
    // constructor
    ElandExtendedReader(void);
    // destructor
    ~ElandExtendedReader(void);
    // returns true if there is another read available
    bool GetNextRead(CasavaRead& cr);
    // set to true if base qualities should be parsed, false otherwise
    void ProvideReadName(bool b);
    // sets the desired reference renaming strategy (contig, reference, or both)
    void SetReferenceRenamingStrategy(casava::kagu::ReferenceRenamingStrategy_t strategy);

private:
    // populates the supplied casava read with the ELAND extended read name
    static void ExtractReadName(CasavaRead& cr, const char* pBegin, const char* pEnd);
    // regex that captures information from a VMF file
    static const boost::regex mPositionsRegex;
    // flags used for the reference renaming strategy (both are false by default)
    bool mUseContigNames;
    bool mUseReferenceNames;
    // flags used if read names should be provided (false by default)
    bool mProvideReadName;
};

}
}
