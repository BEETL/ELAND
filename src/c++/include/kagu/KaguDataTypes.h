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
 ** @file KaguDataTypes.h
 **
 ** @brief Contains data structures used throughout kagu.
 **
 ** @author Michael Stromberg
 **/

#pragma once

namespace casava {
namespace kagu {

// our renaming strategy constants
enum ReferenceRenamingStrategy_t {
    USE_CONTIG_NAME,
    USE_REFERENCE_NAME,
    USE_BOTH_NAMES
};

}
}
