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
 ** @file CasavaRead.hh
 **
 ** @brief This data structure represents a single read with all of the
 **        associated alignments.
 **
 ** @author Michael Stromberg
 **/

#pragma once

#include <stdint.h>
#include <string>
#include <vector>

namespace casava {
namespace common {

struct CasavaAlignment {
    std::string ContigName;
    std::string MatchDescriptor;
    std::string ReferenceName;
    int32_t ReferencePosition; // this is required since ELAND outputs negative positions
    bool IsReverseStrand;

    // constructor
    CasavaAlignment(void)
        : ReferencePosition(0)
        , IsReverseStrand(false)
    {}
};

typedef std::vector<CasavaAlignment> CasavaAlignments;

// mate status
enum MateStatus {
    MS_Unknown              = 0,
    MS_ManyAlignmentsFound  = 1,
    MS_NM                   = 2,
    MS_QC                   = 3,
    MS_Repeat               = 4,
    MS_SingleAlignmentFound = 5
};

#define NUM_MATE_STATUS 5

struct CasavaRead {
    std::string Bases;
    std::string ControlID;
    std::string FlowcellID;
    std::string Index;
    std::string Lane;
    std::string Machine;
    std::string Positions;
    std::string Qualities;
    std::string ReadNumber;
    std::string RunNumber;
    std::string Status;
    std::string Tile;
    std::string XCoord;
    std::string YCoord;
    uint32_t SeedErrors[3];
    CasavaAlignments Alignments;
    uint16_t FragmentAlignmentQuality;
    uint16_t MateAlignmentQuality;
    bool FailedFilters;
    bool IsNm;
    bool IsQc;
    bool IsTmm;
    MateStatus MStatus;

    // constructor
    CasavaRead(void)
        : FragmentAlignmentQuality(0)
        , MateAlignmentQuality(0)
        , FailedFilters(false)
        , IsNm(false)
        , IsQc(false)
        , IsTmm(false)
    {}
};

}
}
