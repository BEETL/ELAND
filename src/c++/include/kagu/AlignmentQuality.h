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
** @file AlignmentQuality.h
**
** @brief This class handles all of the alignment quality calculations.
**
** @author Michael Stromberg
**/

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdint.h>
#include <string>
#include <vector>

#define PHRED_BQ_OFFSET 64

namespace casava {
namespace kagu {

class AlignmentQuality {
public:
    // constructor
    AlignmentQuality(void);
    // destructor
    ~AlignmentQuality(void);
    // creates a new neighborhood std::string containing one entry in the lowest error neighborhood class
    static inline void AdjustNeighborhood(uint32_t* const pSeedErrors);
    // calculates the alignment quality given seeds with 0-2 errors
    uint16_t CalculateAlignmentQualityFromNeighbors(const std::string& qualities, const uint32_t* pSeedErrors, double Pcorrect, double totalPcorrect, uint32_t numAlignments, double baseLnPcorrect, const uint32_t seedLength) const;
    // calculates ln(Pcorrect) assuming that each base matches the reference
    double GetBaseLnPcorrect(const std::string& qualities) const;
    // updates ln(Pcorrect) to account for mismatched bases
    double UpdateLnPcorrect(const std::string& qualities, const std::string& status, const double baseLnPcorrect) const;

private:
    // our look-up tables
    static double mLgPCorrect[100];
    static double mMismatchCorrection[100];
    // retrieves a vector of mismatched base qualities
    void GetMismatchBaseQualities(const std::string& qualities, const std::string& status, std::vector<uint8_t>& mismatchBaseQualities) const;
};

// modifies the neighborhood to contain at most one entry in the lowest error neighborhood class
inline void AlignmentQuality::AdjustNeighborhood(uint32_t* const pSeedErrors) {
    if(pSeedErrors[0] > 0) {
        pSeedErrors[0] = 1;
        pSeedErrors[1] = 0;
        pSeedErrors[2] = 0;
    } else if(pSeedErrors[1] > 0) {
        pSeedErrors[0] = 0;
        pSeedErrors[1] = 1;
        pSeedErrors[2] = 0;
    } else {
        pSeedErrors[0] = 0;
        pSeedErrors[1] = 0;
        pSeedErrors[2] = 1;
    }
}

}
}
