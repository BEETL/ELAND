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
 ** @file AlignmentQuality.cpp
 **
 ** @brief  This class handles all of the alignment quality calculations.
 **
 ** @author Michael Stromberg
 **/

#include "kagu/AlignmentQuality.h"

using namespace std;

namespace casava {
namespace kagu {

// our look-up tables
double AlignmentQuality::mLgPCorrect[100];
double AlignmentQuality::mMismatchCorrection[100];

// constructor
AlignmentQuality::AlignmentQuality(void) {

    // initialize the look-up tables
    for(uint8_t i = 0; i < 100; ++i) {
        const double Pcorrect = 1.0 - pow(10.0, (double)i / -10.0);
        mLgPCorrect[i]         = log(Pcorrect);
        mMismatchCorrection[i] = log((1.0 - Pcorrect) / 3.0) - log(Pcorrect);
    }

    // prevent the logarithmic singularity
    mLgPCorrect[0] = log(1.0 - pow(10.0, 1.0 / -10.0));
}

// destructor
AlignmentQuality::~AlignmentQuality(void) {}

// calculates the alignment quality given seeds with 0-2 errors
uint16_t AlignmentQuality::CalculateAlignmentQualityFromNeighbors(const string& qualities, const uint32_t* pSeedErrors, double Pcorrect, double totalPcorrect, uint32_t numAlignments, double baseLnPcorrect, const uint32_t seedLength) const {

    // extract the neighbor counts from the neighborhood string
    uint32_t errorCounts[3];
    uninitialized_copy(pSeedErrors, pSeedErrors + 3, errorCounts);

    // subtract the # of matched fragments from the neighbor counts
    uint32_t currentError = 0;
    while((numAlignments > 0) && (currentError < 3)) {
        while((errorCounts[currentError] == 0) && (currentError < 3)) ++currentError;
        if(currentError < 3) {
            --errorCounts[currentError];
            --numAlignments;
        }
    }

    // use the two worst base qualities in the seed and the seed error counts
    // to normalize the alignment quality
    string qualitiesCopy = qualities.c_str();
    if((errorCounts[0] + errorCounts[1] + errorCounts[2]) > 0) {

        unsigned char* pQualities = (unsigned char*)qualitiesCopy.data();
        sort(pQualities, pQualities + seedLength);

        double worstPcorrect[3] = { 0, mMismatchCorrection[pQualities[0] - PHRED_BQ_OFFSET], mMismatchCorrection[pQualities[1] - PHRED_BQ_OFFSET] };

        for(uint8_t i = 0; i < 3; ++i) {
            baseLnPcorrect += worstPcorrect[i];
            totalPcorrect  += errorCounts[i] * exp(baseLnPcorrect);
        }
    }

    return (uint16_t)floor(-10.0 * log10(totalPcorrect/(totalPcorrect + Pcorrect)));
}

// calculates ln(Pcorrect) assuming that each base matches the reference
double AlignmentQuality::GetBaseLnPcorrect(const string& qualities) const {

    double baseLnPcorrect = 0.0;
    const char* pQualities = qualities.data();
    for(uint32_t i = 0; i < (uint32_t)qualities.size(); ++i, ++pQualities) {
        baseLnPcorrect += mLgPCorrect[*pQualities - PHRED_BQ_OFFSET];
    }

    return baseLnPcorrect;
}

// retrieves a vector of mismatched base qualities
void AlignmentQuality::GetMismatchBaseQualities(const string& qualities, const string& status, vector<uint8_t>& mismatchBaseQualities) const {

    // initialize
    uint32_t currentPos = 0;
    uint32_t num;
    string match;
    string::const_iterator sEnd;
    const char* pQualities = qualities.data();

    // clear the mismatch vector
    mismatchBaseQualities.clear();

    for(string::const_iterator sIter = status.begin(); sIter != status.end(); ++sIter) {

        if(*sIter == '^') {

            // skip over the INDELs
            ++sIter;
            if(isdigit(*sIter)) {
                sEnd = sIter;
                while((sEnd != status.end()) && (*sEnd != '$')) ++sEnd;
                match = string(sIter, sEnd);
                num = (uint32_t)atoi(match.c_str());
                currentPos += num;
                sIter = sEnd;
            } else {
                sEnd = sIter;
                while((sEnd != status.end()) && (*sEnd != '$')) ++sEnd;
                match = string(sIter, sEnd);
                sIter = sEnd;
            }

        } else if(isdigit(*sIter)) {

            // handle digits
            sEnd = sIter;
            while((sEnd != status.end()) && isdigit(*sEnd)) ++sEnd;
            match = string(sIter, sEnd);
            num = (uint32_t)atoi(match.c_str());
            currentPos += num;
            sIter = sEnd - 1;

        } else {

            // handle characters
            sEnd = sIter;
            while((sEnd != status.end()) && isalpha(*sEnd)) ++sEnd;
            match = string(sIter, sEnd);
            for(uint32_t i = 0; i < (uint32_t)match.size(); ++i, ++currentPos) {
                if(match[i] != 'N') {
                    mismatchBaseQualities.push_back(pQualities[currentPos]);
                }
            }
            sIter = sEnd - 1;
        }
    }
}

// updates ln(Pcorrect) to account for mismatched bases
double AlignmentQuality::UpdateLnPcorrect(const string& qualities, const string& status, const double baseLnPcorrect) const {

    double lnPcorrect = baseLnPcorrect;

    // grab all of the mismatched base qualities
    vector<uint8_t> mismatchBaseQualities;
    GetMismatchBaseQualities(qualities, status, mismatchBaseQualities);

    // apply the mismatch correction
    vector<uint8_t>::const_iterator mmbqIt;
    for(mmbqIt = mismatchBaseQualities.begin(); mmbqIt != mismatchBaseQualities.end(); ++mmbqIt) {
        lnPcorrect += mMismatchCorrection[*mmbqIt - PHRED_BQ_OFFSET];
    }

    return lnPcorrect;
}

}
}
