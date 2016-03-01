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
 ** @file ConfigurationSettings.h
 **
 ** @brief The settings data structure is populated by the program_options
 **        parser and used directly by the FragmentResolver class.
 **
 ** @author Michael Stromberg
 **/

#pragma once

#include <string>
#include <stdint.h>
#include <vector>
#include "kagu/KaguDataTypes.h"

#define DEFAULT_FRAGMENT_LENGTH_THRESHOLD      10000
#define DEFAULT_MIN_FRAGMENT_ALIGNMENT_QUALITY 4
#define DEFAULT_MIN_MATE_ALIGNMENT_QUALITY     4
#define DEFAULT_ELAND_SEED_LENGTH              32

namespace casava {
namespace kagu {

typedef std::vector<std::string> Filenames_t;

struct ConfigurationSettings_t {
    std::string AnomalyFilename;
    std::string ContaminationAlignmentFilename;
    std::string ReferenceSequenceSizeFilename;
    std::string SpliceAlignmentFilename;
    std::string StatisticsFilename;

    bool ForceMinFragmentLength;
    bool ForceMaxFragmentLength;
    bool UseDiscordantFragmentStrategy;
    ReferenceRenamingStrategy_t ReferenceRenamingStrategy;
    std::string CircularReferences;

    // store the mate 1 related info
    std::string Mate1AlignmentFilename;
    Filenames_t Mate1BaseQualityFilenames;
    std::string Mate1ExportFilename;
    uint16_t Mate1SeedLength;

    // store the mate 2 related info
    std::string Mate2AlignmentFilename;
    Filenames_t Mate2BaseQualityFilenames;
    std::string Mate2ExportFilename;
    uint16_t Mate2SeedLength;

    // store the two major alignment models
    uint8_t AlignmentModel1;
    uint8_t AlignmentModel2;

    // store the confidence interval percentages
    double FragmentLengthCIUpperPercent;
    double FragmentLengthCIUpperPercent1Z;
    double FragmentLengthCILowerPercent;
    double FragmentLengthCILowerPercent1Z;

    // store the use bases info
    std::string Mate1UseBases;
    std::string Mate2UseBases;
    uint32_t Mate1TrimmedPrefixBases;
    uint32_t Mate1TrimmedSuffixBases;
    uint32_t Mate2TrimmedPrefixBases;
    uint32_t Mate2TrimmedSuffixBases;

    // forcing the fragment length distribution
    uint32_t MinFragmentLength;
    uint32_t MaxFragmentLength;
    uint32_t FragmentLengthThreshold;

    // minimum alignment qualities
    uint16_t MinFragmentAlignmentQuality;
    uint16_t MinMateAlignmentQuality;

    // define the number of standard deviation equivalents our
    // confidence interval should use
    double NumStandardDeviations;

    // unique pair percentage
    double ConsistentPairsPercent;
    double UniquePairPercent;
};

extern ConfigurationSettings_t ConfigSettings;

}
}
