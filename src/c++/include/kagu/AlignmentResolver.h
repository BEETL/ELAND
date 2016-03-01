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
** @file FragmentResolver.h
**
** @brief This class is responsible for either picking the best alignment in
**        single-end runs or resolving the read fragments in paired-end or
**        mate-pair runs.
**
** @author Michael Stromberg
**/

#pragma once

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/unordered_map.hpp>
#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include "common/StringUtilities.hh"
#include "kagu/AlignmentQuality.h"
#include "kagu/AlignmentReader.h"
#include "kagu/AnomalyWriter.h"
#include "kagu/ConfigurationSettings.h"
#include "kagu/ExportWriter.h"
#include "kagu/Statistics.h"
#include "kagu/Timer.h"
#include "kagu/XmlTree.h"

#define DUMMY_ALIGNMENT_MODEL   100

namespace casava {
namespace kagu {

typedef std::vector<uint32_t> FragmentLengthHistogram;
typedef std::vector<FragmentLengthHistogram> AlignmentModelHistograms;

// stores our fragment length statistics
struct FragmentLengthStatistics {
    uint32_t HighStdDev;
    uint32_t LowStdDev;
    uint32_t Max;
    uint32_t Median;
    uint32_t Min;

    // constructor
    FragmentLengthStatistics(void)
        : HighStdDev(0)
        , LowStdDev(0)
        , Max(0)
        , Median(0)
        , Min(0)
    {}

    // our equality operator
    bool operator==(const FragmentLengthStatistics& fls) const {
        return (Max == fls.Max) && (Median == fls.Median) && (Min == fls.Min)
            && (LowStdDev == fls.LowStdDev) && (HighStdDev == fls.HighStdDev);
    }
};

struct SingleEndStatistics {
    uint32_t NumContaminants;
    uint32_t NumFailAQ;
    uint32_t NumNM;
    uint32_t NumOther;
    uint32_t NumPassAQ;
    uint32_t NumQC;
    uint32_t NumTooManyMatches;

    // constructor
    SingleEndStatistics(void)
        : NumContaminants(0)
        , NumFailAQ(0)
        , NumNM(0)
        , NumOther(0)
        , NumPassAQ(0)
        , NumQC(0)
        , NumTooManyMatches(0)
    {}
};

// stores our alignment model counts
struct AlignmentModel {
    uint8_t ID;
    uint32_t Count;

    // constructor
    AlignmentModel(void)
        : ID(0)
        , Count(0)
    {}

    // our less-than operator
    bool operator<(const AlignmentModel& am) const {
        return am.Count < Count;
    }
};

// stores our reference metadata
struct ReferenceMetadata {
    bool IsCircular;
    bool UseCircularAlignmentModel;
    bool UsedCircularReference;
    uint32_t Length;

    // constructor
    ReferenceMetadata(void)
        : IsCircular(false)
        , UseCircularAlignmentModel(false)
        , UsedCircularReference(false)
        , Length(0)
    {}
};

// our outcome status
enum OutcomeStatus {
    OS_BothAlignButNoFeasiblePair = 1024,
    OS_ManyPairedAlignments       = 2048,
    OS_NoMatchToEither            = 4096,
    OS_NoPairedAlignmentDone      = 8192,
    OS_SingletonRead1             = 16384,
    OS_SingletonRead2             = 32768,
    OS_UniquePairedAlignment      = 65536
};

#define NUM_OUTCOME_STATUS 7

// our secondary status
enum SecondaryStatus {
    SS_AlignmentOK                = 131072,
    SS_AlignmentPoor              = 262144,
    SS_BothAlignButNoFeasiblePair = 524288,
    SS_BothAlignmentsOK           = 1048576,
    SS_NoPairedAlignmentDone      = 2097152,
    SS_None                       = 4194304,
    SS_Read1Poor                  = 8388608,
    SS_Read2Poor                  = 16777216
};

#define NUM_SECONDARY_STATUS 7

typedef std::vector<AlignmentModel> AlignmentModels;

class AlignmentResolver {
public:
    // constructor
    AlignmentResolver(void);
    // destructor
    ~AlignmentResolver(void);
    // closes the input files
    void CloseAlignmentReaders(void);
    // assigns the lower and upper bounds for the desired fragment length confidence interval
    void GetFragmentLengthStatistics(FragmentLengthStatistics& fls);
    // opens the input files and returns true if the readers contain reads
    bool OpenAlignmentReaders(void);
    // resolves the mate-pair or paired-end reads
    void ResolveFragments(const FragmentLengthStatistics& fls);
    // resolves single-end reads
    void ResolveMates(void);
    // resolves single-end reads (RNA mode)
    void ResolveMatesRna(void);
    // sets the use bases for each mate (this should be deprecated)
    void SetUseBases(void);
    // writes the statistics into an XML output file
    void WriteStatistics(const std::string& filename, const FragmentLengthStatistics& fls);

private:
    // calculates the percentages associated with a confidence interval equivalent to n standard deviations
    double CalculateConfidenceIntervalPercentages(void);
    // calculates the fragment length from the reads represented in the iterators
    inline uint32_t CalculateFragmentLength(casava::common::CasavaAlignments::const_iterator& m1It, casava::common::CasavaAlignments::const_iterator& m2It, ReferenceMetadata& metadata);
    // calculates the min, median, and max fragment length given two fragment length std::vectors
    void CalculateFragmentLengthStatistics(FragmentLengthStatistics& fls, const FragmentLengthHistogram& hist1, const FragmentLengthHistogram& hist2);
    // calculates the rest-of-genome correction
    static inline double CalculateRestOfGenomeCorrection(const uint32_t genomeLen, const uint32_t readLen);
    // displays the single-end statistics
    static void DisplaySingleEndStatistics(SingleEndStatistics& s);
    // returns the alignment model associated with ordering and orientation of two mates
    inline uint8_t GetAlignmentModel(const uint32_t mate1Pos, const bool isMate1ReverseStrand, const uint32_t mate2Pos, const bool isMate2ReverseStrand, const bool useCircularAlignmentModel);
    // returns the alignment with the highest alignment quality
    static casava::common::CasavaAlignments::iterator GetBestAlignment(casava::common::CasavaRead& cr, const double baseLnPcorrect, const double rogCorrection, const AlignmentQuality& aq, const uint32_t seedLength, const uint32_t numAlignments);
    // retrieves the reference sequence metadata
    inline void GetReferenceMetadata(const std::string& referenceName, ReferenceMetadata& metadata);
    // returns the aggregate length of the genome represented in the genome size XML file
    uint32_t GetReferenceSequenceLengths(const std::string& filename);
    // parses the circular references command line option and marks each specified reference as being circular
    void MarkCircularReferences(void);
    // updates the read fragment statistics
    void UpdateReadFragmentStatistics(casava::common::CasavaRead& m1, casava::common::CasavaRead& m2, OutcomeStatus outcomeStatus, SecondaryStatus secondaryStatus, bool updateResolvedStats);
    // updates the alignment model and fragment length statistics. Returns true if the mates are resolved.
    bool UpdateAlignmentModelFragmentLengthStatistics(casava::common::CasavaAlignments::const_iterator& m1It, casava::common::CasavaAlignments::const_iterator& m2It, const FragmentLengthStatistics& fls, const bool m1Unique, const bool m2Unique, const bool m1FailedFilter, const bool m2FailedFilter, bool& usedCircularReference);
    // our alignment parser object
    AlignmentReader mMate1Reader;
    AlignmentReader mMate2Reader;
    // our alignment model to circular model conversion array
    static const uint8_t mCircularAlignmentModels[8];
    // stores the statistics required for the XML output file
    Statistics mStatistics;
    // our reference sequence LUT
    boost::unordered_map<std::string, ReferenceMetadata> mReferenceMetadataMap;
    // our mate 1 and mate 2 status LUTs
    static const uint32_t mMate1StatusLUT[6];
    static const uint32_t mMate2StatusLUT[6];
};

// returns the alignment model associated with ordering and orientation of two mates
inline uint8_t AlignmentResolver::GetAlignmentModel(const uint32_t mate1Pos, const bool isMate1ReverseStrand, const uint32_t mate2Pos, const bool isMate2ReverseStrand, const bool useCircularAlignmentModel) {

    // alignment model   orientation
    // =============================
    //       0            FFp -> Fp
    //       1            FRp -> Rp
    //       2            RFp -> Rm
    //       3            RRp -> Fm
    //       4            FFm -> Fm
    //       5            RFm -> Rp
    //       6            FRm -> Rm
    //       7            RRm -> Fp

    uint8_t alignmentModel = DUMMY_ALIGNMENT_MODEL;

    if(mate1Pos < mate2Pos) {

        if(!isMate1ReverseStrand && !isMate2ReverseStrand) alignmentModel = 0;
        if(!isMate1ReverseStrand && isMate2ReverseStrand)  alignmentModel = 1;
        if(isMate1ReverseStrand  && !isMate2ReverseStrand) alignmentModel = 2;
        if(isMate1ReverseStrand  && isMate2ReverseStrand)  alignmentModel = 3;

    } else {

        if(!isMate2ReverseStrand && !isMate1ReverseStrand) alignmentModel = 4;
        if(!isMate2ReverseStrand && isMate1ReverseStrand)  alignmentModel = 5;
        if(isMate2ReverseStrand  && !isMate1ReverseStrand) alignmentModel = 6;
        if(isMate2ReverseStrand  && isMate1ReverseStrand)  alignmentModel = 7;
    }

    // adjust for mates straddling the reference endpoints
    if(useCircularAlignmentModel) alignmentModel = mCircularAlignmentModels[alignmentModel];

    assert(alignmentModel != DUMMY_ALIGNMENT_MODEL);

    return alignmentModel;
}

// calculates the rest-of-genome correction
inline double AlignmentResolver::CalculateRestOfGenomeCorrection(const uint32_t genomeLen, const uint32_t readLen) {
    return exp(log(2.0) + log((double)genomeLen) - (log(4.0) * (double)readLen));
}

}
}
