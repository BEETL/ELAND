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
 ** @file Statistics.h
 **
 ** @brief Stores all of the insane amounts of counting statistics that were
 **        historically used by pickBestPair.pl.
 **
 ** @author Michael Stromberg
 **/

#pragma once

#include <stdint.h>
#include <map>
#include <memory>

namespace casava {
namespace kagu {

typedef std::map<uint32_t, uint32_t> CountingMap;

struct Statistics {

    uint32_t GenomeLength;

    uint32_t NumFragments;
    uint32_t NumUniqueFragmentsPassedFiltering;
    uint32_t NumFragmentsUsedInFragmentLengthDist;
    uint32_t NumOrphans;
    uint32_t NumResolvedFragments;
    uint32_t NumUnresolvedFragments;
    uint32_t NumUniqueFragmentsOnSameRefPerAlignmentModel[8];

    uint32_t NumNominalUniqueFragments;
    uint32_t NumNominalLargeFragmentLengths;
    uint32_t NumNominalSmallFragmentLengths;

    uint32_t NumUU;
    uint32_t NumUM;
    uint32_t NumMM;
    uint32_t NumUUResolved;
    uint32_t NumUMResolved;
    uint32_t NumMMResolved;

    uint32_t NumCircularResolved;

    uint32_t Mate1ReadLength;
    uint32_t Mate2ReadLength;

    CountingMap Counts;

    uint32_t MMM;
    uint32_t MMU;

    uint32_t MUM;
    uint32_t MUU;

    uint32_t UMM;
    uint32_t UMU;

    uint32_t UUM;
    uint32_t UUU;

    // constructor
    Statistics(void)
        : GenomeLength(0)
        , NumFragments(0)
        , NumUniqueFragmentsPassedFiltering(0)
        , NumFragmentsUsedInFragmentLengthDist(0)
        , NumOrphans(0)
        , NumResolvedFragments(0)
        , NumUnresolvedFragments(0)
        , NumNominalUniqueFragments(0)
        , NumNominalLargeFragmentLengths(0)
        , NumNominalSmallFragmentLengths(0)
        , NumUU(0)
        , NumUM(0)
        , NumMM(0)
        , NumUUResolved(0)
        , NumUMResolved(0)
        , NumMMResolved(0)
        , NumCircularResolved(0)
        , Mate1ReadLength(0)
        , Mate2ReadLength(0)
        , MMM(0)
        , MMU(0)
        , MUM(0)
        , MUU(0)
        , UMM(0)
        , UMU(0)
        , UUM(0)
        , UUU(0)
    {
        std::uninitialized_fill(NumUniqueFragmentsOnSameRefPerAlignmentModel, NumUniqueFragmentsOnSameRefPerAlignmentModel + 8, 0);
    }
};

// our xml entry data struct
struct CountingEntry {
    std::string Key;
    uint32_t Value;

    bool operator<(const CountingEntry& ce) const {
        return Key < ce.Key;
    }
};

}
}
