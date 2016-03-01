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
 ** @file FragmentResolver.cpp
 **
 ** @brief This class is responsible for either picking the best alignment in
 **        single-end runs or resolving the read fragments in paired-end or
 **        mate-pair runs.
 **
 ** @author Michael Stromberg
 **/

#include "kagu/AlignmentResolver.h"

using namespace std;
namespace cc = casava::common;

namespace casava {
namespace kagu {


// our alignment model to circular model conversion array
const uint8_t AlignmentResolver::mCircularAlignmentModels[8] = { 4, 6, 5, 7, 0, 2, 1, 3 };

// our mate 1 and mate 2 status LUTs
const uint32_t AlignmentResolver::mMate1StatusLUT[6] = { 0, 1,   2,   4,   8,  16 };
const uint32_t AlignmentResolver::mMate2StatusLUT[6] = { 0, 32, 64, 128, 256, 512 };

// constructor
AlignmentResolver::AlignmentResolver(void) {
    ConfigSettings.UseDiscordantFragmentStrategy = false;
}

// destructor
AlignmentResolver::~AlignmentResolver(void) {
    CloseAlignmentReaders();
}

// calculates the percentages associated with a confidence interval equivalent to n standard deviations
double AlignmentResolver::CalculateConfidenceIntervalPercentages(void) {

    const double fragmentLengthConfidenceInterval   = boost::math::erf(ConfigSettings.NumStandardDeviations / boost::math::constants::root_two<double>());
    const double fragmentLengthConfidenceInterval1Z = boost::math::erf(1.0 / boost::math::constants::root_two<double>());

    const double halfNonConfidenceInterval = (1.0 - fragmentLengthConfidenceInterval) / 2.0;
    const double halfNonConfidenceInterval1Z = (1.0 - fragmentLengthConfidenceInterval1Z) / 2.0;

    ConfigSettings.FragmentLengthCILowerPercent   = halfNonConfidenceInterval;
    ConfigSettings.FragmentLengthCIUpperPercent   = 1.0 - halfNonConfidenceInterval;
    ConfigSettings.FragmentLengthCILowerPercent1Z = halfNonConfidenceInterval1Z;
    ConfigSettings.FragmentLengthCIUpperPercent1Z = 1.0 - halfNonConfidenceInterval1Z;

    return fragmentLengthConfidenceInterval;
}

// calculates the fragment length from the reads represented in the iterators
inline uint32_t AlignmentResolver::CalculateFragmentLength(cc::CasavaAlignments::const_iterator& m1It, cc::CasavaAlignments::const_iterator& m2It, ReferenceMetadata& metadata) {

    // check for the proper fragment length
    int32_t m1Begin = m1It->ReferencePosition;
    int32_t m1End   = m1It->ReferencePosition + mStatistics.Mate1ReadLength - 1;
    int32_t m2Begin = m2It->ReferencePosition;
    int32_t m2End   = m2It->ReferencePosition + mStatistics.Mate2ReadLength - 1;

    // fix coordinates that straddle the reference sequence boundaries (circular genomes)
    if(metadata.IsCircular && (m1Begin < 0)) {
        m1Begin += metadata.Length;
        m1End   += metadata.Length;
        metadata.UsedCircularReference = true;
    }

    if(metadata.IsCircular && (m2Begin < 0)) {
        m2Begin += metadata.Length;
        m2End   += metadata.Length;
        metadata.UsedCircularReference = true;
    }

    uint32_t fragmentLength = (m1Begin < m2Begin ? m2End - m1Begin + 1 : m1End - m2Begin + 1);

    // fix mates that straddle the reference sequence boundaries (circular genomes)
    if(metadata.IsCircular && (fragmentLength > ConfigSettings.FragmentLengthThreshold)) {
        const uint32_t fragmentLengthThreshold = metadata.Length - ConfigSettings.FragmentLengthThreshold;
        if(fragmentLength > fragmentLengthThreshold) {
            fragmentLength = (m1Begin < m2Begin ? m1End + metadata.Length - m2Begin + 1 : m2End + metadata.Length - m1Begin + 1);
            metadata.UseCircularAlignmentModel = true;
            metadata.UsedCircularReference     = true;
        }
    }

    return fragmentLength;
}

// calculates the min, median, and max fragment length given two fragment length vectors
void AlignmentResolver::CalculateFragmentLengthStatistics(FragmentLengthStatistics& fls, const FragmentLengthHistogram& hist1, const FragmentLengthHistogram& hist2) {

    // add the fragments from the best alignment models
    vector<uint32_t> fragmentLengths;
    fragmentLengths.insert(fragmentLengths.end(), hist1.begin(), hist1.end());
    fragmentLengths.insert(fragmentLengths.end(), hist2.begin(), hist2.end());

    // identify the min and max points on our confidence interval
    const unsigned int fragmentVectorLength = (uint32_t)fragmentLengths.size();
    sort(fragmentLengths.begin(), fragmentLengths.end());

    fls.Min    = fragmentLengths[(uint32_t)(fragmentVectorLength * ConfigSettings.FragmentLengthCILowerPercent)];
    fls.Median = fragmentLengths[(uint32_t)(fragmentVectorLength * 0.5)];
    fls.Max    = fragmentLengths[(uint32_t)(fragmentVectorLength * ConfigSettings.FragmentLengthCIUpperPercent)];

    fls.LowStdDev  = fls.Median - fragmentLengths[(uint32_t)(fragmentVectorLength * ConfigSettings.FragmentLengthCILowerPercent1Z)];
    fls.HighStdDev = fragmentLengths[(uint32_t)(fragmentVectorLength * ConfigSettings.FragmentLengthCIUpperPercent1Z)] - fls.Median;
}

// closes the input files
void AlignmentResolver::CloseAlignmentReaders(void) {
    if(mMate1Reader.IsOpen()) mMate1Reader.Close();
    if(mMate2Reader.IsOpen()) mMate2Reader.Close();
}

// displays the single-end statistics
void AlignmentResolver::DisplaySingleEndStatistics(SingleEndStatistics& s) {

    const uint32_t totalReads = s.NumContaminants + s.NumFailAQ + s.NumNM
        + s.NumOther + s.NumPassAQ + s.NumQC + s.NumTooManyMatches;

    const double contamPercent = (double)s.NumContaminants   / (double)totalReads * 100.0;
    const double failAQPercent = (double)s.NumFailAQ         / (double)totalReads * 100.0;
    const double nmPercent     = (double)s.NumNM             / (double)totalReads * 100.0;
    const double otherPercent  = (double)s.NumOther          / (double)totalReads * 100.0;
    const double passAQPercent = (double)s.NumPassAQ         / (double)totalReads * 100.0;
    const double qcPercent     = (double)s.NumQC             / (double)totalReads * 100.0;
    const double tmmPercent    = (double)s.NumTooManyMatches / (double)totalReads * 100.0;

    cout << endl;
    cout << setw(20) << left << "RESULTS" << endl;
    cout << "========================================" << endl;

    if(s.NumPassAQ > 0) {
        cout << setw(20) << left << "passed AQ filter" << right << setw(12)
            << s.NumPassAQ << " (" << fixed << setprecision(1)
            << passAQPercent << "%)" << endl;
    }

    if(s.NumFailAQ > 0) {
        cout << setw(20) << left << "failed AQ filter" << right << setw(12)
            << s.NumFailAQ << " (" << fixed << setprecision(1)
            << failAQPercent << "%)" << endl;
    }

    if(s.NumTooManyMatches > 0) {
        cout << setw(20) << left << "too many matches" << right << setw(12)
            << s.NumTooManyMatches << " (" << fixed << setprecision(1)
            << tmmPercent << "%)" << endl;
    }

    if(s.NumNM > 0) {
        cout << setw(20) << left << "not matched" << right << setw(12)
            << s.NumNM << " (" << fixed << setprecision(1)
            << nmPercent << "%)" << endl;
    }

    if(s.NumQC > 0) {
        cout << setw(20) << left << "failed QC filter" << right << setw(12)
            << s.NumQC << " (" << fixed << setprecision(1)
            << qcPercent << "%)" << endl;
    }

    if(s.NumContaminants > 0) {
        cout << setw(20) << left << "contaminant filter" << right << setw(12)
            << s.NumContaminants << " (" << fixed << setprecision(1)
            << contamPercent << "%)" << endl;
    }

    // ideally this should never be displayed
    if(s.NumOther > 0) {
        cout << setw(20) << left << "other" << right << setw(12)
            << s.NumOther << " (" << fixed << setprecision(1)
            << otherPercent << "%)" << endl;
    }

    cout << "----------------------------------------" << endl;

    cout << setw(20) << left << "total" << right << setw(12)
        << totalReads << endl;
}

// returns the alignment with the highest alignment quality
cc::CasavaAlignments::iterator AlignmentResolver::GetBestAlignment(cc::CasavaRead& cr, const double baseLnPcorrect, const double rogCorrection, const AlignmentQuality& aq, const uint32_t seedLength, const uint32_t numAlignments) {

    // initialize
    cc::CasavaAlignments::iterator alIter, bestIter;
    double totalPcorrect  = rogCorrection;
    double bestLnPcorrect = -DBL_MAX;
    const uint32_t totalAlignments = (uint32_t)cr.Alignments.size(); // this is in contrast to numAlignments which might be spliced or genomic

    // find the alignment with the highest ln(Pcorrect)
    if(totalAlignments == 1) {

        bestIter       = cr.Alignments.begin();
        bestLnPcorrect = aq.UpdateLnPcorrect(cr.Qualities, bestIter->MatchDescriptor, baseLnPcorrect);

    } else {

        // store all of our scores
        vector<double> lnPcorrectScores;
        lnPcorrectScores.resize(totalAlignments);
        vector<double>::iterator dIt          = lnPcorrectScores.begin();
        vector<double>::const_iterator bestIt = dIt;

        for(alIter = cr.Alignments.begin(); alIter != cr.Alignments.end(); ++alIter, ++dIt) {

            *dIt = aq.UpdateLnPcorrect(cr.Qualities, alIter->MatchDescriptor, baseLnPcorrect);

            if(*dIt > bestLnPcorrect) {
                bestIter       = alIter;
                bestIt         = dIt;
                bestLnPcorrect = *dIt;
            }
        }

        // add all of the non-best scores
        // N.B. tried simply subtracting the bestLnPcorrect from the total, but rounding errors proved problematic
        for(dIt = lnPcorrectScores.begin(); dIt != lnPcorrectScores.end(); ++dIt) {
            if(dIt != bestIt) totalPcorrect += exp(*dIt);
        }
    }

    // calculate the mate alignment quality for that alignment
    cr.MateAlignmentQuality = aq.CalculateAlignmentQualityFromNeighbors(cr.Qualities, cr.SeedErrors, exp(bestLnPcorrect), totalPcorrect, numAlignments, baseLnPcorrect, seedLength);

    // return the iterator to the best alignment
    return bestIter;
}

// assigns the lower and upper bounds for the desired fragment length confidence interval
void AlignmentResolver::GetFragmentLengthStatistics(FragmentLengthStatistics& fls) {

    const double confidenceInterval = CalculateConfidenceIntervalPercentages();

    // initialize the histograms (we have 8 possible ways of changing mate orientations and ordering)
    AlignmentModelHistograms histograms;
    AlignmentModels models;
    histograms.resize(8);
    models.resize(8);

    for(uint32_t i = 0; i < 8; ++i) {
        models[i].ID = i;
        histograms[i].reserve(500000);
    }

    // defines how often we should check the fragment length distribution
    const uint32_t reportFrequency = 10000;
    uint32_t currentIteration = 0;

    // keep track of some basic statistics
    uint32_t numTestedFragments       = 0;
    uint32_t numTestedUniqueFragments = 0;

    cout << "- phase 1 of 2: building fragment length distribution... ";
    cout.flush();

    Timer flsBenchmark;
    cc::CasavaRead m1, m2;

    string refName;
    ReferenceMetadata metadata;

    bool initializedFragmentLengthStatistics = false;

    while(mMate1Reader.GetNextRead(m1)) {

        // make sure we also get the next mate 2 read
        if(!mMate2Reader.GetNextRead(m2)) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, "Was able to retrieve the next mate 1 read, but unable to retrieve the next mate 2 read."));
        }

        // figure out which mates are unique
        const bool isMate1Unique = (m1.Alignments.size() == 1 ? true : false);
        const bool isMate2Unique = (m2.Alignments.size() == 1 ? true : false);

        if(!m1.Alignments.empty() && !m2.Alignments.empty()) ++numTestedFragments;

        // process if both are unique
        if(isMate1Unique && isMate2Unique) {

            ++numTestedUniqueFragments;

            cc::CasavaAlignments::const_iterator m1It = m1.Alignments.begin();
            cc::CasavaAlignments::const_iterator m2It = m2.Alignments.begin();

            // skip unique mates occurring on different reference sequences
            if(m1It->ReferenceName != m2It->ReferenceName) continue;

            // lookup the reference name if necessary
            if(m1It->ReferenceName != refName) {
                GetReferenceMetadata(m1It->ReferenceName, metadata);
                refName = m1It->ReferenceName;
            }

            // skip mates that have extend beyond the reference sequence
            int32_t m1Begin = m1It->ReferencePosition;
            int32_t m1End   = m1It->ReferencePosition + mStatistics.Mate1ReadLength - 1;
            int32_t m2Begin = m2It->ReferencePosition;
            int32_t m2End   = m2It->ReferencePosition + mStatistics.Mate2ReadLength - 1;

            if((m1Begin < 1) || (m2Begin < 1) || (m1End > (int)metadata.Length) || (m2End > (int)metadata.Length)) continue;

            // calculate the fragment length and alignment model
            metadata.UseCircularAlignmentModel = false;
            metadata.UsedCircularReference     = false;

            const uint32_t fragmentLength = CalculateFragmentLength(m1It, m2It, metadata);
            if(fragmentLength > ConfigSettings.FragmentLengthThreshold) continue;

            const uint8_t alignmentModel = GetAlignmentModel(m1It->ReferencePosition, m1It->IsReverseStrand, m2It->ReferencePosition, m2It->IsReverseStrand, metadata.UseCircularAlignmentModel);

            ++models[alignmentModel].Count;
            histograms[alignmentModel].push_back(fragmentLength);
            ++currentIteration;

            // calculate the fragment length statistics
            if(currentIteration == reportFrequency) {

                AlignmentModels test(models.begin(), models.end());
                sort(test.begin(), test.end());

                const unsigned int activeModelCountSum = test[0].Count + test[1].Count;
                const unsigned int unusedModelCountSum = test[2].Count + test[3].Count + test[4].Count + test[5].Count + test[6].Count + test[7].Count;
                const double consistentPairsPercentage = (double)activeModelCountSum / (double)(unusedModelCountSum + activeModelCountSum);

                // calculate the fragment length distribution if we're under the threshold
                if(consistentPairsPercentage >= ConfigSettings.ConsistentPairsPercent) {
                    ConfigSettings.AlignmentModel1         = test[0].ID;
                    ConfigSettings.AlignmentModel2         = test[1].ID;

                    FragmentLengthStatistics testFLS;
                    CalculateFragmentLengthStatistics(testFLS, 
                        histograms[ConfigSettings.AlignmentModel1],
                        histograms[ConfigSettings.AlignmentModel2]);

                    // stop if we achieve the same results two iterations in a row
                    if(testFLS == fls) break;
                    fls = testFLS;
                    initializedFragmentLengthStatistics = true;
                }

                currentIteration = 0;
            }
        }
    }

    // try to get some fragment length statistics and alignment models even on small datasets
    if (!initializedFragmentLengthStatistics && currentIteration) {

        // TODO: refactor the copy/paste
        AlignmentModels test(models.begin(), models.end());
        sort(test.begin(), test.end());

        const unsigned int activeModelCountSum = test[0].Count + test[1].Count;
        const unsigned int unusedModelCountSum = test[2].Count + test[3].Count + test[4].Count + test[5].Count + test[6].Count + test[7].Count;
        const double consistentPairsPercentage = (double)activeModelCountSum / (double)(unusedModelCountSum + activeModelCountSum);

        // calculate the fragment length distribution if we're under the threshold
        if(consistentPairsPercentage >= ConfigSettings.ConsistentPairsPercent) {
            ConfigSettings.AlignmentModel1 = test[0].ID;
            ConfigSettings.AlignmentModel2 = test[1].ID;

            CalculateFragmentLengthStatistics(fls, 
                histograms[ConfigSettings.AlignmentModel1],
                histograms[ConfigSettings.AlignmentModel2]);

            initializedFragmentLengthStatistics = true;
        }
    }
    cout << "finished (" << fixed << setprecision(1) << flsBenchmark.GetElapsedWallTime() << " s)." << endl << endl;

    // apply the overrides
    if(ConfigSettings.ForceMinFragmentLength) fls.Min = ConfigSettings.MinFragmentLength;
    if(ConfigSettings.ForceMaxFragmentLength) fls.Max = ConfigSettings.MaxFragmentLength;

    // sanity check: did we have any reads?
    if(numTestedFragments == 0) return;

    // display some statistics
    cout << "Fragment length statistics:" << endl;
    cout << "============================" << endl;
    cout << "Confidence interval: " << (confidenceInterval * 100.0) << " %" << endl;
    cout << "Lower bound:         " << fls.Min    << " bp" << endl;
    cout << "Median:              " << fls.Median << " bp" << endl;
    cout << "Upper bound:         " << fls.Max    << " bp" << endl << endl;

    mStatistics.NumFragmentsUsedInFragmentLengthDist = (uint32_t)histograms[ConfigSettings.AlignmentModel1].size()
        + (uint32_t)histograms[ConfigSettings.AlignmentModel2].size();

    // --------------------------------------------------------
    // determine if we should use single-end resolution instead
    // --------------------------------------------------------

    // calculate the fraction of read fragments that are unique vs unique
    const double numUniquePercentage = (double)numTestedUniqueFragments / (double)numTestedFragments;
    if(numUniquePercentage < ConfigSettings.UniquePairPercent) {
        cout << "- the unique read fragment percentage (" << fixed << setprecision(1)
            << numUniquePercentage * 100.0 << " %) was lower than the configured\n  threshold ("
            << ConfigSettings.UniquePairPercent * 100.0
            << " %). The best alignments will be picked separately for each\n  mate sequence." << endl << endl;
        ConfigSettings.UseDiscordantFragmentStrategy = true;
        return;
    }

    // calculate the fraction of unique vs unique read fragments have the chosen alignment model
    uint32_t chosenModelSum = 0;
    uint32_t modelSum       = 0;
    for(uint8_t i = 0; i < 8; ++i) {
        if((i == ConfigSettings.AlignmentModel1) || (i == ConfigSettings.AlignmentModel2)) {
            chosenModelSum += models[i].Count;
        }
        modelSum += models[i].Count;
    }

    const double consistentPairsPercentage = (double)chosenModelSum / (double)modelSum;
    if(consistentPairsPercentage < ConfigSettings.ConsistentPairsPercent) {
        cout << "- the percentage of read fragments with a consistent alignment model (" << fixed << setprecision(1)
            << consistentPairsPercentage * 100.0 << " %)\n  was lower than the configured threshold ("
            << ConfigSettings.ConsistentPairsPercent * 100.0
            << " %). The best alignments will\n  be picked separately for each mate sequence." << endl << endl;
        ConfigSettings.UseDiscordantFragmentStrategy = true;
    }
}

// retrieves the reference sequence metadata
inline void AlignmentResolver::GetReferenceMetadata(const string& referenceName, ReferenceMetadata& metadata) {

    boost::unordered_map<string, ReferenceMetadata>::const_iterator metadataIter = mReferenceMetadataMap.find(referenceName);

    if(metadataIter == mReferenceMetadataMap.end()) {
        if(referenceName.find('/') != string::npos) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("A read was aligned to the reference (%s), but it was not found in the genome sizes XML file. Maybe the --ucn parameter can help.") % referenceName).str()));
        } else {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("A read was aligned to the reference (%s), but it was not found in the genome sizes XML file. Maybe the --ucn parameter can help.") % referenceName).str()));
        }
    }

    metadata.Length     = metadataIter->second.Length;
    metadata.IsCircular = metadataIter->second.IsCircular;
}

// returns the aggregate length of the genome represented in the genome size XML file
uint32_t AlignmentResolver::GetReferenceSequenceLengths(const string& filename) {

    // initialize
    Attributes_t::const_iterator attribCIter;
    string referenceName;
    bool isReferenceCircular;
    uint32_t aggregateLength = 0;
    uint32_t refLength       = 0;

    const string REFERENCE_NAME_TAG = (ConfigSettings.ReferenceRenamingStrategy == USE_CONTIG_NAME ? "contigName" : "fileName");

    // import the genome sizes XML file
    XmlTree xt;
    xt.Import(filename);

    // fill the unordered map
    Entries_t entries;
    xt.GetElements("sequenceSizes.chromosome", entries);

    for(Entries_t::const_iterator entryIter = entries.begin(); entryIter != entries.end(); ++entryIter) {

        bool foundReferenceName = false;
        bool foundTotalBases    = false;
        bool foundIsCircular    = false;
        isReferenceCircular     = false;

        // extract the reference name, the number of bases, and if the reference is circular
        for(attribCIter = entryIter->Attributes.begin(); attribCIter != entryIter->Attributes.end(); ++attribCIter) {
            if(attribCIter->Name == REFERENCE_NAME_TAG) {
                foundReferenceName = true;
                referenceName = attribCIter->Value;
            } else if(attribCIter->Name == "totalBases") {
                foundTotalBases = true;
                try {
                    refLength = boost::lexical_cast<uint32_t>(attribCIter->Value);
                    aggregateLength += refLength;
                } catch(boost::bad_lexical_cast& blc) {
                    BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("Unable to convert a string into an unsigned integer: [%s]") % blc.what()).str()));
                }
            } else if(attribCIter->Name == "isCircular") {
                foundIsCircular = true;
                isReferenceCircular = (attribCIter->Value == "true" ? true : false);
            }
        }

        // add that information to our LUT
        if(!foundReferenceName) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("Unable to find the reference name field (%s) in one of the chromosome entries in the genome sizes XML file.") % REFERENCE_NAME_TAG).str()));
        }

        if(!foundTotalBases) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, "Unable to find the total bases field in one of the chromosome entries in the genome sizes XML file."));
        }

        ReferenceMetadata metadata;
        metadata.Length     = refLength;
        metadata.IsCircular = (foundIsCircular ? isReferenceCircular : false);
        mReferenceMetadataMap[referenceName] = metadata;
    }

    return aggregateLength;
}

// parses the circular references command line option and marks each specified reference as being circular
void AlignmentResolver::MarkCircularReferences(void) {

    // nothing to do
    if(ConfigSettings.CircularReferences.empty()) return;

    // sanity check
    if(mReferenceMetadataMap.empty()) {
        BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("No reference sequences were found in the genome size xml file (%s), but circular reference sequences were specified (%s).") % ConfigSettings.ReferenceSequenceSizeFilename % ConfigSettings.CircularReferences).str()));
    }

    // initialize
    string lowercaseCircularOption = ConfigSettings.CircularReferences;
    std::transform(lowercaseCircularOption.begin(), lowercaseCircularOption.end(), lowercaseCircularOption.begin(), ::tolower);
    boost::unordered_map<std::string, ReferenceMetadata>::iterator refSeqIter;

    // collect the known reference sequence names just in case we need to throw an exception
    ostringstream refSeqList;
    for(refSeqIter = mReferenceMetadataMap.begin(); refSeqIter != mReferenceMetadataMap.end(); ++refSeqIter) {
        refSeqList << "- " << refSeqIter->first << endl;
    }

    if((lowercaseCircularOption == "y") || (lowercaseCircularOption == "yes")) {
        
        // mark all reference sequences circular
        for(refSeqIter = mReferenceMetadataMap.begin(); refSeqIter != mReferenceMetadataMap.end(); ++refSeqIter) {
            refSeqIter->second.IsCircular = true;
        }

    } else {

        // only specific reference sequences should be circular
        vector<string>::const_iterator rsVIter;
        vector<string> refSeqNames;
        boost::split(refSeqNames, ConfigSettings.CircularReferences, boost::is_any_of(","));

        for(rsVIter = refSeqNames.begin(); rsVIter != refSeqNames.end(); ++rsVIter) {

            refSeqIter = mReferenceMetadataMap.find(*rsVIter);
            if(refSeqIter == mReferenceMetadataMap.end()) {

                // the specified reference sequence doesn't exist in our map
                BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("A circular reference sequence was specified (%s), but that reference sequence was not contained in the genome size xml file (%s). The following reference sequence names were parsed from the genome size xml file: \n%s") % *rsVIter % ConfigSettings.ReferenceSequenceSizeFilename % refSeqList.str()).str()));

            } else refSeqIter->second.IsCircular = true;
        }
    }
}

// opens the input files and returns true if the readers contain reads
bool AlignmentResolver::OpenAlignmentReaders(void) {

    const bool resolveFragments = !ConfigSettings.Mate1AlignmentFilename.empty() && !ConfigSettings.Mate2AlignmentFilename.empty();
    bool doesMate1HaveReads = false;
    bool doesMate2HaveReads = false;
    cc::CasavaRead cr;

    // open the mate 1 alignment reader and get the read length
    if(!ConfigSettings.Mate1AlignmentFilename.empty()) {
        mMate1Reader.Open(ConfigSettings.Mate1AlignmentFilename, ConfigSettings.Mate1BaseQualityFilenames,
            ConfigSettings.Mate1TrimmedPrefixBases, ConfigSettings.Mate1TrimmedSuffixBases,
            ConfigSettings.ReferenceRenamingStrategy);
        mStatistics.Mate1ReadLength = mMate1Reader.GetReadLength();

        doesMate1HaveReads = mMate1Reader.GetNextRead(cr);
        mMate1Reader.Rewind();
    }

    // open the mate 2 alignment reader and get the read length
    if(!ConfigSettings.Mate2AlignmentFilename.empty()) {
        mMate2Reader.Open(ConfigSettings.Mate2AlignmentFilename, ConfigSettings.Mate2BaseQualityFilenames,
            ConfigSettings.Mate2TrimmedPrefixBases, ConfigSettings.Mate2TrimmedSuffixBases,
            ConfigSettings.ReferenceRenamingStrategy);
        mStatistics.Mate2ReadLength = mMate2Reader.GetReadLength();

        doesMate2HaveReads = mMate2Reader.GetNextRead(cr);
        mMate2Reader.Rewind();
    }

    // parse the genome size XML file and mark which reference sequences are circular
    mStatistics.GenomeLength = GetReferenceSequenceLengths(ConfigSettings.ReferenceSequenceSizeFilename);
    MarkCircularReferences();

    // return true if the ELAND extended files contain reads
    bool containsReads = false;
    if(resolveFragments && doesMate1HaveReads && doesMate2HaveReads) containsReads = true;
    else if(!resolveFragments && doesMate1HaveReads)                 containsReads = true;

    return containsReads;
}

// resolves the mate-pair or paired-end reads
void AlignmentResolver::ResolveFragments(const FragmentLengthStatistics& fls) {

    // calculate the rest-of-genome correction
    const double rog_mate1 = CalculateRestOfGenomeCorrection(mStatistics.GenomeLength, mStatistics.Mate1ReadLength);
    const double rog_mate2 = CalculateRestOfGenomeCorrection(mStatistics.GenomeLength, mStatistics.Mate2ReadLength);
    const double rog_total = CalculateRestOfGenomeCorrection(mStatistics.GenomeLength, mStatistics.Mate1ReadLength + mStatistics.Mate2ReadLength);

    cout << "- phase 2 of 2: resolving read fragments... ";
    cout.flush();
    Timer resolveBenchmark;

    // open the export writers
    ExportWriter m1Writer, m2Writer;
    m1Writer.Open(ConfigSettings.Mate1ExportFilename);
    m2Writer.Open(ConfigSettings.Mate2ExportFilename);

    // open the anomaly writer
    const bool writeAnomalies = !ConfigSettings.AnomalyFilename.empty();
    AnomalyWriter anomWriter;
    if(writeAnomalies) anomWriter.Open(ConfigSettings.AnomalyFilename);

    // rewind both readers to the beginning
    mMate1Reader.Rewind();
    mMate2Reader.Rewind();
    mMate1Reader.ProvideBaseQualities(true);
    mMate2Reader.ProvideBaseQualities(true);

    // set the discordance flag
    const bool useDiscordantFragmentStrategy = ConfigSettings.UseDiscordantFragmentStrategy;

    // create an instance of our alignment quality calculator
    AlignmentQuality aq;

    string refName;
    ReferenceMetadata metadata;

    cc::CasavaRead m1, m2;
    cc::CasavaAlignments::const_iterator m1Iter, m2Iter, m1ResIter, m2ResIter;

    while(mMate1Reader.GetNextRead(m1)) {

        // make sure we also get the next mate 2 read
        if(!mMate2Reader.GetNextRead(m2)) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, "Was able to retrieve the next mate 1 read, but unable to retrieve the next mate 2 read."));
        }

        // sanity check: make sure both reads are synchronized
        if((m1.XCoord != m2.XCoord) || (m1.YCoord != m2.YCoord)) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("The mate 1 read name (%s) is not the same as the mate 2 read name (%s).")
                % cc::StringUtilities::GetReadName(m1) % cc::StringUtilities::GetReadName(m2)).str()));
        }

        // update read statistics
        ++mStatistics.NumFragments;
        if(!m1.FailedFilters && !m2.FailedFilters) ++mStatistics.NumUniqueFragmentsPassedFiltering;

        // initialize
        double bestMate1LnPcorrect    = -DBL_MAX;
        double bestMate2LnPcorrect    = -DBL_MAX;
        double bestFragmentLnPcorrect = -DBL_MAX;
        double totalFragmentPcorrect  = rog_total;

        bool bestUseCircular = false;

        // figure out which mates are unique or orphaned
        const bool m1Unique = (m1.Alignments.size() == 1 ? true : false);
        const bool m2Unique = (m2.Alignments.size() == 1 ? true : false);

        const bool m1Empty = m1.Alignments.empty();
        const bool m2Empty = m2.Alignments.empty();

        // calculate the base ln(Pcorrect) for both mates with the assumption
        // that every base matches the reference
        const double m1BaseLnPcorrect = aq.GetBaseLnPcorrect(m1.Qualities);
        const double m2BaseLnPcorrect = aq.GetBaseLnPcorrect(m2.Qualities);

        // update our fragment type statistics
        if(!m1Empty && !m2Empty) {
            if(m1Unique && m2Unique)        ++mStatistics.NumUU;
            else if(!m1Unique && !m2Unique) ++mStatistics.NumMM;
            else                            ++mStatistics.NumUM;
        } else if(m1Empty || m2Empty) {
            ++mStatistics.NumOrphans;
        }

        // ===========================
        // handle discordant data sets
        // ===========================

        if(useDiscordantFragmentStrategy) {

            // choose the best alignment for each mate
            bool isMate1Good = false, isMate2Good = false;

            if(!m1Empty) {
                m1ResIter = GetBestAlignment(m1, m1BaseLnPcorrect, rog_mate1, aq, ConfigSettings.Mate1SeedLength, (uint32_t)m1.Alignments.size());
                isMate1Good = m1.MateAlignmentQuality >= ConfigSettings.MinMateAlignmentQuality;
            }

            if(!m2Empty) {
                m2ResIter = GetBestAlignment(m2, m2BaseLnPcorrect, rog_mate2, aq, ConfigSettings.Mate2SeedLength, (uint32_t)m2.Alignments.size());
                isMate2Good = m2.MateAlignmentQuality >= ConfigSettings.MinMateAlignmentQuality;
            }

            // update the orientation statistics
            if(!m1Empty && !m2Empty) {
                const bool isOnSameRef = (m1ResIter->ReferenceName == m2ResIter->ReferenceName);
                const uint8_t currentModel = GetAlignmentModel(m1ResIter->ReferencePosition, m1ResIter->IsReverseStrand, m2ResIter->ReferencePosition, m2ResIter->IsReverseStrand, false);

                if(m1Unique && m2Unique && isOnSameRef && !m1.FailedFilters && !m2.FailedFilters) {
                    ++mStatistics.NumUniqueFragmentsOnSameRefPerAlignmentModel[currentModel];
                }
            }

            // write the results
            if(isMate1Good) {

                if(isMate2Good) {

                    // both mates are good
                    m1Writer.WriteMate(m1, m1ResIter, m2ResIter);
                    m2Writer.WriteMate(m2, m2ResIter, m1ResIter);
                    UpdateReadFragmentStatistics(m1, m2, OS_NoPairedAlignmentDone, SS_BothAlignmentsOK, true);

                } else {

                    // only mate 1 is good (treat as orphaned mate 1)
                    m1Writer.WriteOrphan(m1, m1ResIter);
                    m2Writer.WriteUnaligned(m2);
                    UpdateReadFragmentStatistics(m1, m2, OS_NoPairedAlignmentDone, SS_Read2Poor, false);
                }

            } else if(isMate2Good) {

                // only mate 2 is good (treat as orphaned mate 2)
                m1Writer.WriteUnaligned(m1);
                m2Writer.WriteOrphan(m2, m2ResIter);
                UpdateReadFragmentStatistics(m1, m2, OS_NoPairedAlignmentDone, SS_Read1Poor, false);

            } else {

                // neither of the mates is good
                m1Writer.WriteUnaligned(m1);
                m2Writer.WriteUnaligned(m2);
                UpdateReadFragmentStatistics(m1, m2, OS_BothAlignButNoFeasiblePair, SS_NoPairedAlignmentDone, false);
            }

            if(writeAnomalies) anomWriter.WriteRead(m1, m2, true);
            continue;
        }

        // =====================
        // handle orphaned reads
        // =====================

        if(m1Empty || m2Empty) {

            if(m1Empty && m2Empty) {

                // neither of the mates is good
                m1Writer.WriteUnaligned(m1);
                m2Writer.WriteUnaligned(m2);
                UpdateReadFragmentStatistics(m1, m2, OS_NoMatchToEither, SS_None, false);

            } else if(m1Empty) {

                // orphaned mate 2
                m2ResIter = GetBestAlignment(m2, m2BaseLnPcorrect, rog_mate2, aq, ConfigSettings.Mate2SeedLength, (uint32_t)m2.Alignments.size());
                m1Writer.WriteUnaligned(m1);

                // check if the read meets the minimum mate alignment quality
                if(m2.MateAlignmentQuality < ConfigSettings.MinMateAlignmentQuality) {
                    m2Writer.WriteUnaligned(m2);
                    UpdateReadFragmentStatistics(m1, m2, OS_SingletonRead2, SS_AlignmentPoor, false);
                } else {
                    m2Writer.WriteOrphan(m2, m2ResIter);
                    UpdateReadFragmentStatistics(m1, m2, OS_SingletonRead2, SS_AlignmentOK, false);
                }

            } else {

                // orphaned mate 1
                m1ResIter = GetBestAlignment(m1, m1BaseLnPcorrect, rog_mate1, aq, ConfigSettings.Mate1SeedLength, (uint32_t)m1.Alignments.size());
                m2Writer.WriteUnaligned(m2);

                // check if the read meets the minimum mate alignment quality
                if(m1.MateAlignmentQuality < ConfigSettings.MinMateAlignmentQuality) {
                    m1Writer.WriteUnaligned(m1);
                    UpdateReadFragmentStatistics(m1, m2, OS_SingletonRead1, SS_AlignmentPoor, false);
                } else {
                    m1Writer.WriteOrphan(m1, m1ResIter);
                    UpdateReadFragmentStatistics(m1, m2, OS_SingletonRead1, SS_AlignmentOK, false);
                }
            }

            if(writeAnomalies) anomWriter.WriteRead(m1, m2, true);
            continue;
        }

        // ================================
        // handle concordant read fragments
        // ================================

        // handle all cases where we have both mates
        uint32_t numResolvedFragments = 0;
        for(m1Iter = m1.Alignments.begin(); m1Iter != m1.Alignments.end(); ++m1Iter) {
            for(m2Iter = m2.Alignments.begin(); m2Iter != m2.Alignments.end(); ++m2Iter) {

                // check that the mates appear on the same reference
                if(m1Iter->ReferenceName == m2Iter->ReferenceName) {

                    // lookup the reference name if necessary
                    if(m1Iter->ReferenceName != refName) {
                        GetReferenceMetadata(m1Iter->ReferenceName, metadata);
                        refName   = m1Iter->ReferenceName;
                    }

                    // check for the proper orientation and order
                    metadata.UseCircularAlignmentModel = false;
                    metadata.UsedCircularReference     = false;
                    const uint32_t fragmentLength = CalculateFragmentLength(m1Iter, m2Iter, metadata);
                    const bool isFragmentLengthOK = ((fragmentLength >= fls.Min) && (fragmentLength <= fls.Max) ? true : false);

                    const uint8_t currentModel = GetAlignmentModel(m1Iter->ReferencePosition, m1Iter->IsReverseStrand, m2Iter->ReferencePosition, m2Iter->IsReverseStrand, metadata.UseCircularAlignmentModel);
                    bool isProperlyOrdered = ((currentModel == ConfigSettings.AlignmentModel1) || (currentModel == ConfigSettings.AlignmentModel2) ? true : false);

                    // update fragment statistics for unique vs unique pairs that passed filtering
                    if(m1Unique && m2Unique && !m1.FailedFilters && !m2.FailedFilters) {
                        ++mStatistics.NumUniqueFragmentsOnSameRefPerAlignmentModel[currentModel];
                        if(isProperlyOrdered) {
                            ++mStatistics.NumNominalUniqueFragments;
                            if(fragmentLength < fls.Min)      ++mStatistics.NumNominalSmallFragmentLengths;
                            else if(fragmentLength > fls.Max) ++mStatistics.NumNominalLargeFragmentLengths;
                        }
                    }

                    // set the iterator to the properly resolved pair with the best fragment alignment score
                    if(isProperlyOrdered && isFragmentLengthOK) {

                        const double m1LnPcorrect = aq.UpdateLnPcorrect(m1.Qualities, m1Iter->MatchDescriptor, m1BaseLnPcorrect);
                        const double m2LnPcorrect = aq.UpdateLnPcorrect(m2.Qualities, m2Iter->MatchDescriptor, m2BaseLnPcorrect);
                        const double fragmentLnPcorrect = m1LnPcorrect + m2LnPcorrect;
                        totalFragmentPcorrect += exp(fragmentLnPcorrect);
                        ++numResolvedFragments;

                        if(fragmentLnPcorrect > bestFragmentLnPcorrect) {
                            m1ResIter              = m1Iter;
                            m2ResIter              = m2Iter;
                            bestFragmentLnPcorrect = fragmentLnPcorrect;
                            bestMate1LnPcorrect    = m1LnPcorrect;
                            bestMate2LnPcorrect    = m2LnPcorrect;
                            bestUseCircular        = metadata.UsedCircularReference;
                        }
                    }
                }
            }
        }

        // -------------------------------------------------------------------------
        // handle resolved read fragments (uu, um, mm) and unresolved read fragments
        // -------------------------------------------------------------------------

        // save the resolved fragment
        if(numResolvedFragments >= 1) {

            // calculate the alignment qualities for each mate
            GetBestAlignment(m1, m1BaseLnPcorrect, rog_mate1, aq, ConfigSettings.Mate1SeedLength, (uint32_t)m1.Alignments.size());
            GetBestAlignment(m2, m2BaseLnPcorrect, rog_mate2, aq, ConfigSettings.Mate2SeedLength, (uint32_t)m2.Alignments.size());

            // calculate the fragment alignment quality
            const bool uniqueFragment = numResolvedFragments == 1;
            if(uniqueFragment) {

                uint32_t m1AdjNeighborhood[3];
                uint32_t m2AdjNeighborhood[3];
                aq.AdjustNeighborhood(m1AdjNeighborhood);
                aq.AdjustNeighborhood(m2AdjNeighborhood);

                const uint16_t mate1AQ = aq.CalculateAlignmentQualityFromNeighbors(m1.Qualities, m1AdjNeighborhood, exp(bestMate1LnPcorrect), rog_mate1, numResolvedFragments, m1BaseLnPcorrect, ConfigSettings.Mate1SeedLength);
                const uint16_t mate2AQ = aq.CalculateAlignmentQualityFromNeighbors(m2.Qualities, m2AdjNeighborhood, exp(bestMate2LnPcorrect), rog_mate2, numResolvedFragments, m2BaseLnPcorrect, ConfigSettings.Mate2SeedLength);
                m1.FragmentAlignmentQuality = m2.FragmentAlignmentQuality = mate1AQ + mate2AQ;

            } else {

                const double bestFragmentPcorrect = exp(bestFragmentLnPcorrect);
                m1.FragmentAlignmentQuality = m2.FragmentAlignmentQuality = (uint16_t)floor(-10.0 * log10(1.0 - (bestFragmentPcorrect < totalFragmentPcorrect ? bestFragmentPcorrect / totalFragmentPcorrect : bestFragmentPcorrect)));
            }

            // write all reads that meet the are above the minimum fragment alignment quality
            const bool failedAQ = m1.FragmentAlignmentQuality < ConfigSettings.MinFragmentAlignmentQuality;
            if(failedAQ) {

                ++mStatistics.NumUnresolvedFragments;
                m1Writer.WriteUnaligned(m1);
                m2Writer.WriteUnaligned(m2);
                if(writeAnomalies) anomWriter.WriteRead(m1, m2, true);

            } else {

                if(bestUseCircular) ++mStatistics.NumCircularResolved;

                ++mStatistics.NumResolvedFragments;

                m1Writer.WriteFragment(m1, m1ResIter, m2ResIter);
                m2Writer.WriteFragment(m2, m2ResIter, m1ResIter);
                if(writeAnomalies) anomWriter.WriteRead(m1, m2, false);
            }

            // update statistics
            if(uniqueFragment) UpdateReadFragmentStatistics(m1, m2, OS_UniquePairedAlignment, SS_None, !failedAQ);
            else               UpdateReadFragmentStatistics(m1, m2, OS_ManyPairedAlignments,  SS_None, !failedAQ);

        } else {

            // none of the fragments could be resolved, treat the reads as single-end fragments
            ++mStatistics.NumUnresolvedFragments;

            m1ResIter = GetBestAlignment(m1, m1BaseLnPcorrect, rog_mate1, aq, ConfigSettings.Mate1SeedLength, (uint32_t)m1.Alignments.size());
            m2ResIter = GetBestAlignment(m2, m2BaseLnPcorrect, rog_mate2, aq, ConfigSettings.Mate2SeedLength, (uint32_t)m2.Alignments.size());

            const bool m1HasGoodAQ = (m1.MateAlignmentQuality < ConfigSettings.MinMateAlignmentQuality ? false : true);
            const bool m2HasGoodAQ = (m2.MateAlignmentQuality < ConfigSettings.MinMateAlignmentQuality ? false : true);

            if(m1HasGoodAQ && m2HasGoodAQ) {
                m1Writer.WriteMate(m1, m1ResIter, m2ResIter);
                m2Writer.WriteMate(m2, m2ResIter, m1ResIter);
                UpdateReadFragmentStatistics(m1, m2, OS_BothAlignButNoFeasiblePair, SS_BothAlignmentsOK, true);
            } else if(m1HasGoodAQ  && !m2HasGoodAQ) {
                m1Writer.WriteOrphan(m1, m1ResIter);
                m2Writer.WriteUnaligned(m2);
                UpdateReadFragmentStatistics(m1, m2, OS_BothAlignButNoFeasiblePair, SS_Read2Poor, false);
            } else if(!m1HasGoodAQ && m2HasGoodAQ) {
                m1Writer.WriteUnaligned(m1);
                m2Writer.WriteOrphan(m2, m2ResIter);
                UpdateReadFragmentStatistics(m1, m2, OS_BothAlignButNoFeasiblePair, SS_Read1Poor, false);
            } else if(!m1HasGoodAQ && !m2HasGoodAQ) {
                m1Writer.WriteUnaligned(m1);
                m2Writer.WriteUnaligned(m2);
                UpdateReadFragmentStatistics(m1, m2, OS_BothAlignButNoFeasiblePair, SS_BothAlignButNoFeasiblePair, false);
            }

            if(writeAnomalies) anomWriter.WriteRead(m1, m2, true);
        }
    }

    // print the processing time
    cout << "finished (" << fixed << setprecision(1) << resolveBenchmark.GetElapsedWallTime() << " s)." << endl;

    // close our files
    m1Writer.Close();
    m2Writer.Close();
    if(writeAnomalies) anomWriter.Close();
    mMate1Reader.Close();
    mMate2Reader.Close();

    // print out some summary statistics to the screen
    const uint32_t totalFragments = mStatistics.NumOrphans
        + mStatistics.NumUU
        + mStatistics.NumUM
        + mStatistics.NumMM;

    const uint32_t totalResolvedFragments = mStatistics.NumUUResolved
        + mStatistics.NumUMResolved
        + mStatistics.NumMMResolved;

    const double uuResolvedPercent = (double)mStatistics.NumUUResolved / (double)mStatistics.NumUU * 100.0;
    const double umResolvedPercent = (double)mStatistics.NumUMResolved / (double)mStatistics.NumUM * 100.0;
    const double mmResolvedPercent = (double)mStatistics.NumMMResolved / (double)mStatistics.NumMM * 100.0;
    const double totResolvedPercent = (double)totalResolvedFragments / (double)totalFragments * 100.0;

    cout << endl;
    cout << setw(20) << left << "FRAGMENT ARRANGEMENT" << right << setw(12) << "ORIGINAL" << setw(12) << "RESOLVED" << endl;
    cout << "=====================================================" << endl;

    if(mStatistics.NumOrphans > 0) {
        cout << setw(20) << left << "orphans" << right << setw(12) << mStatistics.NumOrphans << endl;
    }

    if(mStatistics.NumUU > 0) {
        cout << setw(20) << left << "unique vs unique" << right << setw(12)
            << mStatistics.NumUU << setw(12) << mStatistics.NumUUResolved
            << " (" << fixed << setprecision(1) << uuResolvedPercent
            << "%)" << endl;
    }

    if(mStatistics.NumUM > 0) {
        cout << setw(20) << left << "unique vs multiple" << right << setw(12)
            << mStatistics.NumUM  << setw(12) << mStatistics.NumUMResolved
            << " (" << fixed << setprecision(1) << umResolvedPercent
            << "%)" << endl;
    }

    if(mStatistics.NumMM > 0) {
        cout << setw(20) << left << "multiple vs multiple" << right << setw(12)
            << mStatistics.NumMM << setw(12) << mStatistics.NumMMResolved
            << " (" << fixed << setprecision(1) << mmResolvedPercent
            << "%)" << endl;
    }

    cout << "-----------------------------------------------------" << endl;

    cout << setw(20) << left << "total" << right << setw(12) << totalFragments
        << setw(12) << totalResolvedFragments << " (" << fixed
        << setprecision(1) << totResolvedPercent << "%)" << endl;

    if(mStatistics.NumCircularResolved > 0) {
        cout << endl;
        cout << mStatistics.NumCircularResolved << " fragments were resolved using circular reference sequence logic." << endl;
    }
}

// resolves single-end reads
void AlignmentResolver::ResolveMates(void) {

    // calculate the rest-of-genome correction
    const double rog_mate = CalculateRestOfGenomeCorrection(mStatistics.GenomeLength, mStatistics.Mate1ReadLength);

    cout << "- choosing the best alignments... ";

    cout.flush();
    Timer resolveBenchmark;

    // open the export writers
    ExportWriter writer;
    writer.Open(ConfigSettings.Mate1ExportFilename);

    // rewind both readers to the beginning
    mMate1Reader.Rewind();
    mMate1Reader.ProvideBaseQualities(true);

    // create an instance of our alignment quality calculator
    AlignmentQuality aq;

    cc::CasavaRead cr;
    cr.ReadNumber = "1";
    cc::CasavaAlignments::const_iterator bestIter;
    SingleEndStatistics ses;

    while(mMate1Reader.GetNextRead(cr)) {

        // skip all reads with zero alignments
        if(cr.Alignments.empty()) {
            if(cr.IsNm) ++ses.NumNM;
            else if(cr.IsQc) ++ses.NumQC;
            else if(cr.IsTmm) ++ses.NumTooManyMatches;
            else ++ses.NumOther;
            writer.WriteUnaligned(cr);
            continue;
        }

        // calculate the alignment quality
        const double baseLnPcorrect = aq.GetBaseLnPcorrect(cr.Qualities);
        bestIter = GetBestAlignment(cr, baseLnPcorrect, rog_mate, aq, ConfigSettings.Mate1SeedLength, (uint32_t)cr.Alignments.size());

        if(cr.MateAlignmentQuality < ConfigSettings.MinMateAlignmentQuality) {
            writer.WriteUnaligned(cr);
            ++ses.NumFailAQ;
        } else {
            writer.WriteSingleEndRead(cr, bestIter);
            ++ses.NumPassAQ;
        }
    }

    // print the processing time
    cout << "finished (" << fixed << setprecision(1) << resolveBenchmark.GetElapsedWallTime() << " s)." << endl;

    // close our files
    writer.Close();
    mMate1Reader.Close();

    // print out some summary statistics to the screen
    DisplaySingleEndStatistics(ses);
}

// resolves single-end reads (RNA mode)
void AlignmentResolver::ResolveMatesRna(void) {

    // calculate the rest-of-genome correction
    const double rog_mate = CalculateRestOfGenomeCorrection(mStatistics.GenomeLength, mStatistics.Mate1ReadLength);

    cout << "- choosing the best RNA alignments... ";
    cout.flush();
    Timer resolveBenchmark;

    // open the export writers
    ExportWriter writer;
    writer.Open(ConfigSettings.Mate1ExportFilename);

    // open the contamination and splice readers
    cc::ElandExtendedReader contaminationReader;
    contaminationReader.Open(ConfigSettings.ContaminationAlignmentFilename,
        ConfigSettings.Mate1TrimmedPrefixBases,
        ConfigSettings.Mate1TrimmedSuffixBases);
    contaminationReader.ProvideReadName(true);

    cc::ElandExtendedReader spliceReader;
    spliceReader.Open(ConfigSettings.SpliceAlignmentFilename,
        ConfigSettings.Mate1TrimmedPrefixBases,
        ConfigSettings.Mate1TrimmedSuffixBases);
    spliceReader.ProvideReadName(true);

    // rewind both readers to the beginning
    mMate1Reader.Rewind();
    mMate1Reader.ProvideBaseQualities(true);

    // create an instance of our alignment quality calculator
    AlignmentQuality aq;

    cc::CasavaRead cr, contaminationRead, spliceRead;
    cc::CasavaAlignments::const_iterator bestIter;
    cc::CasavaAlignments::iterator alIter;
    SingleEndStatistics ses;

    while(mMate1Reader.GetNextRead(cr)) {

        // retrieve the alignments from the contamination alignments and splice alignments file
        if(!contaminationReader.GetNextRead(contaminationRead)) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, "An alignment was retrieved by the alignment reader, but a matching alignment was not found in the contamination alignment filename."));
        }

        if(!spliceReader.GetNextRead(spliceRead)) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, "An alignment was retrieved by the alignment reader, but a matching alignment was not found in the splice alignment filename."));
        }

        // sanity check: make sure the read names are the same
        if((cr.XCoord != spliceRead.XCoord) || (cr.YCoord != spliceRead.YCoord)) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("The splice site reads seem to be unsynchronized with the genomic reads. Splice read name: [%s], Genomic read name: [%s]") % cc::StringUtilities::GetReadName(spliceRead) % cc::StringUtilities::GetReadName(cr)).str()));
        }

        if((cr.XCoord != contaminationRead.XCoord) || (cr.YCoord != contaminationRead.YCoord)) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("The contamination reads seem to be unsynchronized with the genomic reads. Contamination read name: [%s], Genomic read name: [%s]") % cc::StringUtilities::GetReadName(contaminationRead) % cc::StringUtilities::GetReadName(cr)).str()));
        }

        // handle reads that aligned to contaminants
        if(contaminationRead.IsQc || !contaminationRead.IsNm) {
            cr.Status = (contaminationRead.IsQc ? "QC" : "RM");
            writer.WriteUnaligned(cr);
            ++ses.NumContaminants;
            continue;
        }

        // get the number of genomic alignments
        const uint32_t numGenomicAlignments = (uint32_t)cr.Alignments.size();

        // skip the rest if we have too many genomic alignments
        if(cr.Alignments.empty() && cr.IsTmm) {
            cr.Status = "RM";
            writer.WriteUnaligned(cr);
            ++ses.NumTooManyMatches;
            continue;
        }

        // add the good spliced alignments
        if(!spliceRead.Alignments.empty()) {

            const int32_t numSplicedBases = (int32_t)spliceRead.Bases.size();

            for(alIter = spliceRead.Alignments.begin(); alIter != spliceRead.Alignments.end(); ++alIter) {
                const int32_t spliceLen = cc::StringUtilities::GetSpliceLength(alIter->ContigName);

                if((alIter->ReferencePosition <= spliceLen) &&
                    ((alIter->ReferencePosition + numSplicedBases) > (spliceLen + 1))) {
                        cr.Alignments.push_back(*alIter);
                }
            }
        }

        // if we don't have any alignments at this point skip the rest
        if(cr.Alignments.empty()) {
            if(cr.IsNm)      ++ses.NumNM;
            else if(cr.IsQc) ++ses.NumQC;
            else             ++ses.NumOther;
            writer.WriteUnaligned(cr);
            continue;
        }

        // adjust the seed error neighborhood
        uint32_t numAlignments = numGenomicAlignments;
        if(numAlignments == 0) {
            uninitialized_copy(spliceRead.SeedErrors, spliceRead.SeedErrors + 3, cr.SeedErrors);
            numAlignments = (uint32_t)cr.Alignments.size();
        }

        // find the best alignment
        bestIter = GetBestAlignment(cr, aq.GetBaseLnPcorrect(cr.Qualities), rog_mate, aq, ConfigSettings.Mate1SeedLength, numAlignments);

        // write all reads that meet the are above the minimum mate alignment quality
        const bool failedAQ = cr.MateAlignmentQuality < ConfigSettings.MinMateAlignmentQuality;

        if(failedAQ) {
            ++ses.NumFailAQ;
            cr.Status = "RM";
            writer.WriteUnaligned(cr);
        } else {
            ++ses.NumPassAQ;
            writer.WriteSingleEndRead(cr, bestIter);
        }
    }

    // print the processing time
    cout << "finished (" << fixed << setprecision(1) << resolveBenchmark.GetElapsedWallTime() << " s)." << endl;

    // close our files
    writer.Close();
    mMate1Reader.Close();
    contaminationReader.Close();
    spliceReader.Close();

    // print out some summary statistics to the screen
    DisplaySingleEndStatistics(ses);
}

// sets the use bases for each mate (this should be deprecated)
void AlignmentResolver::SetUseBases(void) {

    string::const_iterator sIter;

    const bool isPairedEnd = !ConfigSettings.Mate1AlignmentFilename.empty() && !ConfigSettings.Mate2AlignmentFilename.empty();
    const bool hasMate1UseBases = !ConfigSettings.Mate1UseBases.empty();
    const bool hasMate2UseBases = !ConfigSettings.Mate2UseBases.empty();

    // handle the default use-bases scenario (skip the last base)
    const bool ignoreLastBase = (isPairedEnd ? (!hasMate1UseBases || !hasMate2UseBases) : !hasMate1UseBases);
    if(ignoreLastBase) {
        cout << "- ignoring the last fastq base in ";
        if(isPairedEnd && !hasMate1UseBases && !hasMate2UseBases) {
            cout << "both mates." << endl << endl;
        } else if(!hasMate1UseBases) {
            cout << "mate 1." << endl << endl;
        } else if(isPairedEnd && !hasMate2UseBases) {
            cout << "mate 2." << endl << endl;
        }
    }

    // derive the trimmed prefix and suffix bases
    if(!hasMate1UseBases) {

        ConfigSettings.Mate1TrimmedPrefixBases = 0;
        ConfigSettings.Mate1TrimmedSuffixBases = 1;

    } else {

        // sanity check: make sure the use bases string contains only Y's or n's
        bool foundError = false;
        string::const_iterator sCIter;
        for(sCIter = ConfigSettings.Mate1UseBases.begin(); sCIter != ConfigSettings.Mate1UseBases.end(); ++sCIter) {
            if((*sCIter != 'Y') && (*sCIter != 'n')) {
                foundError = true;
                break;
            }
        }

        if(foundError) {
            BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("Found an improperly formatted use bases string (%s) which should consist entirely of Y's or n's. Please check your --ub1 parameter.") % ConfigSettings.Mate1UseBases).str()));
        }

        ConfigSettings.Mate1TrimmedPrefixBases = 0;
        ConfigSettings.Mate1TrimmedSuffixBases = 0;

        for(sIter = ConfigSettings.Mate1UseBases.begin(); sIter != ConfigSettings.Mate1UseBases.end(); ++sIter) {
            if(*sIter == 'n') ++ConfigSettings.Mate1TrimmedPrefixBases;
            else break;
        }

        for(sIter = ConfigSettings.Mate1UseBases.end() - 1; sIter != ConfigSettings.Mate1UseBases.end(); --sIter) {
            if(*sIter == 'n') ++ConfigSettings.Mate1TrimmedSuffixBases;
            else break;
        }
    }

    if(isPairedEnd) {

        if(!hasMate2UseBases) {

            ConfigSettings.Mate2TrimmedPrefixBases = 0;
            ConfigSettings.Mate2TrimmedSuffixBases = 1;

        } else {

            // sanity check: make sure the use bases string contains only Y's or n's
            bool foundError = false;
            string::const_iterator sCIter;
            for(sCIter = ConfigSettings.Mate2UseBases.begin(); sCIter != ConfigSettings.Mate2UseBases.end(); ++sCIter) {
                if((*sCIter != 'Y') && (*sCIter != 'n')) {
                    foundError = true;
                    break;
                }
            }

            if(foundError) {
                BOOST_THROW_EXCEPTION(cc::CasavaException(EINVAL, (boost::format("Found an improperly formatted use bases string (%s) which should consist entirely of Y's or n's. Please check your --ub2 parameter.") % ConfigSettings.Mate2UseBases).str()));
            }

            ConfigSettings.Mate2TrimmedPrefixBases = 0;
            ConfigSettings.Mate2TrimmedSuffixBases = 0;

            for(sIter = ConfigSettings.Mate2UseBases.begin(); sIter != ConfigSettings.Mate2UseBases.end(); ++sIter) {
                if(*sIter == 'n') ++ConfigSettings.Mate2TrimmedPrefixBases;
                else break;
            }

            for(sIter = ConfigSettings.Mate2UseBases.end() - 1; sIter != ConfigSettings.Mate2UseBases.end(); --sIter) {
                if(*sIter == 'n') ++ConfigSettings.Mate2TrimmedSuffixBases;
                else break;
            }
        }
    }
}

// updates the read fragment statistics
void AlignmentResolver::UpdateReadFragmentStatistics(cc::CasavaRead& m1, cc::CasavaRead& m2, OutcomeStatus outcomeStatus, SecondaryStatus secondaryStatus, bool updateResolvedStats) {

    // update our resolved fragment statistics
    const bool m1Unique  = (m1.Alignments.size() == 1 ? true : false);
    const bool m2Unique  = (m2.Alignments.size() == 1 ? true : false);

    if(updateResolvedStats) {
        if(m1Unique && m2Unique)        ++mStatistics.NumUUResolved;
        else if(!m1Unique && !m2Unique) ++mStatistics.NumMMResolved;
        else                            ++mStatistics.NumUMResolved;
    }

    // skip reads that failed filters
    if(m1.FailedFilters || m2.FailedFilters) return;

    // derive the hash code
    const uint32_t hashCode = mMate1StatusLUT[m1.MStatus] | mMate2StatusLUT[m2.MStatus] | outcomeStatus | secondaryStatus;

    // updated the counting statistics map
    map<uint32_t, uint32_t>::iterator countIter = mStatistics.Counts.find(hashCode);
    if(countIter == mStatistics.Counts.end()) {
        mStatistics.Counts[hashCode] = 1;
    } else countIter->second++;
}

// writes the statistics into an XML output file
void AlignmentResolver::WriteStatistics(const string& filename, const FragmentLengthStatistics& fls) {

    // populate the XML property tree with fragment length statistics
    XmlTree out;

    // create our use bases string if one doesn't already exist
    if(ConfigSettings.Mate1UseBases.empty()) {
        const string prefix = (ConfigSettings.Mate1TrimmedPrefixBases > 0 ? string(ConfigSettings.Mate1TrimmedPrefixBases, 'n') : "");
        const string suffix = (ConfigSettings.Mate1TrimmedSuffixBases > 0 ? string(ConfigSettings.Mate1TrimmedSuffixBases, 'n') : "");
        ConfigSettings.Mate1UseBases = prefix + string(mStatistics.Mate1ReadLength, 'Y') + suffix;
    }

    if(ConfigSettings.Mate2UseBases.empty()) {
        const string prefix = (ConfigSettings.Mate2TrimmedPrefixBases > 0 ? string(ConfigSettings.Mate2TrimmedPrefixBases, 'n') : "");
        const string suffix = (ConfigSettings.Mate2TrimmedSuffixBases > 0 ? string(ConfigSettings.Mate2TrimmedSuffixBases, 'n') : "");
        ConfigSettings.Mate2UseBases = prefix + string(mStatistics.Mate2ReadLength, 'Y') + suffix;
    }

    // control parameters
    if(!ConfigSettings.CircularReferences.empty()) out.AddStr("ReadPairProperties.ControlParametersUsed.circular", ConfigSettings.CircularReferences);
    out.AddInt("ReadPairProperties.ControlParametersUsed.max-insert-size", (ConfigSettings.ForceMaxFragmentLength ? ConfigSettings.MaxFragmentLength : -1));
    out.AddInt("ReadPairProperties.ControlParametersUsed.min-insert-size", (ConfigSettings.ForceMinFragmentLength ? ConfigSettings.MinFragmentLength : -1));
    out.AddUInt("ReadPairProperties.ControlParametersUsed.min-paired-read-alignment-score", ConfigSettings.MinFragmentAlignmentQuality);
    out.AddDbl("ReadPairProperties.ControlParametersUsed.min-percent-consistent-pairs", ConfigSettings.ConsistentPairsPercent * 100.0);
    out.AddDbl("ReadPairProperties.ControlParametersUsed.min-percent-unique-pairs", ConfigSettings.UniquePairPercent * 100.0);
    out.AddUInt("ReadPairProperties.ControlParametersUsed.min-single-read-alignment-score", ConfigSettings.MinMateAlignmentQuality);
    out.AddDbl("ReadPairProperties.ControlParametersUsed.num-standard-deviations", ConfigSettings.NumStandardDeviations);
    out.AddStr("ReadPairProperties.ControlParametersUsed.use-bases-1", ConfigSettings.Mate1UseBases);
    out.AddStr("ReadPairProperties.ControlParametersUsed.use-bases-2", ConfigSettings.Mate2UseBases);

    // calculate the number of unique vs unique fragments occurring on the same reference
    uint32_t numUniqueFragmentsOnSameRef = 0;

    for(uint8_t i = 0; i < 8; ++i) {
        numUniqueFragmentsOnSameRef += mStatistics.NumUniqueFragmentsOnSameRefPerAlignmentModel[i];
    }

    const bool hasUniqueFragmentsOnSameRef = (numUniqueFragmentsOnSameRef > 0);

    // fragment length
    if(!ConfigSettings.UseDiscordantFragmentStrategy && hasUniqueFragmentsOnSameRef) {
        out.AddUInt("ReadPairProperties.InsertSize.HighSD", fls.HighStdDev);
        out.AddUInt("ReadPairProperties.InsertSize.LowSD", fls.LowStdDev);
        out.AddUInt("ReadPairProperties.InsertSize.Max", fls.Max);
        out.AddUInt("ReadPairProperties.InsertSize.Median", fls.Median);
        out.AddUInt("ReadPairProperties.InsertSize.Min", fls.Min);
    }

    // seed and read length
    out.AddUInt("ReadPairProperties.Length.Read1.SeedLengthForELAND", ConfigSettings.Mate1SeedLength);
    out.AddUInt("ReadPairProperties.Length.Read1.Total", mStatistics.Mate1ReadLength);
    out.AddUInt("ReadPairProperties.Length.Read2.SeedLengthForELAND", ConfigSettings.Mate2SeedLength);
    out.AddUInt("ReadPairProperties.Length.Read2.Total", mStatistics.Mate2ReadLength);

    if(hasUniqueFragmentsOnSameRef) {

        // deduce the nominal orientation in pbp notation
        const uint32_t numFmOrientation = mStatistics.NumUniqueFragmentsOnSameRefPerAlignmentModel[3] + mStatistics.NumUniqueFragmentsOnSameRefPerAlignmentModel[4];
        const uint32_t numFpOrientation = mStatistics.NumUniqueFragmentsOnSameRefPerAlignmentModel[0] + mStatistics.NumUniqueFragmentsOnSameRefPerAlignmentModel[7];
        const uint32_t numRmOrientation = mStatistics.NumUniqueFragmentsOnSameRefPerAlignmentModel[2] + mStatistics.NumUniqueFragmentsOnSameRefPerAlignmentModel[6];
        const uint32_t numRpOrientation = mStatistics.NumUniqueFragmentsOnSameRefPerAlignmentModel[1] + mStatistics.NumUniqueFragmentsOnSameRefPerAlignmentModel[5];

        string nominalOrientation;
        const string orientationConversion[] = { "Fp", "Rp", "Rm", "Fm", "Fm", "Rp", "Rm", "Fp" };

        if(ConfigSettings.UseDiscordantFragmentStrategy) {

            mStatistics.NumNominalUniqueFragments = max(max(numFmOrientation, numFpOrientation), max(numRmOrientation, numRpOrientation));

            if(numFmOrientation == mStatistics.NumNominalUniqueFragments)      nominalOrientation = "Fm";
            else if(numFpOrientation == mStatistics.NumNominalUniqueFragments) nominalOrientation = "Fp";
            else if(numRmOrientation == mStatistics.NumNominalUniqueFragments) nominalOrientation = "Rm";
            else                                                               nominalOrientation = "Rp";

        } else {

            // modified to show just the majority orientation
            nominalOrientation = orientationConversion[ConfigSettings.AlignmentModel1];;
        }

        const double nominalOrientationButLargeInsertPercent = (double)mStatistics.NumNominalLargeFragmentLengths / (double)mStatistics.NumNominalUniqueFragments * 100.0;
        const double nominalOrientationButSmallInsertPercent = (double)mStatistics.NumNominalSmallFragmentLengths / (double)mStatistics.NumNominalUniqueFragments * 100.0;
        const double nominalOrientationPercent               = (double)mStatistics.NumNominalUniqueFragments / (double)numUniqueFragmentsOnSameRef * 100.0;

        if(numFmOrientation > 0)                           out.AddUInt("ReadPairProperties.Orientation.Fm", numFmOrientation);
        if(numFpOrientation > 0)                           out.AddUInt("ReadPairProperties.Orientation.Fp", numFpOrientation);
        if(!nominalOrientation.empty())                    out.AddStr("ReadPairProperties.Orientation.Nominal", nominalOrientation);
        if(mStatistics.NumNominalLargeFragmentLengths > 0) out.AddUInt("ReadPairProperties.Orientation.NominalOrientationButLargeInsert", mStatistics.NumNominalLargeFragmentLengths);
        if(nominalOrientationButLargeInsertPercent > 0.0)  out.AddDbl("ReadPairProperties.Orientation.NominalOrientationButLargeInsertPercent", nominalOrientationButLargeInsertPercent);
        if(mStatistics.NumNominalSmallFragmentLengths > 0) out.AddUInt("ReadPairProperties.Orientation.NominalOrientationButSmallInsert", mStatistics.NumNominalSmallFragmentLengths);
        if(nominalOrientationButSmallInsertPercent > 0.0)  out.AddDbl("ReadPairProperties.Orientation.NominalOrientationButSmallInsertPercent", nominalOrientationButSmallInsertPercent);
        if(nominalOrientationPercent > 0.0)                out.AddDbl("ReadPairProperties.Orientation.NominalOrientationPercent", nominalOrientationPercent);
        if(numRmOrientation > 0)                           out.AddUInt("ReadPairProperties.Orientation.Rm", numRmOrientation);
        if(numRpOrientation > 0)                           out.AddUInt("ReadPairProperties.Orientation.Rp", numRpOrientation);

        // pairs
        out.AddUInt("ReadPairProperties.Pairs.ClustersPassedFiltering", mStatistics.NumUniqueFragmentsPassedFiltering);
        out.AddUInt("ReadPairProperties.Pairs.ClustersTotal", mStatistics.NumFragments);
        out.AddUInt("ReadPairProperties.Pairs.ClustersUsedToComputeInsert", mStatistics.NumFragmentsUsedInFragmentLengthDist);
        out.AddDbl("ReadPairProperties.Pairs.InitialUniquePairsPercent", (double)numUniqueFragmentsOnSameRef / (double)mStatistics.NumFragments * 100.0);

        // reads
        vector<CountingEntry> countingEntries;
        CountingMap::const_iterator cmCIter;
        ostringstream sb;

        for(cmCIter = mStatistics.Counts.begin(); cmCIter != mStatistics.Counts.end(); ++cmCIter) {

            // clear the string builder
            sb.str("");
            sb << "ReadPairProperties.Reads.";
            const uint32_t hashCode = cmCIter->first;

            // test the mate 1 status bits
            if(hashCode & mMate1StatusLUT[cc::MS_ManyAlignmentsFound])       sb << "Read1ManyAlignmentsFound.";
            else if(hashCode & mMate1StatusLUT[cc::MS_NM])                   sb << "Read1NM.";
            else if(hashCode & mMate1StatusLUT[cc::MS_QC])                   sb << "Read1QC.";
            else if(hashCode & mMate1StatusLUT[cc::MS_Repeat])               sb << "Read1Repeat.";
            else if(hashCode & mMate1StatusLUT[cc::MS_SingleAlignmentFound]) sb << "Read1SingleAlignmentFound.";

            // test the mate 2 status bits
            if(hashCode & mMate2StatusLUT[cc::MS_ManyAlignmentsFound])       sb << "Read2ManyAlignmentsFound.";
            else if(hashCode & mMate2StatusLUT[cc::MS_NM])                   sb << "Read2NM.";
            else if(hashCode & mMate2StatusLUT[cc::MS_QC])                   sb << "Read2QC.";
            else if(hashCode & mMate2StatusLUT[cc::MS_Repeat])               sb << "Read2Repeat.";
            else if(hashCode & mMate2StatusLUT[cc::MS_SingleAlignmentFound]) sb << "Read2SingleAlignmentFound.";

            // test the outcome status bits
            if(hashCode & OS_BothAlignButNoFeasiblePair) sb << "BothAlignButNoFeasiblePair";
            else if(hashCode & OS_ManyPairedAlignments)  sb << "ManyPairedAlignments";
            else if(hashCode & OS_NoMatchToEither)       sb << "NoMatchToEither";
            else if(hashCode & OS_NoPairedAlignmentDone) sb << "NoPairedAlignmentDone";
            else if(hashCode & OS_SingletonRead1)        sb << "SingletonRead1";
            else if(hashCode & OS_SingletonRead2)        sb << "SingletonRead2";
            else if(hashCode & OS_UniquePairedAlignment) sb << "UniquePairedAlignment";

            // test the secondary status bits
            if(hashCode & SS_AlignmentOK)                     sb << ".AlignmentOK";
            else if(hashCode & SS_AlignmentPoor)              sb << ".AlignmentPoor";
            else if(hashCode & SS_BothAlignButNoFeasiblePair) sb << ".BothAlignButNoFeasiblePair";
            else if(hashCode & SS_BothAlignmentsOK)           sb << ".BothAlignmentsOK";
            else if(hashCode & SS_NoPairedAlignmentDone)      sb << ".NoPairedAlignmentDone";
            else if(hashCode & SS_Read1Poor)                  sb << ".Read1Poor";
            else if(hashCode & SS_Read2Poor)                  sb << ".Read2Poor";

            // add the counting entry to the XML property tree
            CountingEntry ce;
            ce.Key   = sb.str();
            ce.Value = cmCIter->second;
            countingEntries.push_back(ce);
        }

        sort(countingEntries.begin(), countingEntries.end());
        vector<CountingEntry>::const_iterator ceCIter;
        for(ceCIter = countingEntries.begin(); ceCIter != countingEntries.end(); ++ceCIter) {
            out.AddUInt(ceCIter->Key, ceCIter->Value);
        }
    }

    // write the XML file
    out.Write(filename);
}

}
}
