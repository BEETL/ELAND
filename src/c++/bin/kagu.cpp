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
** @file kagu.cpp
**
** @brief Contains the command-line parsing and prerequisite checking logic
**        for kagu.
**
** @author Michael Stromberg
**/

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <zlib.h>
#include "kagu/AlignmentResolver.h"
#include "kagu/ConfigurationSettings.h"
#include "kagu/Timer.h"

using namespace std;
using namespace casava::common;
using namespace casava::kagu;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

#define DEFAULT_ELAND_EXTENDED_MATE1_FILENAME "reanalysis_1_eland_extended.txt"
#define DEFAULT_ELAND_EXTENDED_MATE2_FILENAME "reanalysis_2_eland_extended.txt"
#define DEFAULT_EXPORT_MATE1_FILENAME         "reanalysis_1_export.txt.gz"
#define DEFAULT_EXPORT_MATE2_FILENAME         "reanalysis_2_export.txt.gz"
#define DEFAULT_REFERENCE_SIZE_FILENAME       "reanalysis_genomesize.xml"
#define DEFAULT_STATISTICS_FILENAME           "reanalysis_pair.xml"

#define DEFAULT_UNIQUE_PAIR_PERCENT        "0.10"
#define DEFAULT_NUM_STANDARD_DEVIATIONS    "3.0"
#define DEFAULT_CONSISTENT_PAIR_PERCENT    "0.70"

// function prototypes
void AppendFilenameExtension(string& filename);
void CreateEmptyGzipFile(string& filename);

int main(int argc, char* argv[]) {

    cout << "---------------------------------------------------------------------------" << endl;
    cout << "Kagu - chooses the best alignments from single-end reads and read fragments" << endl;
    cout << "Casava 1.8.0                                        (C) 2010 Illumina, Inc." << endl;
    cout << "---------------------------------------------------------------------------" << endl;

    // ==============================
    // parse our command line options
    // ==============================

    string minPercentageConsistentAlignmentModel;
    string minPercentageUniqueFragments;
    string numStandardDeviations;

    po::options_description commonOptions("Common options");
    commonOptions.add_options()
        ("ie1", po::value<string>(&ConfigSettings.Mate1AlignmentFilename),
        "the ELAND extended filename for the mate 1 reads")

        ("if1", po::value<Filenames_t>(&ConfigSettings.Mate1BaseQualityFilenames)->multitoken(),
        "the fastq filenames for the mate 1 reads (separated by a space)")

        ("irs", po::value<string>(&ConfigSettings.ReferenceSequenceSizeFilename)->default_value(DEFAULT_REFERENCE_SIZE_FILENAME),
        "the reference size XML filename")

        ("mmaq", po::value<uint16_t>(&ConfigSettings.MinMateAlignmentQuality)->default_value(DEFAULT_MIN_MATE_ALIGNMENT_QUALITY),
        "the alignment quality threshold")

        ("oe1", po::value<string>(&ConfigSettings.Mate1ExportFilename)->default_value(DEFAULT_EXPORT_MATE1_FILENAME),
        "the export filename for the mate 1 reads")

        ("sl1", po::value<uint16_t>(&ConfigSettings.Mate1SeedLength)->default_value(DEFAULT_ELAND_SEED_LENGTH),
        "the ELAND seed length for the mate 1 reads")

        ("ub1", po::value<string>(&ConfigSettings.Mate1UseBases), "specifies which mate 1 bases should be used")

        ("ucn", "use contig names rather than the reference filenames");

    po::options_description pairedEndOptions("Paired-end and mate-pair options");
    pairedEndOptions.add_options()
        ("circular,c", po::value<string>(&ConfigSettings.CircularReferences),
        "instructs the resolver which references are circular. Multiple references can be specified using commas. e.g -c chrM,phix")

        ("flt", po::value<uint32_t>(&ConfigSettings.FragmentLengthThreshold)->default_value(DEFAULT_FRAGMENT_LENGTH_THRESHOLD),
        "fragments longer than this value will be ignored when calculating the fragment length distribution")

        ("ie2", po::value<string>(&ConfigSettings.Mate2AlignmentFilename),
        "the ELAND extended filename for the mate 2 reads")

        ("if2", po::value<Filenames_t>(&ConfigSettings.Mate2BaseQualityFilenames)->multitoken(),
        "the fastq filenames for the mate 2 reads (separated by a space)")

        ("maxfl", po::value<uint32_t>(&ConfigSettings.MaxFragmentLength),
        "the upper bounds of the fragment length")

        ("minfl", po::value<uint32_t>(&ConfigSettings.MinFragmentLength),
        "the lower bounds of the fragment length")

        ("mfaq", po::value<uint16_t>(&ConfigSettings.MinFragmentAlignmentQuality)->default_value(DEFAULT_MIN_FRAGMENT_ALIGNMENT_QUALITY),
        "the fragment alignment quality threshold")

        ("mcf",  po::value<string>(&minPercentageConsistentAlignmentModel)->default_value(DEFAULT_CONSISTENT_PAIR_PERCENT),
        "the minimum percentage of unique fragments that should have the same orientation")

        ("muf", po::value<string>(&minPercentageUniqueFragments)->default_value(DEFAULT_UNIQUE_PAIR_PERCENT),
        "the minimum percentage of fragments that should be unique")

        ("oa",  po::value<string>(&ConfigSettings.AnomalyFilename),
        "the anomaly filename")

        ("oe2", po::value<string>(&ConfigSettings.Mate2ExportFilename)->default_value(DEFAULT_EXPORT_MATE2_FILENAME),
        "the export filename for the mate 2 reads")

        ("os",  po::value<string>(&ConfigSettings.StatisticsFilename)->default_value(DEFAULT_STATISTICS_FILENAME),
        "the statistics XML filename")

        ("sl2", po::value<uint16_t>(&ConfigSettings.Mate2SeedLength)->default_value(DEFAULT_ELAND_SEED_LENGTH),
        "the ELAND seed length for the mate 2 reads")

        ("std", po::value<string>(&numStandardDeviations)->default_value(DEFAULT_NUM_STANDARD_DEVIATIONS),
        "used to calculate the confidence interval in the fragment length distribution")

        ("ub2", po::value<string>(&ConfigSettings.Mate2UseBases), "specifies which mate 2 bases should be used");

    po::options_description rnaFilenameOptions("RNA options");
    rnaFilenameOptions.add_options()
        ("ic", po::value<string>(&ConfigSettings.ContaminationAlignmentFilename),
        "the contamination alignments filename")

        ("is", po::value<string>(&ConfigSettings.SpliceAlignmentFilename),
        "the splice alignments filename");

    po::options_description helpOptions("Help");
    helpOptions.add_options()
        ("help,h", "shows this help text");

    po::options_description all;
    all.add(commonOptions).add(pairedEndOptions).add(rnaFilenameOptions).add(helpOptions);

    // parse the command line options
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, all), vm);
    } catch(po::error& e) {
        cerr << "ERROR: An error occurred while parsing the command line: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    po::notify(vm);

    // display the help menu
    if((argc == 1) || vm.count("help")) {
        all.print(cout);
        exit(EXIT_FAILURE);
    } else cout << endl;

    // =================================
    // check for our required parameters
    // =================================

    bool foundErrors = false;
    ostringstream parsingErrors;

    bool resolveFragments  = vm.count("ie1") && vm.count("ie2");
    const bool isSingleEnd = vm.count("ie1") && !vm.count("ie2");
    const bool useRnaMode  = vm.count("ic")  || vm.count("is");

    // ElandExtendedMate1Filename
    if(!vm.count("ie1")) {
        parsingErrors << "ERROR: An ELAND extended file was not supplied for the mate 1 reads. Please use the --ie1 parameter." << endl << endl;
        foundErrors = true;
    } else {
        if(!fs::exists(ConfigSettings.Mate1AlignmentFilename)) {
            parsingErrors << "ERROR: The ELAND extended file for the mate 1 reads (" << ConfigSettings.Mate1AlignmentFilename << ") could not be found." << endl << endl;
            foundErrors = true;
        }
    }

    // ElandExtendedMate2Filename
    if(!vm.count("ie2")) {
        if(!isSingleEnd) {
            parsingErrors << "ERROR: An ELAND extended file was not supplied for the mate 2 reads. Please use the --ie2 parameter." << endl << endl;
            foundErrors = true;
        }
    } else {
        if(!fs::exists(ConfigSettings.Mate2AlignmentFilename)) {
            parsingErrors << "ERROR: The ELAND extended file for the mate 2 reads (" << ConfigSettings.Mate2AlignmentFilename << ") could not be found." << endl << endl;
            foundErrors = true;
        }
    }

    // ElandExtendedMate1Filename && ElandExtendedMate2Filename
    if(ConfigSettings.Mate1AlignmentFilename == ConfigSettings.Mate2AlignmentFilename) {
        parsingErrors << "ERROR: The ELAND extended filenames are the same for both mates 1 and 2. Please review the --ie1 and --ie2 parameters." << endl << endl;
        foundErrors = true;
    }

    // Mate1UseBases
    if(!vm.count("ub1")) {
        parsingErrors << "ERROR: A use-bases string was not supplied for the mate 1 reads. Please use the --ub1 parameter." << endl << endl;
        foundErrors = true;
    } else {
        replace(ConfigSettings.Mate1UseBases.begin(), ConfigSettings.Mate1UseBases.end(), 'N', 'n');
        replace(ConfigSettings.Mate1UseBases.begin(), ConfigSettings.Mate1UseBases.end(), 'y', 'Y');
        if(ConfigSettings.Mate1UseBases.find_first_not_of("Yn") != string::npos) {
            parsingErrors << "ERROR: Only 'Y' and 'n' characters are allowed in the use-bases string (--ub1). Found: " << ConfigSettings.Mate1UseBases << endl << endl;
            foundErrors = true;
        }
    }

    // Mate2UseBases
    if(!vm.count("ub2") && resolveFragments) {
        parsingErrors << "ERROR: A use-bases string was not supplied for the mate 2 reads. Please use the --ub2 parameter." << endl << endl;
        foundErrors = true;
    } else {
        replace(ConfigSettings.Mate2UseBases.begin(), ConfigSettings.Mate2UseBases.end(), 'N', 'n');
        replace(ConfigSettings.Mate2UseBases.begin(), ConfigSettings.Mate2UseBases.end(), 'y', 'Y');
        if(ConfigSettings.Mate2UseBases.find_first_not_of("Yn") != string::npos) {
            parsingErrors << "ERROR: Only 'Y' and 'n' characters are allowed in the use-bases string (--ub2). Found: " << ConfigSettings.Mate2UseBases << endl << endl;
            foundErrors = true;
        }
    }

    // FastqMate1Filenames
    if(!vm.count("if1")) {
        parsingErrors << "ERROR: A fastq file was not supplied for the mate 1 reads. Please use the --if1 parameter." << endl << endl;
        foundErrors = true;
    } else {
        Filenames_t::const_iterator fIter;
        for(fIter = ConfigSettings.Mate1BaseQualityFilenames.begin(); fIter != ConfigSettings.Mate1BaseQualityFilenames.end(); ++fIter) {
            if(!fs::exists(*fIter)) {
                parsingErrors << "ERROR: A fastq file for the mate 1 reads (" << *fIter << ") could not be found." << endl << endl;
                foundErrors = true;
            }
        }
    }

    // FastqMate2Filenames
    if(!vm.count("if2")) {
        if(!isSingleEnd) {
            parsingErrors << "ERROR: A fastq file was not supplied for the mate 2 reads. Please use the --if2 parameter." << endl << endl;
            foundErrors = true;
        }
    } else {
        Filenames_t::const_iterator fIter;
        for(fIter = ConfigSettings.Mate2BaseQualityFilenames.begin(); fIter != ConfigSettings.Mate2BaseQualityFilenames.end(); ++fIter) {
            if(!fs::exists(*fIter)) {
                parsingErrors << "ERROR: A fastq file for the mate 2 reads (" << *fIter << ") could not be found." << endl << endl;
                foundErrors = true;
            }
        }
    }

    // FastqMate1Filenames & FastqMate2Filenames
    if(vm.count("if1") && vm.count("if2") && (ConfigSettings.Mate1BaseQualityFilenames.size() != ConfigSettings.Mate2BaseQualityFilenames.size())) {
        parsingErrors << "ERROR: A different number of fastq files were supplied for the mate 1 and mate 2 reads. Please check the --iq1 and --iq2 parameters." << endl << endl;
        foundErrors = true;
    }

    // ReferenceSequenceSizesFilename
    if(!vm.count("irs")) {
        parsingErrors << "ERROR: A reference sequence size file was not supplied. Please use the --irs parameter." << endl << endl;
        foundErrors = true;
    } else {
        if(!fs::exists(ConfigSettings.ReferenceSequenceSizeFilename)) {
            parsingErrors << "ERROR: The reference sequence size file (" << ConfigSettings.ReferenceSequenceSizeFilename << ") could not be found." << endl << endl;
            foundErrors = true;
        }
    }

    // OutputExportMate1Filename
    if(!vm.count("oe1")) {
        parsingErrors << "ERROR: A filename was not provided for the export output file for the mate 1 reads. Please use the --oe1 parameter." << endl << endl;
        foundErrors = true;
    } else {
        // add the .gz suffix to the export filename if needed
        AppendFilenameExtension(ConfigSettings.Mate1ExportFilename);
    }

    // OutputExportMate2Filename
    if(!vm.count("oe2")) {
        parsingErrors << "ERROR: A filename was not provided for the export output file for the mate 2 reads. Please use the --oe2 parameter." << endl << endl;
        foundErrors = true;
    } else {
        // add the .gz suffix to the export filename if needed
        AppendFilenameExtension(ConfigSettings.Mate2ExportFilename);
    }

    // OutputExportMate1Filename && OutputExportMate2Filename
    if(ConfigSettings.Mate1ExportFilename == ConfigSettings.Mate2ExportFilename) {
        parsingErrors << "ERROR: The export filenames are the same for both mates 1 and 2. Please review the --oe1 and --oe2 parameters." << endl << endl;
        foundErrors = true;
    }

    // OutputStatisticsFilename
    if(!vm.count("os") && !isSingleEnd) {
        parsingErrors << "ERROR: A filename was not provided for the statistics output file. Please use the --os parameter." << endl << endl;
        foundErrors = true;
    }

    // MinPercentageConsistentAlignmentModel
    ConfigSettings.ConsistentPairsPercent = 0.0;
    try {
        ConfigSettings.ConsistentPairsPercent = boost::lexical_cast<double>(minPercentageConsistentAlignmentModel);
    } catch(boost::bad_lexical_cast& blc) {
        parsingErrors << "ERROR: Unable to convert the minimum percentage of consistent fragments parameter (" << minPercentageConsistentAlignmentModel << ") to a floating point number between 0 and 1: " << blc.what() << endl << endl;
        foundErrors = true;
    }

    if((ConfigSettings.ConsistentPairsPercent < 0.0) || (ConfigSettings.ConsistentPairsPercent > 1.0)) {
        parsingErrors << "ERROR: The minimum percentage of consistent fragments parameter should be a floating point number between 0 and 1." << endl << endl;
        foundErrors = true;
    }

    // MinPercentageUniqueFragments
    ConfigSettings.UniquePairPercent = 0.0;
    try {
        ConfigSettings.UniquePairPercent = boost::lexical_cast<double>(minPercentageUniqueFragments);
    } catch(boost::bad_lexical_cast& blc) {
        parsingErrors << "ERROR: Unable to convert the minimum percentage of unique fragments parameter (" << ConfigSettings.UniquePairPercent << ") to a floating point number between 0 and 1: " << blc.what() << endl << endl;
        foundErrors = true;
    }

    if((ConfigSettings.UniquePairPercent < 0.0) || (ConfigSettings.UniquePairPercent > 1.0)) {
        parsingErrors << "ERROR: The minimum percentage of unique fragments parameter should be a floating point number between 0 and 1." << endl << endl;
        foundErrors = true;
    }

    // MinPercentageUniqueFragments
    ConfigSettings.NumStandardDeviations = 0.0;
    try {
        ConfigSettings.NumStandardDeviations = boost::lexical_cast<double>(numStandardDeviations);
    } catch(boost::bad_lexical_cast& blc) {
        parsingErrors << "ERROR: Unable to convert the number of standard deviations parameter (" << ConfigSettings.NumStandardDeviations << ") to a floating point number larger than 0.0: " << blc.what() << endl << endl;
        foundErrors = true;
    }

    if(ConfigSettings.NumStandardDeviations <= 0.0) {
        parsingErrors << "ERROR: The number of standard deviations parameter should be a floating point number larger than 0.0." << endl << endl;
        foundErrors = true;
    }

    if(useRnaMode) {

        // ContaminationFilename
        if(!vm.count("ic")) {
            parsingErrors << "ERROR: An contamination file was not supplied, but is required when processing RNA data. Please use the --ic parameter." << endl << endl;
            foundErrors = true;
        } else {
            if(!fs::exists(ConfigSettings.ContaminationAlignmentFilename)) {
                parsingErrors << "ERROR: The ELAND extended contamination file (" << ConfigSettings.ContaminationAlignmentFilename << ") could not be found." << endl << endl;
                foundErrors = true;
            }
        }

        // SpliceFilename
        if(!vm.count("is")) {
            parsingErrors << "ERROR: A splice file was not supplied, but is required when processing RNA data. Please use the --is parameter." << endl << endl;
            foundErrors = true;
        } else {
            if(!fs::exists(ConfigSettings.SpliceAlignmentFilename)) {
                parsingErrors << "ERROR: The ELAND extended splice file (" << ConfigSettings.SpliceAlignmentFilename << ") could not be found." << endl << endl;
                foundErrors = true;
            }
        }
    }

    // dump the errors
    if(foundErrors) {
        cerr << parsingErrors.str();
        exit(EXIT_FAILURE);
    }

    // =======================================
    // configure and run the fragment resolver
    // =======================================

    Timer benchmark;

    try {

        AlignmentResolver ar;

        ConfigSettings.ForceMinFragmentLength    = (vm.count("minfl") ? true : false);
        ConfigSettings.ForceMaxFragmentLength    = (vm.count("maxfl") ? true : false);
        ConfigSettings.ReferenceRenamingStrategy = (vm.count("ucn") ? USE_CONTIG_NAME : USE_REFERENCE_NAME);
        ar.SetUseBases();

        // open our alignment readers
        FragmentLengthStatistics fls;
        const bool containsReads = ar.OpenAlignmentReaders();

        // decide if we should resolve read fragments or pick the best alignments
        if(containsReads) {
            if(resolveFragments) {
                ar.GetFragmentLengthStatistics(fls);
                ar.ResolveFragments(fls);
            } else if(useRnaMode) {
                ar.ResolveMatesRna();
            } else {
                ar.ResolveMates();
            }
        }

        // serialize the statistics into the supplied XML filename
        if(resolveFragments) {
            ar.WriteStatistics(ConfigSettings.StatisticsFilename, fls);
        }

        // close our alignment readers
        ar.CloseAlignmentReaders();

        // display a warning message if no reads were found
        if(!containsReads) {
            cerr << "WARNING: No reads were found in the supplied ELAND extended ";
            if(resolveFragments) cerr << "files: " << ConfigSettings.Mate1AlignmentFilename << " & " << ConfigSettings.Mate2AlignmentFilename << endl;
            else cerr << "file: " << ConfigSettings.Mate1AlignmentFilename << endl;

            // create empty export files
            if(vm.count("ie1")) CreateEmptyGzipFile(ConfigSettings.Mate1ExportFilename);
            if(vm.count("ie2")) CreateEmptyGzipFile(ConfigSettings.Mate2ExportFilename);
        }

    } catch(const ExceptionData& ed) {
        cerr << "ERROR: " << ed.getMessage() << endl
            << ed.getContext() << endl;
        exit(EXIT_FAILURE);
    }

    cout << endl << "Kagu elapsed time: " << benchmark.GetElapsedTime() << endl;

    return EXIT_SUCCESS;
}

// appends a .gz filename extension if missing
void AppendFilenameExtension(string& filename) {

    bool needExtension = false;
    const string::size_type dotPos = filename.rfind('.');

    if(dotPos != string::npos) {
        const string extension = filename.substr(dotPos);
        if(extension != ".gz") needExtension = true;
    } else needExtension = true;

    if(needExtension) filename.append(".gz");
}

// creates an empty gzipped file
void CreateEmptyGzipFile(string& filename) {
    gzFile empty = gzopen(filename.c_str(), "wb1");

    if(!empty) {
        BOOST_THROW_EXCEPTION(IoException(EINVAL, (boost::format("Unable to create an empty gzip file (%s).") % filename).str()));
    }

    gzclose(empty);
}
