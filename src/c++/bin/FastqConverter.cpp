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
** @file FastqConverter
**
** @brief
**
** @author Michael Stromberg
**/

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <string>
#include "common/Exceptions.hh"
#include "common/FileConversion.hh"
#include "kagu/Timer.h"

using namespace std;
namespace bi = boost::iostreams;
namespace cc = casava::common;
namespace ck = casava::kagu;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

struct ConfigurationSettings {
    string BarcodeSequence;
    string FlowcellID;
    string InputFilename;
    string InputFilenameType;
    string OutputFilename;
    string ReadNum;
    string RunID;
    int16_t BaseQuality;
};

const string DEFAULT_FLOWCELL_ID   = "0";
const string DEFAULT_RUN_ID        = "0";
const string DEFAULT_READ_NUM      = "1";

// function prototypes
void AppendFilenameExtension(string& filename);

int main(int argc, char* argv[]) {

    // cout << "---------------------------------------------------------------------------" << endl;
    // cout << "FastqConverter - converts sequence files into compressed FASTQ" << endl;
    // cout << "Casava 1.8.0                                        (C) 2010 Illumina, Inc." << endl;
    // cout << "---------------------------------------------------------------------------" << endl;

    // ==============================
    // parse our command line options
    // ==============================

    ConfigurationSettings cs;

    po::options_description requiredOptions("Required");
    requiredOptions.add_options()
        ("in", po::value<string>(&cs.InputFilename),
        "input filename")

        ("out", po::value<string>(&cs.OutputFilename),
        "output compressed FASTQ filename");

    po::options_description optionalOptions("Optional");
    optionalOptions.add_options()
        ("bq", po::value<int16_t>(&cs.BaseQuality),
        "sets the base quality for FASTA files [0-99]")

        ("bs", po::value<string>(&cs.BarcodeSequence),
        "barcode sequence")

        ("fc", po::value<string>(&cs.FlowcellID)->default_value(DEFAULT_FLOWCELL_ID),
        "flowcell ID")

        ("it", po::value<string>(&cs.InputFilenameType),
        "type of input file: fastq, fasta, export or qseq")

        ("read", po::value<string>(&cs.ReadNum)->default_value(DEFAULT_READ_NUM),
        "the read/mate number [1/2] used for FASTA/FASTQ files")

        ("run", po::value<string>(&cs.RunID)->default_value(DEFAULT_RUN_ID),
        "run ID")

        ("no-compression","don't compress fastq output");

    po::options_description helpOptions("Help");
    helpOptions.add_options()
        ("help,h", "shows this help text");

    po::options_description all;
    all.add(requiredOptions).add(optionalOptions).add(helpOptions);

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
        all.print(cerr);
        exit(EXIT_FAILURE);
    } else cerr << endl;

    // ====================
    // check our parameters
    // ====================

    bool foundErrors = false;
    ostringstream parsingErrors;

    // InputFilename
    bool hasGoodInputFilename = false;
    if(!vm.count("in")) {
        parsingErrors << "ERROR: An input filename was not supplied. Please use the --in parameter." << endl << endl;
        foundErrors = true;
    } else {
        if(!fs::exists(cs.InputFilename)) {
            parsingErrors << "ERROR: The input filename (" << cs.InputFilename << ") could not be found." << endl << endl;
            foundErrors = true;
        } else {
            hasGoodInputFilename = true;
        }
    }

    // OutputFilename
    if(!vm.count("out")) {
        parsingErrors << "ERROR: An output compressed FASTQ filename was not supplied. Please use the --out parameter." << endl << endl;
        foundErrors = true;
    }

    // InputFilenameType
    cc::SeqFormat fmt = cc::SeqFormat_UNKNOWN;
    if(vm.count("it")) {

        boost::to_upper(cs.InputFilenameType);
        if(cs.InputFilenameType == "FASTA") {
            fmt = cc::SeqFormat_FASTA;
        } else if(cs.InputFilenameType == "FASTQ") {
            fmt = cc::SeqFormat_FASTQ;
        } else if(cs.InputFilenameType == "EXPORT") {
            fmt = cc::SeqFormat_EXPORT;
        } else if(cs.InputFilenameType == "QSEQ") {
            fmt = cc::SeqFormat_QSEQ;
        } else {
            parsingErrors << "ERROR: An invalid input format (" << cs.InputFilenameType << ") was supplied. Please use one of the following values: 'fasta', 'fastq', 'export', or 'qseq'." << endl << endl;
            foundErrors = true;
        }

    } else if(hasGoodInputFilename) {

        cerr << "- autodetecting input format: ";
        cerr.flush();

        fmt = cc::FileConversion::CheckInputFormat(cs.InputFilename);

        switch(fmt) {
        case cc::SeqFormat_FASTA:
            cerr << "fasta" << endl;
            break;
        case cc::SeqFormat_FASTQ:
            cerr << "fastq" << endl;
            break;
        case cc::SeqFormat_EXPORT:
            cerr << "export" << endl;
            break;
        case cc::SeqFormat_QSEQ:
            cerr << "qseq" << endl;
            break;
        default:
            cerr << "unknown" << endl;
            parsingErrors << "ERROR: The input format could not be autodetected. Please check the input file or use the --it parameter." << endl << endl;
            foundErrors = true;
            break;
        }
    }

    // BaseQuality
    if(vm.count("bq")) {
        if((cs.BaseQuality < 1) || (cs.BaseQuality > 99)) {
            parsingErrors << "ERROR: An invalid FASTA base quality was supplied. Please supply a base quality in the range [1 - 99]." << endl << endl;
            foundErrors = true;
        }
    } else {
        if(fmt == cc::SeqFormat_FASTA) {
            parsingErrors << "ERROR: A FASTA file was supplied but a default base quality was not provided. Please use the --bq parameter." << endl << endl;
            foundErrors = true;
        }
    }

    // read number
    if(vm.count("read")) {
        try {

            uint32_t readNum = boost::lexical_cast<uint32_t>(cs.ReadNum);

            if((readNum == 0) || (readNum > 4)) {
                parsingErrors << "ERROR: Read numbers should be in the range [1, 4]." << endl << endl;
                foundErrors = true;
            }

        } catch(boost::bad_lexical_cast&) {
            parsingErrors << "ERROR: The read number could not be converted to an integer." << endl << endl;
            foundErrors = true;
        }
    }

    const bool isCompressedOutput(0==vm.count("no-compression"));

    // dump the errors
    if(foundErrors) {
        cerr << parsingErrors.str();
        exit(EXIT_FAILURE);
    }

    // =====================================
    // configure and run the FASTQ converter
    // =====================================

    ck::Timer benchmark;

    try {

        // append the .gz filename extension
        if(isCompressedOutput) {
            AppendFilenameExtension(cs.OutputFilename);
        }

        cc::FileConversion fc(cs.BarcodeSequence, cs.FlowcellID, cs.RunID, cs.ReadNum,isCompressedOutput);

        switch(fmt) {
        case cc::SeqFormat_FASTA:
            fc.FastaToFastq(cs.InputFilename, cs.OutputFilename, static_cast<char>(cs.BaseQuality));
            break;
        case cc::SeqFormat_FASTQ:
            fc.FastqToFastq(cs.InputFilename, cs.OutputFilename);
            break;
        case cc::SeqFormat_EXPORT:
            fc.ExportToFastq(cs.InputFilename, cs.OutputFilename);
            break;
        case cc::SeqFormat_QSEQ:
            fc.QseqToFastq(cs.InputFilename, cs.OutputFilename);
            break;
        default:
            cerr << "ERROR: Could not parse unknown input file format." << endl;
        }

    } catch(const cc::ExceptionData& ed) {
        cerr << "ERROR: " << ed.getMessage() << endl
            << ed.getContext() << endl;
        exit(EXIT_FAILURE);
    }

    cerr << endl << "FastqConverter elapsed time: " << benchmark.GetElapsedTime() << endl;

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

    if(needExtension) {
        filename.append(".gz");
        cerr << "- appending .gz to the output filename (" << filename << ")" << endl;
    }
}
