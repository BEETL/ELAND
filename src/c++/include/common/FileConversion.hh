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

#pragma once

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/regex.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace casava {
namespace common {

#define FASTQ_BQ_OFFSET 33

    // data structures
    enum SeqFormat {
        SeqFormat_EXPORT,
        SeqFormat_FASTA,
        SeqFormat_FASTQ,
        SeqFormat_QSEQ,
        SeqFormat_UNKNOWN
    };

    enum FastqFormat {
        FastqFormat_CASAVA17_FC_INDEX,
        FastqFormat_CASAVA17_INDEX,
        FastqFormat_CASAVA17_FC,
        FastqFormat_CASAVA17,
        FastqFormat_CASAVA18,
        FastqFormat_EXTERNAL,
        FastqFormat_UNKNOWN
    };

    struct HeaderData_t {
        std::string Machine;
        std::string RunNumber;
        std::string FlowcellID;
        std::string Lane;
        std::string Tile;
        std::string XCoord;
        std::string YCoord;
        std::string ReadNumber;
        std::string IsFiltered;
        std::string ControlID;
        std::string BarcodeSequence;
    };

    class FileConversion {
    public:
        // constructor
        FileConversion(const std::string& barcodeSequence, const std::string& flowcellID, const std::string& runID, const std::string& readNum,
                       const bool isCompressedOutput=true);
        // destructor
        ~FileConversion(void);
        // determines what the input file format might be
        static SeqFormat CheckInputFormat(const std::string& inputFilename);
        // converts a EXPORT file to a FASTQ file
        void ExportToFastq(const std::string& inputFilename, const std::string& outputFilename);
        // converts a FASTA file to a FASTQ file
        void FastaToFastq(const std::string& inputFilename, const std::string& outputFilename, const char bq);
        // converts a FASTQ file to a FASTQ file
        void FastqToFastq(const std::string& inputFilename, const std::string& outputFilename);
        // converts a QSEQ file to a FASTQ file
        void QseqToFastq(const std::string& inputFilename, const std::string& outputFilename);

    private:
        // closes the I/O streams
        void CloseInput(void);
        void CloseOutput(void);
        // extracts the metadata from the FASTA or FASTQ headers
        void ExtractHeaderData(HeaderData_t& data, const std::string& s, FastqFormat& headerStyle) const;
        // returns the number of columns in the specified string
        static uint32_t GetNumColumns(const std::string& s, const char delimiter);
        // returns true if the character is a nucleotide
        static inline bool IsNucleotide(char c);
        // opens the I/O streams
        void OpenInput(const std::string& inputFilename);
        void OpenOutput(const std::string& outputFilename);
        // converts Phred+64 to Phred+33
        static inline char Phred64ToPhred33(char c);
        // replaces a dot with an N
        static inline char RemoveDots(char c);
        // our I/O file stream variables
        bool mIsInputOpen;
        bool mIsOutputOpen;
        std::ifstream mInStream;
        std::ofstream mOutStream;
        boost::iostreams::filtering_istream mInFilterStream;
        boost::iostreams::filtering_ostream mOutFilterStream;
        // our metadata
        std::string mBarcodeSequence;
        std::string mControlID;
        std::string mFlowcellID;
        std::string mReadNum;
        std::string mRunID;
        //
        bool mIsCompressedOutput;
        // our regular expressions
        static const boost::regex mExternalHeaderRegex;
        static const boost::regex mCasava18HeaderRegex;
        static const boost::regex mCasava17HeaderRegex;
        static const boost::regex mCasava17FcHeaderRegex;
        static const boost::regex mCasava17IdxHeaderRegex;
        static const boost::regex mCasava17FcIdxHeaderRegex;
    };

    // returns true if the character is a nucleotide
    inline bool FileConversion::IsNucleotide(char c) {
        c = toupper(c);
        if((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T') || (c == 'N')) return true;
        return false;
    }

    // converts Phred+64 to Phred+33
    inline char FileConversion::Phred64ToPhred33(char c) {
        return c - 31;
    }

    // replaces a dot with an N
    inline char FileConversion::RemoveDots(char c) {
        return (c == '.' ? 'N' : c);
    }

}
}
