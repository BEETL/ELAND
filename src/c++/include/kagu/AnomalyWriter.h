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
 ** @file AnomalyWriter.h
 **
 ** @brief This class is responsible for creating the anomaly output files.
 **
 ** @author Michael Stromberg
 **/

#pragma once

#include <boost/exception/all.hpp>
#include <boost/format.hpp>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "common/CasavaRead.hh"
#include "common/Exceptions.hh"

namespace casava {
namespace kagu {

class AnomalyWriter {
public:
    // constructor
    AnomalyWriter(void);
    // destructor
    ~AnomalyWriter(void);
    // closes the file stream
    void Close();
    // opens the export file for the associated mate
    void Open(const std::string& filename);
    // writes an anomalous read entry to disk
    void WriteRead(const casava::common::CasavaRead& mate1, const casava::common::CasavaRead& mate2, const bool isAnomalous);

private:
    // toggles the state of the writer
    bool mIsOpen;
    // our output streams
    std::ofstream mOutStream;
};

}
}
