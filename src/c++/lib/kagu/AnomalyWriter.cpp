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
 ** @file AnomalyWriter.cpp
 **
 ** @brief This class is responsible for creating the anomaly output files.
 **
 ** @author Michael Stromberg
 **/

#include "kagu/AnomalyWriter.h"

using namespace std;
namespace cc = casava::common;

namespace casava {
namespace kagu {

// constructor
AnomalyWriter::AnomalyWriter(void) {}

// destructor
AnomalyWriter::~AnomalyWriter(void) {
    if(mIsOpen) Close();
}

// closes the file streams
void AnomalyWriter::Close() {
    mIsOpen = false;
    mOutStream.close();
}

// opens the export file for the associated mate
void AnomalyWriter::Open(const string& filename) {

    // open the mate 1 file
    mOutStream.open(filename.c_str(), ios::binary);

    if(mOutStream.fail()) {
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, (boost::format("Unable to open the anomaly file (%s) for writing.") % filename).str()));
    }

    // toggle the writer state
    mIsOpen = true;
}

// writes an anomalous read entry to disk
void AnomalyWriter::WriteRead(const cc::CasavaRead& mate1, const cc::CasavaRead& mate2, const bool isAnomalous) {

    // write the read name header
    mOutStream << '>' << mate1.Machine << '_' << setfill('0') << setw(4) << mate1.RunNumber << ':'
        << mate1.Lane << ':' << mate1.Tile << ':' << mate1.XCoord << ':'
        << mate1.YCoord << '#' << mate1.Index << '\t';

    if(isAnomalous) {

        // save details about the anomalous read
        mOutStream << mate1.Bases << '\t'
            << mate2.Bases << '\t'
            << mate1.Qualities << '\t'
            << mate2.Qualities << '\t'
            << (mate1.Positions[0] == '-' ? mate1.Status : mate1.Positions) << '\t'
            << (mate2.Positions[0] == '-' ? mate2.Status : mate2.Positions) << endl;

    } else {

        // this is a good read, we just wanted to waste some more disk space
        mOutStream << "OK" << endl;
    }
}

}
}
