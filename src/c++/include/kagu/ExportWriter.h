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
 ** @file ExportWriter.h
 **
 ** @brief This class is responsible for writing the export files.
 **
 ** @author Michael Stromberg
 **/

#pragma once

#include <boost/exception/all.hpp>
#include <boost/format.hpp>
#include <cstdlib>
#include <iostream>
#include <string>
#include <zlib.h>
#include "common/CasavaRead.hh"
#include "common/Exceptions.hh"

namespace casava {
namespace kagu {

class ExportWriter {
public:
    // constructor
    ExportWriter(void);
    // destructor
    ~ExportWriter(void);
    // closes the file stream
    void Close();
    // opens the export file for the associated mate
    void Open(const std::string& filename);
    // writes a resolved fragment entry to disk
    void WriteFragment(const casava::common::CasavaRead& cr, casava::common::CasavaAlignments::const_iterator& alIt, casava::common::CasavaAlignments::const_iterator& mateIt);
    // writes a mate entry to disk
    void WriteMate(const casava::common::CasavaRead& cr, casava::common::CasavaAlignments::const_iterator& alIt, casava::common::CasavaAlignments::const_iterator& mateIt);
    // writes an orphan mate entry to disk
    void WriteOrphan(const casava::common::CasavaRead& cr, casava::common::CasavaAlignments::const_iterator& alIt);
    // writes a single end read entry to disk
    void WriteSingleEndRead(const casava::common::CasavaRead& cr, casava::common::CasavaAlignments::const_iterator& alIt);
    // writes an unaligned read entry to disk
    void WriteUnaligned(const casava::common::CasavaRead& cr);

private:
    // toggles the state of the writer
    bool mIsOpen;
    // our output streams
    gzFile mOutStream;
    // our export filename
    std::string mFilename;
    // check that we have opened the output file stream
    inline void CheckOpen(void) const;
    // writes the alignment info for the current entry
    void WriteAlignmentInfo(const casava::common::CasavaRead& cr, casava::common::CasavaAlignments::const_iterator& alIt);
    // writes the header info for the current entry
    void WriteHeader(const casava::common::CasavaRead& cr);
};

// check that we have opened the output file stream
inline void ExportWriter::CheckOpen(void) const {
    if(!mIsOpen) {
        BOOST_THROW_EXCEPTION(casava::common::IoException(EINVAL, "An attempt was made to write to the export file without opening it first."));
    }
}

}
}
