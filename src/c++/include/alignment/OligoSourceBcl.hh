// -*- mode: c++; tab-width: 8; -*-
// GlobalUtilities.h
// Common definitions and generally useful functions go here
/*
PROJECT: IMPALA (Inexact Matching Program ALlowing Ambiguity)
MODULE: OligoSourceBcl.h
AUTHOR: A. J. Cox

 * Copyright (c) 2003-2006 Solexa
 *
 ** This software is covered by the "Illumina Genome Analyzer Software
 ** License Agreement" and the "Illumina Source Code License Agreement",
 ** and certain third party copyright/licenses, and any user of this
 ** source file is bound by the terms therein (see accompanying files
 ** Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
 ** Illumina_Source_Code_License_Agreement.pdf and third party
 ** copyright/license notices).
 */

#ifndef CASAVA_ALIGNMENT_OLIGO_SOURCE_BCL_HH
#define CASAVA_ALIGNMENT_OLIGO_SOURCE_BCL_HH

#include "GlobalUtilities.hh"

namespace casava
{
namespace alignment
{


/*****************************************************************************/
// Read oligos from the BCL, LOC and filter files for a list of tiles
// Note: the tiles must be from the same lane. The bclDirectoryList is the
// list of "cycle" directories.
class OligoSourceBcl : public OligoSource
{
  public:
    OligoSourceBcl(const std::vector<fs::path> &bclDirectoryList,
                   const std::vector<fs::path> &barcodeDirectoryList,
                   const fs::path &positionsDirectory,
                   const fs::path &filterDirectory,
                   const boost::format &positionsFileNameFormat,
                   const std::string &machineName,
                   const unsigned int runNumber,
                   const unsigned int lane,
                   const std::vector<unsigned int> tileList,
                   const unsigned int readNumber);
    ~OligoSourceBcl() {delete bclReader_; delete barcodeReader_; delete positionsReader_; delete filtersReader_;}

    // Returns reference to next Sequence (supersedes getNextOligo).
    // isValid will be false if there are no sequences left.
    virtual const casava::common::Sequence& getNextSequenceSelect(bool& isValid,
                                                                  const bool isProvideHeaders,
                                                                  const bool isProvideQualities);

    // Returns reference to last Sequence fetched (supersedes getLastOligo).
    // isValid will be false if there are no sequences left.
    virtual const casava::common::Sequence& getLastSequence(bool& isValid) const;

    // Returns pointer to ASCII sequence of next oligo, or null if at end
    virtual const char* getNextOligo( void );

    // Returns pointer to ASCII sequence of last oligo fetched
    virtual const char* getLastOligo( void ) const;

    // Returns pointer to ASCII name of last oligo read
    virtual const char* getLastName( void );

    // Rewind - all streams
    virtual void rewind ( void );

  private:
    const bool format_;
    const std::vector<fs::path> bclDirectoryList_;
    const std::vector<fs::path> barcodeDirectoryList_;
    const fs::path &positionsDirectory_;
    const fs::path &filterDirectory_;
    const boost::format &positionsFileNameFormat_;
    const unsigned int lane_;
    const std::vector<unsigned int> tileList_;
    std::vector<unsigned int>::const_iterator currentTile_;
    unsigned int currentCluster_;
    unsigned int currentClusterInTile_;
    casava::alignment::BclReader *bclReader_;
    casava::alignment::BclReader *barcodeReader_;
    casava::alignment::PositionsReader *positionsReader_;
    casava::alignment::FiltersReader *filtersReader_;
    casava::common::Sequence sequence_;
    std::string sequenceName_;
    void initializeNewTile();
    bool getCluster();
};

} //namespace alignment
} //namespace casava

#endif //CASAVA_ALIGNMENT_OLIGO_SOURCE_BCL_HH
