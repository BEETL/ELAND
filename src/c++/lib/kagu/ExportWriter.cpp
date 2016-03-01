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
 ** @file ExportWriter.cpp
 **
 ** @brief This class is responsible for writing the export files.
 **
 ** @author Michael Stromberg
 **/

#include "kagu/ExportWriter.h"

using namespace std;
namespace cc = casava::common;

namespace casava {
namespace kagu {

// constructor
ExportWriter::ExportWriter(void) {}

// destructor
ExportWriter::~ExportWriter(void) {
    if(mIsOpen) Close();
}

// closes the file streams
void ExportWriter::Close() {

    // toggle the writer state
    mIsOpen = false;

    // close our files
    gzclose(mOutStream);
}

// opens the export file for the associated mate
void ExportWriter::Open(const string& filename) {

    // open the export file
    mOutStream = gzopen(filename.c_str(), "wb1");

    if(!mOutStream) {
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, (boost::format("Unable to open the export file (%s) for writing.") % filename).str()));
    }

    // toggle the writer state
    mFilename = filename;
    mIsOpen   = true;
}

// writes the alignment info for the current entry
void ExportWriter::WriteAlignmentInfo(const cc::CasavaRead& cr, cc::CasavaAlignments::const_iterator& alIt) {
    const int numBytesWritten = gzprintf(mOutStream, "%s\t%s\t%d\t%c\t%s\t%u\t",
        alIt->ReferenceName.c_str(),
        alIt->ContigName.c_str(),
        alIt->ReferencePosition,
        (alIt->IsReverseStrand ? 'R' : 'F'),
        alIt->MatchDescriptor.c_str(),
        cr.MateAlignmentQuality);

    if(numBytesWritten == 0) {
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, (boost::format("Unable to write to the export file (%s).") % mFilename).str()));
    }
}

// writes a resolved fragment entry to disk
void ExportWriter::WriteFragment(const cc::CasavaRead& cr, cc::CasavaAlignments::const_iterator& alIt, cc::CasavaAlignments::const_iterator& mateIt) {

    CheckOpen();
    WriteHeader(cr);
    WriteAlignmentInfo(cr, alIt);

    const bool isOnSameRef    = (alIt->ReferenceName == mateIt->ReferenceName ? true : false);
    const bool isOnSameContig = (alIt->ContigName    == mateIt->ContigName    ? true : false);
    const int64_t fragmentOffset = (isOnSameRef ? mateIt->ReferencePosition - alIt->ReferencePosition : mateIt->ReferencePosition);

    const int numBytesWritten = gzprintf(mOutStream, "%u\t%s\t%s\t%lld\t%c\t%c\n",
        cr.FragmentAlignmentQuality,
        (isOnSameRef ? "" : mateIt->ReferenceName.c_str()),
        (isOnSameRef && !isOnSameContig ? mateIt->ContigName.c_str() : ""),
        fragmentOffset,
        (mateIt->IsReverseStrand ? 'R' : 'F'),
        (cr.FailedFilters ? 'N' : 'Y'));

    if(numBytesWritten == 0) {
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, (boost::format("Unable to write to the export file (%s).") % mFilename).str()));
    }
}

// writes the header info for the current entry
void ExportWriter::WriteHeader(const cc::CasavaRead& cr) {
    const int numBytesWritten = gzprintf(mOutStream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
        cr.Machine.c_str(),
        cr.RunNumber.c_str(),
        cr.Lane.c_str(),
        cr.Tile.c_str(),
        cr.XCoord.c_str(),
        cr.YCoord.c_str(),
        cr.Index.c_str(),
        cr.ReadNumber.c_str(),
        cr.Bases.c_str(),
        cr.Qualities.c_str());

    if(numBytesWritten == 0) {
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, (boost::format("Unable to write to the export file (%s).") % mFilename).str()));
    }
}

// writes a mate entry to disk
void ExportWriter::WriteMate(const cc::CasavaRead& cr, cc::CasavaAlignments::const_iterator& alIt, cc::CasavaAlignments::const_iterator& mateIt) {

    CheckOpen();
    WriteHeader(cr);
    WriteAlignmentInfo(cr, alIt);

    const bool isOnSameRef = (alIt->ReferenceName == mateIt->ReferenceName ? true : false);
    const int64_t fragmentOffset = (isOnSameRef ? mateIt->ReferencePosition - alIt->ReferencePosition : mateIt->ReferencePosition);

    int numBytesWritten = 0;
    if(!isOnSameRef) {

        numBytesWritten = gzprintf(mOutStream, "0\t%s\t\t%d\t%c\t%c\n",
            mateIt->ReferenceName.c_str(),
            mateIt->ReferencePosition,
            (mateIt->IsReverseStrand ? 'R' : 'F'),
            (cr.FailedFilters ? 'N' : 'Y'));

    } else {

        numBytesWritten = gzprintf(mOutStream, "0\t\t\t%lld\t%c\t%c\n",
            fragmentOffset,
            (mateIt->IsReverseStrand ? 'R' : 'F'),
            (cr.FailedFilters ? 'N' : 'Y'));
    }

    if(numBytesWritten == 0) {
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, (boost::format("Unable to write to the export file (%s).") % mFilename).str()));
    }
}

// writes an orphan mate entry to disk
void ExportWriter::WriteOrphan(const cc::CasavaRead& cr, cc::CasavaAlignments::const_iterator& alIt) {

    CheckOpen();
    WriteHeader(cr);
    WriteAlignmentInfo(cr, alIt);

    const int numBytesWritten = gzprintf(mOutStream, "0\t\t\t0\tN\t%c\n",
        (cr.FailedFilters ? 'N' : 'Y'));

    if(numBytesWritten == 0) {
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, (boost::format("Unable to write to the export file (%s).") % mFilename).str()));
    }
}

// writes a single end read entry to disk
void ExportWriter::WriteSingleEndRead(const cc::CasavaRead& cr, cc::CasavaAlignments::const_iterator& alIt) {

    CheckOpen();
    WriteHeader(cr);
    WriteAlignmentInfo(cr, alIt);

    const int numBytesWritten = gzprintf(mOutStream, "\t\t\t\t\t%c\n",
        (cr.FailedFilters ? 'N' : 'Y'));

    if(numBytesWritten == 0) {
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, (boost::format("Unable to write to the export file (%s).") % mFilename).str()));
    }
}

// writes an unaligned read entry to disk
void ExportWriter::WriteUnaligned(const cc::CasavaRead& cr) {

    CheckOpen();
    WriteHeader(cr);

    const int numBytesWritten = gzprintf(mOutStream, "%s\t\t\t\t\t\t\t\t\t\t\t%c\n",
        cr.Status.c_str(),
        (cr.FailedFilters ? 'N' : 'Y'));

    if(numBytesWritten == 0) {
        BOOST_THROW_EXCEPTION(cc::IoException(EINVAL, (boost::format("Unable to write to the export file (%s).") % mFilename).str()));
    }
}

}
}
