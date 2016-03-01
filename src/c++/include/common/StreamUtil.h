// -*- mode: c++; indent-tabs-mode: nil; c-basic-offset: 4; -*-
//
// Copyright 2008 Illumina, Inc
// Author: Chris Saunders
//
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//

/// \file

#ifndef STREAMUTIL_H
#define STREAMUTIL_H

/*****************************************************************************/

#include <iosfwd>
#include <string>
#include <vector>
#include <cstdio>
#include <boost/filesystem.hpp>

#ifdef _MSC_VER
#include <windows.h>
#endif

/*****************************************************************************/

typedef std::vector<std::string> StrVec;
typedef StrVec::iterator StrVecIter;
typedef StrVec::const_iterator StrVecCIter;

/*****************************************************************************/

/// \brief check that an istream is valid and non-empty
///
void
check_nonempty_istream(std::istream& is,
                       const char* stream_label);

/*****************************************************************************/

/// \brief Consume strings checking that they match a specified list in order.

bool eat_speced_strs(std::istream& is, const StrVec& str_vec);

/*****************************************************************************/

/// \brief returns a file handle to a temporary file (cross-platform)

FILE* casava_tmpfile(void);

/*****************************************************************************/

#endif
