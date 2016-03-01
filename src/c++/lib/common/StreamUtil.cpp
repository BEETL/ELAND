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

#include <iostream>
#include <cstdlib>

#include "common/StreamUtil.h"

void
check_nonempty_istream(std::istream& is,
                       const char* stream_label){

  // rudimentary stream check
  if( ! is ) {
      std::cerr << "ERROR: Invalid input stream: " << stream_label << std::endl;
      exit(EXIT_FAILURE);
  }
  is.get();
  if( is.eof() ){
      std::cerr << "ERROR: Input stream is empty: " << stream_label << std::endl;
      exit(EXIT_FAILURE);
  }
  is.unget();
}

/*****************************************************************************/

bool eat_speced_strs(std::istream& is, const StrVec& str_vec)
{
    if (!is) {
        return false;
    }

    std::string found_str;

    for (StrVecCIter str_citer(str_vec.begin());
         str_citer != str_vec.end(); ++str_citer) {
        is >> found_str;

        if (found_str != *str_citer) {
            return false;
        }
    }

    return true;
}

/*****************************************************************************/

FILE* casava_tmpfile(void)
{
#ifdef _MSC_VER

    // Windows doesn't put temporary files in the user's temp directory,
    // instead it dumps them to c:\ when using tmpfile() *lame*. This 
    // method reproduces the expected tmpfile() behavior.

    // find a unique temporary filename
    char* tempNam = NULL;
    std::string tempFilename;

    char prefixBuffer[255];
    sprintf_s(prefixBuffer, 255, "%u_", GetCurrentProcessId());

    do {
        tempNam = _tempnam(NULL, prefixBuffer);
        tempFilename = tempNam;
        if(tempNam) free(tempNam);
        tempFilename.append(".tmp");
    } while(boost::filesystem::exists(tempFilename));

    if (tempFilename.empty())
    {
        std::cerr << "ERROR: Unable to create a temporary filename." << std::endl;
        exit(EXIT_FAILURE);
    }

    FILE* fp = fopen(tempFilename.c_str(), "w+bD");

    if (!fp)
    {
        std::cerr << "ERROR: Unable to open the temporary file (%s) for writing." << std::endl;
        exit(EXIT_FAILURE);
    }

    return fp;

#else

    // all POSIX-compliant operating systems
    return tmpfile();

#endif
}

/*****************************************************************************/
