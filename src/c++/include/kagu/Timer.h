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
 ** @file Timer.h
 **
 ** @brief Provides simple elapsed time support. Implemented because
 **        boost::timer does not give actual walltime.
 **
 ** @author Michael Stromberg
 **/

#pragma once

#include <ctime>
#include <iomanip>
#include <sstream>
#include <stdint.h>
#include <string>

#ifdef _MSC_VER
#include <windows.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

#define NUM_MICROSECONDS_IN_SECOND 1000000.0
#define NUM_MICROSECONDS_IN_100NS  10000000.0

namespace casava {
namespace kagu {

class Timer {
public:
    // constructor
    Timer(void);
    // destructor
    ~Timer(void);
    // returns the elapsed CPU time
    double GetElapsedCpuTime(void) const;
    // returns a string containing both the elapsed wall time and CPU time
    std::string GetElapsedTime(void) const;
    // returns the elapsed wall time
    double GetElapsedWallTime(void) const;
    // restarts the internal timer
    void Restart(void);

private:
    // the start time
#ifdef _MSC_VER
    FILETIME mWallStartTime;
    FILETIME mUserStartTime;
    FILETIME mKernelStartTime;
#else
    struct timeval mWallStartTime;
    struct rusage mStartResUsage;
#endif
};

}
}
