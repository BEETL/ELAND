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
 **/

#include "kagu/Timer.h"

using namespace std;

namespace casava {
namespace kagu {

// constructor
Timer::Timer(void) {
    Restart();
}

// destructor
Timer::~Timer(void) {}

// returns the elapsed CPU time
double Timer::GetElapsedCpuTime(void) const {
#ifdef _MSC_VER
    FILETIME CreationTime, ExitTime, KernelTime, UserTime;
    const bool status = GetProcessTimes(GetCurrentProcess(), &CreationTime, &ExitTime, &KernelTime, &UserTime) != 0;
    uint64_t kernelStart, kernelStop, userStart, userStop;

    memcpy((char*)&kernelStart, (char*)&mKernelStartTime, 8);
    memcpy((char*)&kernelStop,  (char*)&KernelTime,       8);
    memcpy((char*)&userStart,   (char*)&mUserStartTime,   8);
    memcpy((char*)&userStop,    (char*)&UserTime,         8);

    const uint64_t kernelElapsed = (kernelStop - kernelStart);
    const uint64_t userElapsed   = (userStop - userStart);
    return (double)(kernelElapsed + userElapsed) / NUM_MICROSECONDS_IN_100NS;
#else
    struct rusage stopResUsage;
    getrusage(RUSAGE_SELF, &stopResUsage);
    const double user_start_usec   = (double)mStartResUsage.ru_utime.tv_sec * NUM_MICROSECONDS_IN_SECOND + (double)mStartResUsage.ru_utime.tv_usec;
    const double user_stop_usec    = (double)stopResUsage.ru_utime.tv_sec * NUM_MICROSECONDS_IN_SECOND + (double)stopResUsage.ru_utime.tv_usec;
    const double system_start_usec = (double)mStartResUsage.ru_stime.tv_sec * NUM_MICROSECONDS_IN_SECOND + (double)mStartResUsage.ru_stime.tv_usec;
    const double system_stop_usec  = (double)stopResUsage.ru_stime.tv_sec * NUM_MICROSECONDS_IN_SECOND + (double)stopResUsage.ru_stime.tv_usec;
    return ((user_stop_usec - user_start_usec) + (system_stop_usec - system_start_usec)) / NUM_MICROSECONDS_IN_SECOND;
#endif
}

// returns a string containing both the elapsed wall time and CPU time
string Timer::GetElapsedTime(void) const {
    ostringstream sb;
    const double cpuTime  = GetElapsedCpuTime();
    const double wallTime = GetElapsedWallTime();

    double ioPercentage = (wallTime - cpuTime) / wallTime * 100.0;
    if(ioPercentage < 0.0) ioPercentage = 0.0;

    sb << "CPU: " << fixed << setprecision(1) << cpuTime
        << " s, wall: " << wallTime << " s (I/O: " << ioPercentage << " %)";
    return sb.str();
}

// returns the elapsed wall time
double Timer::GetElapsedWallTime(void) const {
#ifdef _MSC_VER
    FILETIME wallStopTime;
    GetSystemTimeAsFileTime(&wallStopTime);
    unsigned long long start, stop;
    memcpy((char*)&start, (char*)&mWallStartTime, 8);
    memcpy((char*)&stop, (char*)&wallStopTime, 8);
    return (double)(stop - start) / NUM_MICROSECONDS_IN_100NS;
#else
    struct timeval wallStopTime;
    gettimeofday(&wallStopTime, (struct timezone *)0);
    double start_usec = (double)mWallStartTime.tv_sec * NUM_MICROSECONDS_IN_SECOND + (double)mWallStartTime.tv_usec;
    double stop_usec  = (double)wallStopTime.tv_sec   * NUM_MICROSECONDS_IN_SECOND + (double)wallStopTime.tv_usec;
    return (stop_usec - start_usec) / NUM_MICROSECONDS_IN_SECOND;
#endif
}

// restarts the internal timer
void Timer::Restart(void) {
#ifdef _MSC_VER
    GetSystemTimeAsFileTime(&mWallStartTime);
    FILETIME CreationTime, ExitTime;
    GetProcessTimes(GetCurrentProcess(), &CreationTime, &ExitTime, &mKernelStartTime, &mUserStartTime);
#else
    gettimeofday(&mWallStartTime, (struct timezone *)0);
    getrusage(RUSAGE_SELF, &mStartResUsage);
#endif
}

}
}
