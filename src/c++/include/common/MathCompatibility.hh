/**
 ** Copyright (c) 2007-2009 Illumina, Inc.
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
 ** \file MathCompatibility.hh
 **
 ** \brief Compatibility layer for math-related constructs.
 **
 ** \author Come Raczy
 **/

#ifndef CASAVA_COMMON_MATH_COMPATIBILITY_HH
#define CASAVA_COMMON_MATH_COMPATIBILITY_HH

#if defined(__BORLANDC__)
#include <math.h>
#else
#include <cmath>
#endif

#include "config.h"

#ifndef HAVE_FLOORF
inline float floorf(float x)
{
    return static_cast<float> (floor(x));
}
#define HAVE_FLOORF
#endif

#ifndef HAVE_ROUND
inline double round(double x)
{
    return (x - floor(x) < 0.5) ? floor(x) : ceil(x);
}
#define HAVE_ROUND
#endif

#ifndef HAVE_ROUNDF
inline float roundf(float x)
{
    return (x - floorf(x) < 0.5f ? floorf(x) : ceil(x));
}
#define HAVE_ROUNDF
#endif

#ifndef HAVE_POWF
inline float powf(float x, float y)
{
    return static_cast<float> (pow(x, y));
}
#define HAVE_POWF
#endif

#ifndef HAVE_ERF
#include <boost/math/special_functions/erf.hpp>
inline double erf(double x)
{
    return boost::math::erf(x);
}
#define HAVE_ERF
#endif

#ifndef HAVE_ERFF
#include <boost/math/special_functions/erf.hpp>
inline float erff(float x)
{
    return boost::math::erf(x);
}
#define HAVE_ERFF
#endif

#ifndef HAVE_ERFC
#include <boost/math/special_functions/erf.hpp>
inline double erfc(double x)
{
    return boost::math::erfc(x);
}
#define HAVE_ERFC
#endif

#ifndef HAVE_ERFCF
#include <boost/math/special_functions/erf.hpp>
inline float erfc(float x)
{
    return boost::math::erfc(x);
}
#define HAVE_ERFCF
#endif

#endif // #ifndef CASAVA_COMMON_MATH_COMPATIBILITY_HH
