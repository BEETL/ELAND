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

/* c++/include/config.h.cmake. Manually edited */
/* c++/include/config.h.in.  Generated from configure.ac by autoheader.  */

/* Helper macros to get the version string */
/* const std::string version(EXPAND(VERSION));  */

#define STRINGIFY(x) #x
#define EXPAND(x) STRINGIFY(x)

/* Pipeline install path prefix */
/* #undef CASAVA_PIPELINE_DIR */

/* Endianness of the architecture */
/* #undef CASAVA_IS_BIG_ENDIAN */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the `floorf' function. */
#define HAVE_FLOORF 1

/* Define to 1 if you have the `round' function. */
#define HAVE_ROUND 1

/* Define to 1 if you have the `roundf' function. */
#define HAVE_ROUNDF 1

/* Define to 1 if you have the `powf' function. */
#define HAVE_POWF 1

/* Define to 1 if you have the `erf' function. */
#define HAVE_ERF 1

/* Define to 1 if you have the `erff' function. */
#define HAVE_ERFF 1

/* Define to 1 if you have the `erfc' function. */
#define HAVE_ERFC 1

/* Define to 1 if you have the `erfcf' function. */
#define HAVE_ERFCF 1

/* Define to 1 if you have the `zlib' library */
#define HAVE_ZLIB 1

/* Define to 1 if you have the `bzip2' library */
#define HAVE_BZIP2 1
#define HAVE_BZLIB 1

/* Define to 1 if you have the `fftw3f' library */
/* #undef HAVE_FFTW3F */

/* Define to 1 if you have the `cpgplot' library */
/* #undef HAVE_CPGPLOT */

/* Define to 1 if you have the `pgplot' library */
/* #undef HAVE_PGPLOT */

/* Define to 1 if you have the `X11' library */
#define HAVE_X11 1

/* Define to 1 if you have the `g2c' library */
/* #undef HAVE_G2C */

/* Define to 1 if you have the `boost_xxx_yyy' library
   (-lboost_xxx_yyy). */
#define HAVE_LIBBOOST_DATE_TIME 1
#define HAVE_LIBBOOST_FILESYSTEM 1
#define HAVE_LIBBOOST_IOSTREAMS 1
#define HAVE_LIBBOOST_PROGRAM_OPTIONS 1
/* #undef HAVE_LIBBOOST_PYTHON */
#define HAVE_LIBBOOST_REGEX 1
#define HAVE_LIBBOOST_SERIALIZATION 1
#define HAVE_LIBBOOST_SYSTEM 1

/* Define to 1 if you have the `cppunit' library (-lcppunit). */
/* #undef HAVE_CPPUNIT */

/* Name of package */
/* #undef PACKAGE */

/* Top level namespace */
/* #undef NAMESPACE */

/* Define to the address where bug reports for this package should be sent. */
/* #undef PACKAGE_BUGREPORT casava_bugs@illumina.com */

/* Define to the full name of this package. */
/* #undef PACKAGE_NAME */

/* Define to the full name and version of this package. */
/* #undef PACKAGE_STRING */

/* Define to the one symbol short name of this package. */
#undef PACKAGE_TARNAME

/* Define to the version of this package. */
/* #undef PACKAGE_VERSION */

/* Version number of package */
/* #undef VERSION */

/* Define to empty if `const' does not conform to ANSI C. */
#undef const
