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
 ** \file eland_ms/MatchTable.cpp
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#include "eland_ms/Hasher.hh"

namespace casava
{
namespace eland_ms
{

template<> void Hasher<0, 32>::operator()( const Oligo& ol, Oligo& out ) const
{
    do_interspersed( ol,out );
}

template<> void Hasher<1, 32>::operator()( const Oligo& ol, Oligo& out ) const
{
    Oligo r;
    do_interspersed( ol,r );

    out.ui[0]=r.ui[0]&HasherBase<32>::fragMaskB_; // out.ui[0] = 0B
    out.ui[0]|=(r.ui[1]&HasherBase<32>::fragMaskC_);// out.ui[0] = CB
    out.ui[1]=r.ui[1]&HasherBase<32>::fragMaskD_; // out.ui[1] = 0D
//  out.ui[1]|=(ol.ui[1]&HasherBase<32>::fragMaskA_); // out.ui[1] = AD - oops! AC 6.6.3
    out.ui[1]|=(r.ui[0]&HasherBase<32>::fragMaskA_); // out.ui[1] = AD
}

template<> void Hasher<2, 32>::operator()( const Oligo& ol, Oligo& out ) const
{
    Oligo r;
    do_interspersed( ol,r );

    out.ui[0]=r.ui[0]&HasherBase<32>::fragMaskB_; // out.ui[0] = 0B
    out.ui[0]<<=HasherBase<32>::hashShift2_; // shift B if necessary so D fits next to it
    out.ui[0]|=r.ui[1]>>(numBitsPerBase*ElandConstants<32>::fragLengthC); // out.ui[0]=DB

    // ordering of A and C swapped in second fragment TC 03.07.03
    // out.ui[1]&=fragMaskC_; // out.ui[1]=C0
    // out.ui[1]|=((ol.ui[0]&fragMaskA_)<<(numBitsPerBase*fragLengthC));
    out.ui[1]=r.ui[1]&HasherBase<32>::fragMaskC_;
    out.ui[1]<<=(numBitsPerBase*ElandConstants<32>::fragLengthA);
    out.ui[1]|=(r.ui[0]&HasherBase<32>::fragMaskA_);

    //  printWord(out.ui[1],16);
    //  cout << '-';
    //  printWord(out.ui[0],16);
    //  cout << endl;
}

} //namespace eland_ms
} //namespace casava
