/**
 ** Copyright (c) 2003-2006 Solexa Limited. All rights reserved.
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
 ** \file eland_ms/TableEntry.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 ** Contains interface functionality from ELAND - specifically, stuff to do
 ** with scoring and storage of alignments
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_TABLE_ENTRY_H
#define CASAVA_ELAND_MS_TABLE_ENTRY_H

#include "ElandConstants.hh"

namespace casava
{
namespace eland_ms
{

typedef unsigned short PrefixType;



// separate out portion of data structure without template dependencies:
//
struct TableEntryData
{
  union SuffixType
  {
    Word ui;
    uchar uc[sizeof(Word)];
  }; // ~union SuffixType

  typedef unsigned short MaskTableEntry;

  PrefixType prefix;
  MaskTableEntry mask;
  SuffixType suffix;
  OligoNumber position;
} __attribute__ ((packed)); // ~struct TableEntryData



// TC 18.2.8 - in the split prefix case (where the prefix size exceeds the
// hash key size) we will store the extra bits of separately in 'prefix.'
// This removes the restriction on number of oligos per batch at expense of
// extra 2 bytes per TableEntry. However may have some processing advantages.

// Each Partition HashTable consists of a look up table of pointers into a list
// of table entries. Templatized because < and == operators need to change
// in split prefix mode
template <bool useSplitPrefix> struct TableEntry : public TableEntryData
{
  bool operator<(const TableEntry<useSplitPrefix>& rhs) const;
  bool operator==(const TableEntry<useSplitPrefix>& rhs) const;
  bool operator!=(const TableEntry<useSplitPrefix>& rhs) const;

} __attribute__ ((packed)); // ~struct TableEntry


template<> inline bool TableEntry<false>::operator<(const TableEntry<false>& rhs) const
{
  return ((mask<rhs.mask)
          ||((mask==rhs.mask)
             &&((suffix.ui<rhs.suffix.ui)
                ||((suffix.ui==rhs.suffix.ui)&&(position<rhs.position)))));
}

template<> inline bool TableEntry<false>::operator==(const TableEntry<false>& rhs) const
{
  return ( (suffix.ui==rhs.suffix.ui)
           &&(mask==rhs.mask)
           &&(position==rhs.position));
}

template<> inline bool TableEntry<false>::operator!=(const TableEntry<false>& rhs) const
{
  return ( (suffix.ui!=rhs.suffix.ui)
           ||(mask!=rhs.mask)
           ||(position!=rhs.position));
}

template<> inline bool TableEntry<true>::operator<(const TableEntry<true>& rhs) const
{
  //  return ( ((position&splitPrefixMaskHigh)<(rhs.position&splitPrefixMaskHigh))
  //       ||(((position&splitPrefixMaskHigh)==(rhs.position&splitPrefixMaskHigh))
  return ( (prefix<rhs.prefix)
           ||((prefix==rhs.prefix)
              &&((mask<rhs.mask)
                 ||((mask==rhs.mask)
                    &&((suffix.ui<rhs.suffix.ui)
                       ||((suffix.ui==rhs.suffix.ui)
                          &&(position<rhs.position)))))));
} // ~bool operator<(const TableEntry& rhs)

template<> inline bool TableEntry<true>::operator==(const TableEntry<true>& rhs) const
{
  return ( (suffix.ui==rhs.suffix.ui)
           &&(mask==rhs.mask)
           &&(prefix==rhs.prefix)
           &&(position==rhs.position));

  //       &&(position==rhs.position));

}

template<> inline bool TableEntry<true>::operator!=(const TableEntry<true>& rhs) const
{
  return ( (suffix.ui!=rhs.suffix.ui)
           ||(mask!=rhs.mask)
           ||(prefix!=rhs.prefix)
           ||(position!=rhs.position));

  //       ||(position!=rhs.position));

}

} // namespace eland_ms
} // namespace casava



#endif // CASAVA_ELAND_MS_TABLE_ENTRY_H
