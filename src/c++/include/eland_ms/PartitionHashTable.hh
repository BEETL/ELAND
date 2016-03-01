/**
 ** copyright (c) 2007-2010 Illumina, Inc.
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
 ** \file eland_ms/PartitionHashTable.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_PARTITION_HASH_TABLE_H
#define CASAVA_ELAND_MS_PARTITION_HASH_TABLE_H

#include <boost/format.hpp>

#include "common/Exceptions.hh"

#include "pht/HelperFwd.hh"
#include "pht/HelperRvrs.hh"

namespace casava
{
namespace eland_ms
{

// class PartitionHashTable
// This stores all the oligo information for one of the six partitions
// The partitions are searched two at a time, thus three runs through the
// genome are required in total
template< bool useSplitPrefix, int PASS, bool isFwd> class PartitionHashTable
{
public:
  PartitionHashTable( HashTableDataStore<useSplitPrefix>& data, MatchTable& results ) :
    data_(data), sps_(data, results), results_(results)
  {} // ~ctor

//  void setTable( int prefixLength,
//		 int lowerFragSize,
//		 Word lowerFragMask,
//		 const FragmentErrorType* lowerFragScore,
//		 const FragmentErrorType* upperFragScore );
//
//  void makePointerArray( void );
//  //  void removeRepeatedEntries( vector<MatchPosition>& matchPosition_ );__
//  void removeRepeatedEntries( MatchTable& results );

  // private:
  // data storage for hash table
  int numBits_;

  HashTableDataStore<useSplitPrefix>& data_;
  // Moved to HashTableDataStore
  //  vector<TablePointer> entryPointer_;
  //  vector<TableEntry<useSplitPrefix> > hashRem_;
  //  vector<OligoNumber> position_;

  PHTHelper<PASS, isFwd, useSplitPrefix> sps_;
  //  TablePointer* pCount_;
  //  int lowerFragSize_;
  //  Word lowerFragMask_;
  //  const FragmentErrorType* lowerFragScore_;
  //  const FragmentErrorType* upperFragScore_;

  MatchTable& results_;


void setTable
( int prefixLength,
  int lowerFragSize,
  Word lowerFragMask,
  const FragmentErrorType* lowerFragScore,
  const FragmentErrorType* upperFragScore )
{
  if ((2*prefixLength)>MAX_HASH_BITS)
  {
    numBits_= MAX_HASH_BITS;
    sps_.splitPrefixShift_= (2*prefixLength)-MAX_HASH_BITS;
    cerr << "partition hash table will run in split prefix mode, split is "
	 << numBits_ << "/" << sps_.splitPrefixShift_ << endl;
    sps_.splitPrefixMask_=hashmask(sps_.splitPrefixShift_);
  } // ~if
  else
  {
    numBits_=2*prefixLength;
    // shift and mask still need to be explicitly set to zero because this
    // branch can still be executed when useSplitPrefix is true: if
    // oligo length is 26 and max bits is 25, then on pass 2 there is one
    // 25 bit table and 1 24 bit table TC 27.02.04
    sps_.splitPrefixShift_=0;
    sps_.splitPrefixMask_=0;
  } // ~else
  const bool needToClear(! data_.entryPointer_.empty());
  data_.entryPointer_.resize( 2 + (1<<numBits_) );
  if (needToClear)
    memset(&data_.entryPointer_[0], '\0',
	   sizeof(TablePointer)*data_.entryPointer_.size());

  sps_.lowerFragSize_ = lowerFragSize;
  sps_.lowerFragMask_ = lowerFragMask;
  sps_.lowerFragScore_ = lowerFragScore;
  sps_.upperFragScore_ = upperFragScore;
  sps_.pCount_ = &data_.entryPointer_[2];

  cerr << "Setting partition hash table size to "
       << numBits_ << " bits" << endl;
} // ~void setTable( int numBits )



void makePointerArray(MaskMapType& maskMap)
{
  //
  // transform entryPointer table entries from individual to sub-total counts:
  //
  {
    typedef std::vector<TablePointer>::iterator titer;
    titer i(data_.entryPointer_.begin()+3), i_end(data_.entryPointer_.end());
    for(;i!=i_end;++i) { *i += *(i-1); }
  }
  cerr << "will place " << data_.entryPointer_.back() << " entries in table"<< endl;

  const bool needToClear(! data_.hashRem_.empty());
  data_.hashRem_.resize( data_.entryPointer_.back());

  if (needToClear)
  {
    memset(&data_.hashRem_[0], '\0', sizeof(TableEntry<useSplitPrefix>)*data_.hashRem_.size());
  } // ~if

  sps_.pCount_--; // now pCount_ == data_.entryPointer_ + 1

  //
  // take all of the masks found and put them into a faster look-up:
  //

  // guard against the unlikely event that there were no unmasked suffices:
  // (the 0 state is required for the current hash lup function to work correctly)
  maskMap[0] = 0;

  assert(maskMap.size()<65535);
  sps_.maskTable_.reserve(maskMap.size());
  {
    MaskMapType::iterator i(maskMap.begin()), i_end(maskMap.end());
    for(uint j(0); i!=i_end; ++i, ++j)
    {
      i->second=j;
      sps_.maskTable_.push_back(i->first);
    }
  }
} // ~void PartitionHashTable::makePointerArray( void )

void removeRepeatedEntries ( MatchTable& results )
{
  sps_.pCount_--; // now sps_.pCount_ == data_.entryPointer_

  if (data_.hashRem_.empty()) {
      // If there are no entries there can be no repeated entries.
      // (Not bailing until after the pCount_ decrement.)
      return;
  }

  int numSaved(0);

  TableEntry<useSplitPrefix> lastEntry;
  //  Word lastEntry;

  uint k,l(0);
  OligoNumber existingOligo, newOligo;

  // remove repeats from table

  // Sort entries for each hash value into order then eliminate any repeats.
  // This bit is based on code in IndexGenome.cpp

  const uint32_t tableSize(1<<numBits_);
  for ( uint32_t i(0) ; i < tableSize ; i++ )
  {

    //    cerr << "A" << endl;

    if (sps_.pCount_[i+1]-sps_.pCount_[i]>1)
    {
      //      cerr << "### " <<  (sps_.pCount_[i+1]-sps_.pCount_[i]) << endl;
      sort( data_.hashRem_.begin()+sps_.pCount_[i], data_.hashRem_.begin() + sps_.pCount_[i+1] );
      //   mirrorSort( hashRem_.begin(),
      //	  hashRem_.begin()+pCount_[i],
      //	  hashRem_.begin() + pCount_[i+1],
      //	  swapper );

    } // ~if

    k = sps_.pCount_[i];
    sps_.pCount_[i] = l;
    //    cout << i << ": " << l << " from " << k << endl;
    //    lastEntry = hashRem_[k].suffix.ui;
    //   lastEntry ^= ~0; // ensures false on first iteration
    if (k < data_.hashRem_.size())
    {
        lastEntry = data_.hashRem_[k];
    }
    else
    {
        using boost::format;
        using casava::common::CasavaException;
        const std::string message = (format("unexpected sps_.pCount[%i] == %i (data_.hashRem_.size() == %d)") % i % k % data_.hashRem_.size()).str();
        //BOOST_THROW_EXCEPTION(CasavaException(EINVAL, message));
        // get the last element of the table instead of reading past the end
        //std::cerr << "WARNING: " << message << std::endl;
        lastEntry = data_.hashRem_.back();
    }

    lastEntry.suffix.ui ^= ~0; // ensures false on first iteration


    for ( ; k!= sps_.pCount_[i+1] ; k++ )
    {
      if (data_.hashRem_[k].mask==lastEntry.mask) numSaved++;
      if (data_.hashRem_[k] != lastEntry)
      { // no repeat, copy as normal
	//	entry_[l]=entry_[k];
	data_.hashRem_[l]=data_.hashRem_[k];
	//	position_[l]=position_[k];
	//	cout << i << ": "<< hashRem_[l].suffix.ui << " to "<< l<< " from "<< k<< endl;
	++l;

      } // ~if
      else
      { // repeat

	--l; // shift back to first repeated entry
	//	type_[l]=repeatInBatch;
	//	entry_[l]=entry_[k];
	//	entry_[l].position_ = repeatFlag + repeatInGenome + 1;
	do
	{

	  // flag current entry as being same as that entry
	  //	  oligos_.flagRepeat( oligoNum_[l]&(~reverseComplementFlag),
	  //	      oligoNum_[k]&(~reverseComplementFlag));
	  //cout<<((hashRem_[k].position&splitPrefixMaskLow)&(~isReverseOligo))
	    //	       << ((hashRem_[k].position&isReverseOligo)?'R':'F')
	  //   << " is the same as "
	  // << ((hashRem_[l].position&splitPrefixMaskLow)&(~isReverseOligo))
	  //   << ((hashRem_[l].position&isReverseOligo)?'R':'F')
	  //   << endl;
	  existingOligo
	    = (data_.hashRem_[l].position)&(~isReverseOligo);
	  newOligo
	    = (data_.hashRem_[k].position)&(~isReverseOligo);
          results.setSameAs__( newOligo, existingOligo );
#ifdef ZZZZ
	  if (  (matchPosition_[newOligo]<blockRepeat)
		&& (matchPosition_[existingOligo]!=blockRepeat+newOligo) )
	  {
	    matchPosition_[newOligo] = blockRepeat+existingOligo;
	  } // ~if
#endif
	    //	  matchPosition_[position_[k]&(~isReverseOligo)]
	    //   =blockRepeat+(position_[l]&(~isReverseOligo));

	  k++;
	} // ~do
	while( (data_.hashRem_[k]==lastEntry)&&(k!=sps_.pCount_[i+1]) );
	k--; // need to look at next entry again
	l++; // and shift destination pointer
      } // ~else
      lastEntry = data_.hashRem_[k];
    } // ~if
  } // ~for i

  cerr << sps_.pCount_[tableSize] << " entries reduced to " << l << " entries"
       << endl;

  //  numOligos_ = pCount_[tableSize];
  sps_.pCount_[tableSize] = l;

  data_.hashRem_.resize(l);
  //  position_.resize(l);

  // add first prefix into the top bits of the table entry pointers 
  //
  sps_.setTopPrefix(tableSize);
}

};

} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_PARTITION_HASH_TABLE_H
