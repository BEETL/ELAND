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
 ** \file eland_ms/QueryGenerator.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_QUERY_GENERATOR_H
#define CASAVA_ELAND_MS_QUERY_GENERATOR_H

#include "ElandDefines.hh"
#include "ReverseShifter.hh"

namespace casava
{
namespace eland_ms
{

// Oligo format, e.g. for 21 mers, 1 character = 1 bit:
// ui0.............................ui1.............................
// uc0.....uc1.....uc2.....uc3.....uc4.....uc5.....uc6.....uc7.....
// 2.1.1.1.1.1.1.1.1.1.1.9.8.7.6.5.4.3.2.1.0.  <-base numbers 0 20
// 0 9 8 7 6 5 4 3 2 1 0

// QueryGenerator: generates a set of query sequences for an oligo, allowing
// for:
// i)  reverse complements
// ii) Ns
// iii) 2 base ambiguity codes (TBD)

template <int OLIGO_LEN> class QueryGenerator
{
public:

  QueryGenerator()
  {}

//  // generate a set of query sequences from ASCII sequence data
//  int operator()( const char* buf, const uint oligoNum,
//		  vector<Oligo>& queryOligo,
//		  vector<Oligo>& queryMask,
//		  vector<OligoNumber>&  queryOligoNum );
//
//  // convert a single oligo from ASCII to binary.
//  // Fail if the oligo contains Ns or other funny stuff
//  void operator()( const char* buf, Oligo& oligo );
//
//  // convert an oligo from ASCII to binary and produce a
//  // mask for the Ns - return number of Ns found
//  int encodeOligo( const char* buf, Oligo& oligo, Oligo& mask );

  void reverseOligo( const Oligo& ol, Oligo& rc ) const
  {

    rc.uc[0]=reverseChar[ol.uc[7]];
    rc.uc[1]=reverseChar[ol.uc[6]];
    rc.uc[2]=reverseChar[ol.uc[5]];
    rc.uc[3]=reverseChar[ol.uc[4]];
    rc.uc[4]=reverseChar[ol.uc[3]];
    rc.uc[5]=reverseChar[ol.uc[2]];
    rc.uc[6]=reverseChar[ol.uc[1]];
    rc.uc[7]=reverseChar[ol.uc[0]];
    reverseShift_( ol, rc );
  } // ~reverseOligo

protected:
  char tempBuf1[OLIGO_LEN+1];
  char tempBuf2[OLIGO_LEN+2];
  ReverseShifter<ElandConstants<OLIGO_LEN>::prefixLength-ElandConstants<OLIGO_LEN>::suffixLength, OLIGO_LEN> reverseShift_;

//  int convert_( const char* buf, const uint oligoNum,
//		vector<Oligo>& queryOligo,
//		vector<Oligo>& queryMask,
//		vector<OligoNumber>&  queryOligoNum, const short& seed_no );




public:
// convert an oligo from ASCII to binary and produce a
// mask for the Ns - return number of Ns found
// template< int oligoLength >
int encodeOligo( const char* buf, Oligo& oligo, Oligo& mask )
{
  int numNs(0);
  bool isN;

  for (int i(0);i<ElandConstants<OLIGO_LEN>::prefixLength;++i, ++buf)
  {
    isN=isBlank(*buf);
    //    isN=( (*buf=='n')||(*buf=='N') );
    oligo.ui[1]<<=numBitsPerBase;
    oligo.ui[1]|=(isN?0:whichBase[(uint)*buf]);
    mask.ui[1]<<=numBitsPerBase;
    mask.ui[1]|=0x3*isN;
    numNs+=1*isN;
  } // ~for i

  for (int i(0);i<ElandConstants<OLIGO_LEN>::suffixLength;i++, ++buf) // was max bases per word
  {
    isN=isBlank(*buf);
    //    isN=( (*buf=='n')||(*buf=='N') );
    oligo.ui[0]<<=numBitsPerBase;
    oligo.ui[0]|=(isN?0:whichBase[(uint)*buf]);
    mask.ui[0]<<=numBitsPerBase;
    mask.ui[0]|=0x3*isN;
    numNs+=1*isN;
  } // ~for i

  return numNs;
} // ~int encodeOligo( const char* buf, Oligo& oligo, Oligo& mask )


// convert a single oligo from ASCII to binary.
// Fail if the oligo contains Ns or other funny stuff
// template< int oligoLength >
void operator()( const char* buf, Oligo& oligo )
{
  Oligo dummy;
  oligo.ui[0]=0;
  oligo.ui[1]=0;
  if (encodeOligo( buf, oligo, dummy)!=0)
  {
    cerr << "Error: did not expect to find Ns in oligo!" << endl;
  }
} // ~void operator()( const char* buf, Oligo& oligo )

#ifdef ORIGINAL_CODE
// generate a set of query sequences from ASCII sequence data
// template< int oligoLength >
int QueryGenerator::operator()
( const char* buf, const uint oligoNum,
  vector<Oligo>& queryOligo,
  vector<Oligo>& queryMask,
  vector<OligoNumber>&  queryOligoNum )
{
  int numNs(0), tailSize(0), headSize(0);
  //  char tempBuf1[OLIGO_LENGTH+1], tempBuf2[OLIGO_LENGTH+1];
  tempBuf1[OLIGO_LENGTH]='\0';
  tempBuf2[OLIGO_LENGTH]='\0';
  queryOligo.clear();
  queryMask.clear();
  queryOligoNum.clear();

  //  cout << endl << buf << " - initial " << endl;

  for (int i(0);i<OLIGO_LENGTH;i++)
    numNs+=isBlank(buf[i]);

  //  if (numNs>4)
  //  {
  //    cout << "Too many Ns, rejecting" << endl;
  //   return numNs;
  //  } // ~if


  //  while ( ((buf[headSize]=='N')||(buf[headSize]=='n'))
  //	  &&(headSize<OLIGO_LENGTH) )
  while ( isBlank(buf[headSize])&&(headSize<OLIGO_LENGTH) )
    headSize++;

  // Catch the all Ns case - otherwise head and tail will fully overlap
  // producing bogus results.
  if (headSize == OLIGO_LENGTH) {
      // cout << "Seed is all Ns, rejecting" << endl;
      return OLIGO_LENGTH;
  }

  //  while( ((buf[OLIGO_LENGTH-tailSize-1]=='N')
  //	  ||(buf[OLIGO_LENGTH-tailSize-1]=='n'))&&(tailSize<OLIGO_LENGTH))
  while( (isBlank(buf[OLIGO_LENGTH-tailSize-1])&&(tailSize<OLIGO_LENGTH)))
    tailSize++;

  // copy to buffer, shifting any head->tail

  for ( int i(0);i<OLIGO_LENGTH-headSize;i++)
    tempBuf1[i]=buf[i+headSize];

  for ( int i(OLIGO_LENGTH-headSize);i<OLIGO_LENGTH;i++)
    tempBuf1[i]='N';

  tailSize += headSize; // can now forget about headSize

  int numInternalNs(numNs-tailSize);

  // Changes to N handling - previously max number of Ns had a hard-coded stop
  // at four bases. Fundamental limitation is that at most two of the four
  // fragments each read is split into can contain either errors or Ns.
  // So the possibilities are:
  // i. No Ns - will find any combination of two substitution errors
  // ii. Two internal Ns
  // ii. Trailing Ns fit inside one fragment: one N or error will be found
  // iii. Trailing Ns fit inside two fragments: read must match exactly
  if ((numInternalNs>2)
      ||((numInternalNs==2)&&(tailSize!=0))
      ||((numInternalNs==1)&&(tailSize>fragLengthA))
      ||((numInternalNs==0)&&(tailSize>fragLengthA+fragLengthD)))
  {
    // cout << "Too many separate Ns, rejecting" << endl;
    return numNs;
  }

  // OK, output oligo as is
  //  cout << tempBuf1 << " - shifted " << endl;
  queryOligo.push_back( Oligo() );
  queryMask.push_back( Oligo() );
  queryOligoNum.push_back(oligoNum);

  assert( numNs==encodeOligo(tempBuf1, queryOligo.back(), queryMask.back()));
  //    printWord(queryOligo.back().ui[1],prefixLength_);
  //    printWord(queryOligo.back().ui[0],16);
  //   cout << " - encoded oligo" << endl;
  //    printWord(queryMask.back().ui[1],prefixLength_);
  //    printWord(queryMask.back().ui[0],16);
  //  cout << " - encoded mask" << endl;

#ifdef ALLOW_N_TO_BE_DELETION
  // if >=1 error, copy removing first N, send to output list
  if (numInternalNs>=1)
  {
    int i, firstN, secondN;

    tempBuf2[OLIGO_LENGTH-1]='n';
    //    for (i=0;((i<OLIGO_LENGTH-tailSize)
    //      &&(tempBuf1[i]!='n')
    //      &&(tempBuf1[i]!='N'));i++)
    for (i=0;((i<OLIGO_LENGTH-tailSize)&&(!isBlank(tempBuf1[i])));i++)
    {
      tempBuf2[i]=tempBuf1[i];
    } // ~for i
    assert(i!=OLIGO_LENGTH-tailSize);
    firstN=i;
    i++;

    for (;(i<OLIGO_LENGTH);i++)
    {
      tempBuf2[i-1]=tempBuf1[i];
    } // ~for i

    //     cout << "After removal of first internal N: " << endl;
    //   cout << tempBuf2 << " - modified buffer" << endl;
    // OK, output oligo
    queryOligo.push_back( Oligo() );
    queryMask.push_back( Oligo() );
    queryOligoNum.push_back(oligoNum);

    assert( numNs==encodeOligo(tempBuf2, queryOligo.back(), queryMask.back()));
    assert( (queryMask.back().ui[0]&0x3)==0x3 );
    //  queryOligo.back().ui[0] |= nonIncorpFirstN;
    queryOligo.back().ui[0] |= (nonIncorpFirstN&queryMask.back().ui[0]);
    queryOligo.back().ui[1] |= (nonIncorpFirstN&queryMask.back().ui[1]);

    //      printWord(queryOligo.back().ui[1],prefixLength_);
    //   printWord(queryOligo.back().ui[0],16);
    //   cout << " - encoded ol " << endl;
    //     printWord(queryMask.back().ui[1],prefixLength_);
    //    printWord(queryMask.back().ui[0],16);
    //   cout << endl;

    if (numInternalNs==2)
    {
      //      for (i=firstN+1;((i<OLIGO_LENGTH-tailSize)
      //	       &&(tempBuf1[i]!='n')
      //	       &&(tempBuf1[i]!='N'));i++);
      for (i=firstN+1;((i<OLIGO_LENGTH-tailSize)&(!isBlank(tempBuf1[i])));i++);

      assert(i!=OLIGO_LENGTH-tailSize);
      secondN=i;

      for (i=0;i<firstN;i++) tempBuf2[i]=tempBuf1[i];
      for (i=firstN+1;i<secondN;i++) tempBuf2[i-1]=tempBuf1[i];
      for (i=secondN+1;i<OLIGO_LENGTH;i++) tempBuf2[i-2]=tempBuf1[i];
      tempBuf2[OLIGO_LENGTH-2]='n';
      tempBuf2[OLIGO_LENGTH-1]='n';
      //      cout << "After removal of both internal Ns: " << endl;
      //  cout << tempBuf2 << " - buffer" << endl;
      queryOligo.push_back( Oligo() );
      queryMask.push_back( Oligo() );
      queryOligoNum.push_back(oligoNum);

      assert( numNs
	      ==encodeOligo(tempBuf2, queryOligo.back(), queryMask.back()));
      assert( (queryMask.back().ui[0]&0x3)==0x3 );
      //      queryOligo.back().ui[0] |= nonIncorpBothNs;
      queryOligo.back().ui[0] |= (nonIncorpBothNs&queryMask.back().ui[0]);
      queryOligo.back().ui[1] |= (nonIncorpBothNs&queryMask.back().ui[1]);


      //         printWord(queryOligo.back().ui[1],prefixLength_);
      //    printWord(queryOligo.back().ui[0],16);
      //   cout << " - encoded " << endl;
      //    printWord(queryMask.back().ui[1],prefixLength_);
      //   printWord(queryMask.back().ui[0],16);
      //  cout << endl;

      if (secondN!=firstN+1)
      {
	for (i=0;i<secondN;i++) tempBuf2[i]=tempBuf1[i];
	for (i=secondN+1;i<OLIGO_LENGTH;i++) tempBuf2[i-1]=tempBuf1[i];
	tempBuf2[OLIGO_LENGTH-1]='n';
	//	 cout << "After removal of second internal N:" << endl;
	//	 cout << tempBuf2 << " - buffer" << endl;
	queryOligo.push_back( Oligo() );
	queryMask.push_back( Oligo() );
	queryOligoNum.push_back(oligoNum);

	assert( numNs
		==encodeOligo(tempBuf2, queryOligo.back(), queryMask.back()));
	assert( (queryMask.back().ui[0]&0x3)==0x3 );
	//	queryOligo.back().ui[0] |= nonIncorpSecondN;
	queryOligo.back().ui[0] |= (nonIncorpSecondN&queryMask.back().ui[0]);
	queryOligo.back().ui[1] |= (nonIncorpSecondN&queryMask.back().ui[1]);

	//		printWord(queryOligo.back().ui[1],prefixLength_);
	//		printWord(queryOligo.back().ui[0],16);
	//	 cout << " - encoded" << endl;
	//	printWord(queryMask.back().ui[1],prefixLength_);
	//	printWord(queryMask.back().ui[0],16);
	//	 cout << endl;

      } // ~if
      else
      {
	//	 cout << "Internal Ns are adjacent" << endl;
      }
   } // ~if
  } // ~if
#endif // ~ifdef ALLOW_N_TO_BE_DELETION

#ifndef DONT_SEARCH_REVERSE_STRAND
  // add reverse to all oligos in pile
  const int numToDo(queryOligo.size());
  for (int i(0);i<numToDo;++i)
  {
    queryOligo.push_back( Oligo() );
    queryMask.push_back( Oligo() );
    queryOligoNum.push_back(queryOligoNum[i]|isReverseOligo);
    reverseOligo( queryOligo[i], queryOligo.back());

    // check for self-complementary sequences. These are rare but can
    // cause the oligo concerned to be spuriously flagged as a repeat
    if (queryOligo[i]!=queryOligo.back())
    { // not self-complementary, carry on
      reverseOligo( queryMask[i], queryMask.back());
      //      queryMask.back().ui[0] ^= ((uint)~0); // TBD
      // changed 03.06.03 TC
      // rhs of next line is supposed to contain all ones in its first
      // 2*prefixLength_ bits - does not work for prefixLength_=16
      //	queryMask.back().ui[1] ^= ~( ((uint)~0)<< (prefixLength_<<1) );
      //      queryMask.back().ui[1] ^=
      //	( ((uint)~0) >> ((maxBasesPerWord-prefixLength_)<<1) );
      // Changed to match new oligo format TC 05.08.03

      queryMask.back().ui[0] ^=
	( ((uint)~0) >> ((maxBasesPerWord-suffixLength)<<1) );
      queryMask.back().ui[1] ^=
	( ((uint)~0) >> ((maxBasesPerWord-prefixLength)<<1) );

    } // ~if
    else
    { // self-complementary
      queryOligo.pop_back();
      queryMask.pop_back();
      queryOligoNum.pop_back();
      //	cout << oligoNum << " is self-complementary" << endl;
    } // ~else

  } // ~for
#endif // ~ifndef DONT_SEARCH_REVERSE_STRAND

  return numNs;

} // ~int QueryGenerator::operator()
#endif // ORIGINAL_CODE


// generate a set of query sequences from ASCII sequence data
// template< int oligoLength >
int operator()
( const char* buf, const uint oligoNum,
  vector<Oligo>& queryOligo,
  vector<Oligo>& queryMask,
  vector<OligoNumber>&  queryOligoNum )
{

  // set the stage
  queryOligo.clear();
  queryMask.clear();
  queryOligoNum.clear();

  // _convert now encapsulates the functionality below
  return convert_( buf,
		   oligoNum,
		   queryOligo,
		   queryMask,
		   queryOligoNum,0 );

} // ~int QueryGenerator::operator()

protected:
// private method that encapsulates the main functionality of converting
// buf into oligoNum
int convert_( const char* buf, const uint oligoNum,
			      vector<Oligo>& queryOligo,
			      vector<Oligo>& queryMask,
			      vector<OligoNumber>&  queryOligoNum, const short& seed_no )
{

  // the following should be encapsulated
  int numNs(0), tailSize(0), headSize(0);
  //  char tempBuf1[OLIGO_LENGTH+1], tempBuf2[OLIGO_LENGTH+1];

  tempBuf1[OLIGO_LEN]='\0';
  tempBuf2[OLIGO_LEN]='\0';

  //  cout << endl << buf << " - initial " << endl;

  for (int i(0);i<OLIGO_LEN;i++)
    numNs+=isBlank(buf[i]);

  //  if (numNs>4)
  //  {
  //    cout << "Too many Ns, rejecting" << endl;
  //   return numNs;
  //  } // ~if


  //  while ( ((buf[headSize]=='N')||(buf[headSize]=='n'))
  //	  &&(headSize<OLIGO_LEN) )
  while ( isBlank(buf[headSize])&&(headSize<OLIGO_LEN) )
    headSize++;

  // Catch the all Ns case - otherwise head and tail will fully overlap
  // producing bogus results.
  if (headSize == OLIGO_LEN) {
      // cout << "Seed is all Ns, rejecting" << endl;
      return OLIGO_LEN;
  }

  //  while( ((buf[OLIGO_LEN-tailSize-1]=='N')
  //	  ||(buf[OLIGO_LEN-tailSize-1]=='n'))&&(tailSize<OLIGO_LEN))
  while( (isBlank(buf[OLIGO_LEN-tailSize-1])&&(tailSize<OLIGO_LEN)))
    tailSize++;

  // copy to buffer, shifting any head->tail

  for ( int i(0);i<OLIGO_LEN-headSize;i++)
    tempBuf1[i]=buf[i+headSize];

  for ( int i(OLIGO_LEN-headSize);i<OLIGO_LEN;i++)
    tempBuf1[i]='N';

  tailSize += headSize; // can now forget about headSize

  int numInternalNs(numNs-tailSize);

  // Changes to N handling - previously max number of Ns had a hard-coded stop
  // at four bases. Fundamental limitation is that at most two of the four
  // fragments each read is split into can contain either errors or Ns.
  // So the possibilities are:
  // i. No Ns - will find any combination of two substitution errors
  // ii. Two internal Ns
  // ii. Trailing Ns fit inside one fragment: one N or error will be found
  // iii. Trailing Ns fit inside two fragments: read must match exactly
  if ((numInternalNs>2)
      ||((numInternalNs==2)&&(tailSize!=0))
      ||((numInternalNs==1)&&(tailSize>ElandConstants<OLIGO_LEN>::fragLengthA))
      ||((numInternalNs==0)&&(tailSize>ElandConstants<OLIGO_LEN>::fragLengthA+ElandConstants<OLIGO_LEN>::fragLengthD)))
  {
    // cout << "Too many separate Ns, rejecting" << endl;
    return numNs;
  }

  // OK, output oligo as is
  //  cout << tempBuf1 << " - shifted " << endl;
  queryOligo.push_back( Oligo() );
  queryMask.push_back( Oligo() );
  queryOligoNum.push_back( ((oligoNum|seed_bits[seed_no])) );

  assert( numNs==encodeOligo(tempBuf1, queryOligo.back(), queryMask.back()));
  //    printWord(queryOligo.back().ui[1],prefixLength_);
  //    printWord(queryOligo.back().ui[0],16);
  //   cout << " - encoded oligo" << endl;
  //    printWord(queryMask.back().ui[1],prefixLength_);
  //    printWord(queryMask.back().ui[0],16);
  //  cout << " - encoded mask" << endl;

#ifdef ALLOW_N_TO_BE_DELETION
  // if >=1 error, copy removing first N, send to output list
  if (numInternalNs>=1)
  {
    int i, firstN, secondN;

    tempBuf2[OLIGO_LEN-1]='n';
    //    for (i=0;((i<OLIGO_LENGTH-tailSize)
    //      &&(tempBuf1[i]!='n')
    //      &&(tempBuf1[i]!='N'));i++)
    for (i=0;((i<OLIGO_LENGTH-tailSize)&&(!isBlank(tempBuf1[i])));i++)
    {
      tempBuf2[i]=tempBuf1[i];
    } // ~for i
    assert(i!=OLIGO_LENGTH-tailSize);
    firstN=i;
    i++;

    for (;(i<OLIGO_LENGTH);i++)
    {
      tempBuf2[i-1]=tempBuf1[i];
    } // ~for i

    //     cout << "After removal of first internal N: " << endl;
    //   cout << tempBuf2 << " - modified buffer" << endl;
    // OK, output oligo
    queryOligo.push_back( Oligo() );
    queryMask.push_back( Oligo() );
    queryOligoNum.push_back(oligoNum|seed_bits[seed_no]);

    assert( numNs==encodeOligo(tempBuf2, queryOligo.back(), queryMask.back()));
    assert( (queryMask.back().ui[0]&0x3)==0x3 );
    //  queryOligo.back().ui[0] |= nonIncorpFirstN;
    queryOligo.back().ui[0] |= (nonIncorpFirstN&queryMask.back().ui[0]);
    queryOligo.back().ui[1] |= (nonIncorpFirstN&queryMask.back().ui[1]);

    //      printWord(queryOligo.back().ui[1],prefixLength_);
    //   printWord(queryOligo.back().ui[0],16);
    //   cout << " - encoded ol " << endl;
    //     printWord(queryMask.back().ui[1],prefixLength_);
    //    printWord(queryMask.back().ui[0],16);
    //   cout << endl;

    if (numInternalNs==2)
    {
      //      for (i=firstN+1;((i<OLIGO_LENGTH-tailSize)
      //	       &&(tempBuf1[i]!='n')
      //	       &&(tempBuf1[i]!='N'));i++);
      for (i=firstN+1;((i<OLIGO_LENGTH-tailSize)&(!isBlank(tempBuf1[i])));i++);

      assert(i!=OLIGO_LENGTH-tailSize);
      secondN=i;

      for (i=0;i<firstN;i++) tempBuf2[i]=tempBuf1[i];
      for (i=firstN+1;i<secondN;i++) tempBuf2[i-1]=tempBuf1[i];
      for (i=secondN+1;i<OLIGO_LENGTH;i++) tempBuf2[i-2]=tempBuf1[i];
      tempBuf2[OLIGO_LENGTH-2]='n';
      tempBuf2[OLIGO_LENGTH-1]='n';
      //      cout << "After removal of both internal Ns: " << endl;
      //  cout << tempBuf2 << " - buffer" << endl;
      queryOligo.push_back( Oligo() );
      queryMask.push_back( Oligo() );
      queryOligoNum.push_back(oligoNum|seed_bits[seed_no]);

      assert( numNs
	      ==encodeOligo(tempBuf2, queryOligo.back(), queryMask.back()));
      assert( (queryMask.back().ui[0]&0x3)==0x3 );
      //      queryOligo.back().ui[0] |= nonIncorpBothNs;
      queryOligo.back().ui[0] |= (nonIncorpBothNs&queryMask.back().ui[0]);
      queryOligo.back().ui[1] |= (nonIncorpBothNs&queryMask.back().ui[1]);


      //         printWord(queryOligo.back().ui[1],prefixLength_);
      //    printWord(queryOligo.back().ui[0],16);
      //   cout << " - encoded " << endl;
      //    printWord(queryMask.back().ui[1],prefixLength_);
      //   printWord(queryMask.back().ui[0],16);
      //  cout << endl;

      if (secondN!=firstN+1)
      {
	for (i=0;i<secondN;i++) tempBuf2[i]=tempBuf1[i];
	for (i=secondN+1;i<OLIGO_LENGTH;i++) tempBuf2[i-1]=tempBuf1[i];
	tempBuf2[OLIGO_LENGTH-1]='n';
	//	 cout << "After removal of second internal N:" << endl;
	//	 cout << tempBuf2 << " - buffer" << endl;
	queryOligo.push_back( Oligo() );
	queryMask.push_back( Oligo() );
	queryOligoNum.push_back(oligoNum|seed_bits[seed_no]);

	assert( numNs
		==encodeOligo(tempBuf2, queryOligo.back(), queryMask.back()));
	assert( (queryMask.back().ui[0]&0x3)==0x3 );
	//	queryOligo.back().ui[0] |= nonIncorpSecondN;
	queryOligo.back().ui[0] |= (nonIncorpSecondN&queryMask.back().ui[0]);
	queryOligo.back().ui[1] |= (nonIncorpSecondN&queryMask.back().ui[1]);

	//		printWord(queryOligo.back().ui[1],prefixLength_);
	//		printWord(queryOligo.back().ui[0],16);
	//	 cout << " - encoded" << endl;
	//	printWord(queryMask.back().ui[1],prefixLength_);
	//	printWord(queryMask.back().ui[0],16);
	//	 cout << endl;

      } // ~if
      else
      {
	//	 cout << "Internal Ns are adjacent" << endl;
      }
   } // ~if
  } // ~if
#endif // ~ifdef ALLOW_N_TO_BE_DELETION

#ifndef DONT_SEARCH_REVERSE_STRAND
  // add reverse to all oligos in pile
  const int numToDo(queryOligo.size());
  for (int i(0);i<numToDo;++i)
  {
    queryOligo.push_back( Oligo() );
    queryMask.push_back( Oligo() );
    queryOligoNum.push_back(queryOligoNum[i]|(isReverseOligo|seed_bits[seed_no]));
    reverseOligo( queryOligo[i], queryOligo.back());

    // check for self-complementary sequences. These are rare but can
    // cause the oligo concerned to be spuriously flagged as a repeat
    if (queryOligo[i]!=queryOligo.back())
    { // not self-complementary, carry on
      reverseOligo( queryMask[i], queryMask.back());
      //      queryMask.back().ui[0] ^= ((uint)~0); // TBD
      // changed 03.06.03 TC
      // rhs of next line is supposed to contain all ones in its first
      // 2*prefixLength_ bits - does not work for prefixLength_=16
      //	queryMask.back().ui[1] ^= ~( ((uint)~0)<< (prefixLength_<<1) );
      //      queryMask.back().ui[1] ^=
      //	( ((uint)~0) >> ((maxBasesPerWord-prefixLength_)<<1) );
      // Changed to match new oligo format TC 05.08.03

      queryMask.back().ui[0] ^=
	( ((uint)~0) >> ((maxBasesPerWord-ElandConstants<OLIGO_LEN>::suffixLength)<<1) );
      queryMask.back().ui[1] ^=
	( ((uint)~0) >> ((maxBasesPerWord-ElandConstants<OLIGO_LEN>::prefixLength)<<1) );

    } // ~if
    else
    { // self-complementary
      queryOligo.pop_back();
      queryMask.pop_back();
      queryOligoNum.pop_back();
      //	cout << oligoNum << " is self-complementary" << endl;
    } // ~else

  } // ~for
#endif // ~ifndef DONT_SEARCH_REVERSE_STRAND

  return numNs;


}

};



// ======================================================================
// ========================== MULTI-SEED ================================
// ======================================================================
// now let's create a multi-seed query generator
template<int OLIGO_LEN> class MultiSeedQueryGenerator : public QueryGenerator<OLIGO_LEN>
{
 public:
  // c'tor
    MultiSeedQueryGenerator(bool single, const vector<int>& seedOffsets ):single_(single),seedOffsets_( seedOffsets )
    {
      // nothing to do
    }

  // d'tor
  ~MultiSeedQueryGenerator(void)
    {
      // nothing to do so far
    }

  // generate a set of query sequences from ASCII sequence data
  int operator()( const char* buf, const uint oligoNum,
		  vector<Oligo>& queryOligo,
		  vector<Oligo>& queryMask,
		  vector<OligoNumber>&  queryOligoNum,
		  vector<uint>& queryCnt )
    {
      queryOligo.clear();
      queryMask.clear();
      queryOligoNum.clear();
      queryCnt.clear();

      // by default we pack as many seeds into the read as possible if
      // single_ is true (ie the user specified only one seed per
      // read, set no_of_seeds accordingly
      int startIdx = 1;
      uint no_of_seeds = 4;

      if( single_ == true )
      {
          no_of_seeds = 1;
          startIdx = 0;
      }

      // in the multiseed phase, take overlapping seeds to be as
      // sensitive as possible; if we are in multiseed mode, we know
      // that the first 32bp either did not match (NM), or are a
      // hypermatch (255:255:255). Therefore, start the first seed at
      // base 16.
      //
      // set up an array with the starting positions for all the
      // seeds; consider those starting positions when inserting the
      // matches into StateMachine


      if( single_ == false )
      {
          assert( no_of_seeds == seedOffsets_.size() );
          no_of_seeds = 4;
      }

      int total_count = 0;

      // create a set of vectors for each vector to keep the
      // intermediate results
      //
      // the reason for not passing on queryOligo, queryMask or
      // queryOligoNum is that in the convert_ method we iterate over
      // the entire vector (create, for instance, the reverse
      // oligos,...)


      // cut the buffer into different seeds
      for(uint i=startIdx;i<no_of_seeds;i++ )
      {

          vector<Oligo> tmp_queryOligo;
          vector<Oligo> tmp_queryMask;
          vector<OligoNumber> tmp_queryOligoNum;
          char seed_buf[1024];
          strncpy( seed_buf,(buf+seedOffsets_[i]),1024 );

          total_count += QueryGenerator<OLIGO_LEN>::convert_( seed_buf,
                                     oligoNum,
                                     tmp_queryOligo,
                                     tmp_queryMask,
                                     tmp_queryOligoNum,(short)i );

          // append the temporary vectors
          queryOligo.insert( queryOligo.end(),tmp_queryOligo.begin(),tmp_queryOligo.end() );
          queryMask.insert( queryMask.end(),tmp_queryMask.begin(),tmp_queryMask.end() );
          queryOligoNum.insert( queryOligoNum.end(),tmp_queryOligoNum.begin(),tmp_queryOligoNum.end() );
          queryCnt.push_back( tmp_queryMask.size() );
      }

      // put everything together
      return total_count;
    }


 private:
  bool single_;
  vector<int> seedOffsets_;

};


} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_QUERY_GENERATOR_H
