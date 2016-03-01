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
 ** \file eland_ms/ElandThread.hh
 **
 ** \brief Part of ELAND
 **
 ** Part of ELAND
 **
 ** \author Tony Cox
 **/

#ifndef CASAVA_ELAND_MS_ELAND_THREAD_H
#define CASAVA_ELAND_MS_ELAND_THREAD_H

#include "MatchTable.hh"
#include "QueryGenerator.hh"
#include "Hasher.hh"
#include "OligoHashTable.hh"
#include "ElandDefines.hh"

namespace casava
{
namespace eland_ms
{
// code to scan multiple chromosomes using threads
#if (NUM_THREADS>1)

template <int PASS> void* scan( void* pv );
//void* scan( void* pv );

// ThreadInfo: contains data specific to an individual thread, plus
// static functions & data to maintain a count of the total number of threads

template<int PASS> struct ThreadInfo
{
  ThreadInfo( OligoHashTable<PASS>* pHash,
	      FileReader* pFile,
	      MatchPosition currentBlock, const char* name ) :
    pHash_(pHash), pFile_(pFile),
    currentBlock_(currentBlock), name_(name) {}

  ~ThreadInfo()
  {
    delete pFile_;
  } // ~dtor

  OligoHashTable<PASS>* pHash_;
  FileReader* pFile_;
  MatchPosition currentBlock_;
  const char* name_;

  static int threadCount_;
  static pthread_mutex_t threadCountLock_;
  static pthread_cond_t threadFreeSignal_;

  // start thread, waiting for others to complete if necessary
  static void startThread( ThreadInfo* pt )
  {
    pthread_t p;
    pthread_mutex_lock(&threadCountLock_);
    while (threadCount_==NUM_THREADS)
    {
      cout << "Waiting for thread to become free..." << endl;
      pthread_cond_wait(&threadFreeSignal_, &threadCountLock_);
    } // ~while
    cout << "OK, got a free thread..." << endl;
    ++threadCount_;
    pthread_mutex_unlock(&threadCountLock_);
    pthread_create(&p, NULL, scan<PASS>,  (void*) pt);
    pthread_detach(p);

  } // ~startThread

  // decrement thread count, and signal to say a thread is free
  static void freeThread( void )
  {
    cout << "Freeing thread..."  << endl;
    pthread_mutex_lock(&threadCountLock_);
    --threadCount_;

    pthread_cond_signal(&threadFreeSignal_);

    pthread_mutex_unlock(&threadCountLock_);

  } // ~freeThread

  // wait for all outstanding threads to finish
  static void waitForEnd( void )
  {
    pthread_mutex_lock(&threadCountLock_);
    while (threadCount_>0)
    {
      cout << "Waiting for " << threadCount_ << " threads to become free..." << endl;
      pthread_cond_wait(&threadFreeSignal_, &threadCountLock_);
    } // ~while
    cout << "All threads done, OK to finish" << endl;
    pthread_mutex_unlock(&threadCountLock_);
  } // ~waitForEnd

}; // ~struct ThreadInfo

template<int PASS>
int ThreadInfo<PASS>::threadCount_=0;

template<int PASS>
pthread_mutex_t ThreadInfo<PASS>::threadCountLock_ = PTHREAD_MUTEX_INITIALIZER;

template<int PASS>
pthread_cond_t ThreadInfo<PASS>::threadFreeSignal_ = PTHREAD_COND_INITIALIZER;

template <int PASS> void* scan( void* pv )
{
  Timer t;
  ThreadInfo<PASS>* pt((ThreadInfo<PASS>*)pv);
  cerr << "Starting thread to scan " << pt->name_ << ": " << t << endl;
  //  currentBlock = hashTable.scan(thisFile, currentBlock);
  (pt->pHash_)->scan(*(pt->pFile_), pt->currentBlock_ );
  cerr << "Finished thread to scan " << pt->name_ << ": " << t << endl;
  ThreadInfo<PASS>::freeThread();
  delete pt;
  return NULL;
} // ~scan

// scanAll: complete a single pass through all the chromosomes
template< int PASS > void  scanAll
(
 OligoSource* pOligos,
 const string& directoryName,
 const vector<string>& chromNames,
 vector<MatchPosition>& blockStarts,
 OligoHashTable<PASS>& hashTable,
 Timer& timer,
 bool singleseed=true/*,
 bool firstRun=true*/
)
{
  // reset block counter at start of each pass
  MatchPosition currentBlock(blockSize);

  string fullChromName;

  ThreadInfo<PASS>* pThread;

  if (PASS==0) {
    blockStarts.push_back(currentBlock);
  }
  else {
    assert(blockStarts[1]==currentBlock);
  }

  cerr << "About to build hash tables for pass " << PASS << ": " << timer << endl;

  if (hashTable.buildTable( *pOligos,singleseed )==false)
  {
    cerr << "No oligos to hash, returning" << endl;
    return;
  }

  cerr << "Built hash tables: " << timer << endl;

  for (int j(1);j<chromNames.size();j++)
  {
    //    if (PASS==0) blockStarts.push_back(currentBlock);

    fullChromName = directoryName+chromNames[j];
    cerr << "Scanning file " << fullChromName << ": " << timer << endl;

    pThread = new ThreadInfo<PASS>(&hashTable, NULL, currentBlock, chromNames[j].c_str());

    cerr << "Starting block: " << (currentBlock>>blockShift) << endl;

    pThread->pFile_=new FileReader(fullChromName.c_str());

    cerr << "Last valid base in file: " << pThread->pFile_->getLastValidBase() << endl;

    currentBlock+=(((pThread->pFile_->getLastValidBase()>>blockShift)+1) <<blockShift);

    if (PASS==0) {
      blockStarts.push_back(currentBlock);
    }
    else {
      assert(blockStarts[j+1]==currentBlock);
    }

    ThreadInfo<PASS>::startThread(pThread);

    cerr << "Finishing block: " << (currentBlock>>blockShift) << endl;

  }

  ThreadInfo<PASS>::waitForEnd();
    //  if (PASS==0) blockStarts.push_back(currentBlock);
  assert(blockStarts.size()==chromNames.size()+1);
} // ~scanAll

// single threaded scanAll
#else

// scanAll: complete a single pass through all the chromosomes
template< int PASS, int OLIGO_LEN> void  scanAll
(OligoSource* pOligos,
 const string& directoryName,
 const vector<string>& chromNames,
 vector<MatchPosition>& blockStarts,
 OligoHashTable<PASS, OLIGO_LEN>& hashTable,
 Timer& timer,
 bool singleseed=true/*,
 bool firstRun=true*/
)
{
  // reset block counter at start of each pass
  MatchPosition currentBlock(blockSize);

  string fullChromName;

  //  OligoHashTable hashTable(oligoLength, numBits, part1, part2, scoreTable, results);
    // set base to extract
  //  hashTable.setPass( i );

  cerr << "About to build hash tables for pass " << PASS << ": " << timer << endl;

  if (hashTable.buildTable( *pOligos,singleseed )==false)
  {
    cerr << "No oligos to hash, returning" << endl;
    return;
  }
  //  hashTable.buildTable( *pOligos );

  cerr << "Built hash tables: " << timer << endl;

  for (uint j(1);j<chromNames.size();j++)
  {
    if (PASS==0) blockStarts.push_back(currentBlock);


    fullChromName = directoryName+chromNames[j];
    cerr << "Scanning file " << fullChromName << ": " << timer << endl;

    cerr << "Starting block: " << (currentBlock>>24) << endl;
    FileReader thisFile(fullChromName.c_str());

    currentBlock = hashTable.scan(thisFile, currentBlock);
    cerr << "Finishing block: " << (currentBlock>>24) << endl;

    cerr << "... done " << timer << endl;
  } // ~for j
  if (PASS==0) blockStarts.push_back(currentBlock);
  assert(blockStarts.size()==chromNames.size()+1);
}

#endif


} //namespace eland_ms
} //namespace casava

#endif // CASAVA_ELAND_MS_ELAND_THREAD_H
