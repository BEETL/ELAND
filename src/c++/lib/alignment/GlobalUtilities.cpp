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
 ** \file GlobalUtilities.cpp
 **
 ** \brief Commonly used functions and definitions
 **
 ** \author Tony Cox
 **/

#include <iostream>
#include <ctime>
#include <dirent.h>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "alignment/GlobalUtilities.hh"
#include "eland_ms/ElandConstants.hh"

void printWord(Word w, int l)
{
    for (int i(2* (l -1)); i>=0; i-=2 )
    cout << baseNames[ (w >> i)&0x3 ];
}  // ~void printWord( Word w, int l )

char reverseComp(const char c)
{
    return (c == 'A') ? 'T' : ((c == 'T') ? 'A' : ((c == 'C') ? 'G' : ((c
            == 'G') ? 'C' : '!')));
} // ~reverseComp;

//RP: thread-unsafe, commented-out
//const char* getPrefixString(Word w, uint prefixLength)
//{
//    static char buf[10];
//    char* pNext(buf + prefixLength);
//    *pNext = '\0';
//    do
//    {
//        *(--pNext) = baseNames[w & 0x3];
//        w >>= numBitsPerBase;
//    } while (pNext != buf);
//
//    return buf;
//} // ~getPrefixString( Word w, uint prefixLength )

/*****************************************************************************/

Timer::Timer(void)
{
    if (gettimeofday(&lastTime_, NULL) != 0)
        exit(-1);
    if (getrusage(RUSAGE_SELF, &lastUsage_) != 0)
        exit(-1);
} // ~Timer::Timer( void )


ostream& Timer::print(ostream& os)
{
    if (gettimeofday(&thisTime_, NULL) != 0)
        exit(-1);
    if (getrusage(RUSAGE_SELF, &thisUsage_) != 0)
        exit(-1);

    static double elapsedActual, elapsedUser, elapsedSystem;

    elapsedUser = thisUsage_.ru_utime.tv_sec - lastUsage_.ru_utime.tv_sec
            + (thisUsage_.ru_utime.tv_usec
                    - (double) lastUsage_.ru_utime.tv_usec) / 1000000;

    elapsedSystem = thisUsage_.ru_stime.tv_sec - lastUsage_.ru_stime.tv_sec
            + (thisUsage_.ru_stime.tv_usec
                    - (double) lastUsage_.ru_stime.tv_usec) / 1000000;

    elapsedActual = thisTime_.tv_sec - lastTime_.tv_sec + (thisTime_.tv_usec
            - (double) lastTime_.tv_usec) / 1000000;

    os << "User: " << elapsedUser << "s System: " << elapsedSystem
            << "s Actual: " << elapsedActual << "s Efficiency: "
            << ((elapsedActual == 0) ? 0 : ((elapsedUser + elapsedSystem)
                    * 100.0 / elapsedActual)) << '\045'; // only way to print '%' liked by both Intel and gcc!
    lastTime_ = thisTime_;
    lastUsage_ = thisUsage_;
    return os;
} // ~ostream& Timer::print( ostream& os )

ostream& operator<<(ostream& os, Timer& timer)
{
    return timer.print(os);
} // ~ostream& operator<<( ostream& os, Timer& timer )

const char* Timer::timeNow(void) const
{
    time_t tt(time( NULL));
    return ctime(&tt);
} // ~const char* Timer::timeNow( void ) const


/*****************************************************************************/
// ExpandedTranslationTable function definitions

ExpandedTranslationTable::ExpandedTranslationTable(int oligoLength) :
    prefixLength_((oligoLength - 16) >> 1), isOddLength_((oligoLength & 0x1)
            != 0), reverseShift1_(2* ( (2*maxBasesPerWord)-oligoLength)),
    reverseShift2_(numBitsPerWord-reverseShift1_)
    {
        for (int i(0);i<numPossibleChars;i++)
        for (int j(0);j<numPossibleChars;j++)
        {
            if ((whichBase[i]==nv)||(whichBase[j]==nv))
            t_[(i<<8)+j]=nv;
            else
            t_[(i<<8)+j]=((whichBase[j]<<numBitsPerBase)|whichBase[i]);
        } // ~for
        } // ~ExpandedTranslationTable::ctor

void ExpandedTranslationTable::translate(const char* buf, Oligo& ol, Oligo& rc) const
{
    ol.ui[0] = 0;
    ol.ui[1] = 0;
    if (isOddLength_)
    {
        ol.ui[1] |= whichBase[(uint)(*buf++)];
    } // ~if isOddLength_

    ushort* ps((ushort*) buf);

    for (int i(0); i < prefixLength_; ++i)
    {
        ol.ui[1] <<= (2* numBitsPerBase );
        ol.ui[1] |= t_[*(ps++)];
    }

    //  ol.ui[1]<<=(2*numBitsPerBase); // 0
    ol.ui[0] |= t_[*(ps++)];

    ol.ui[0] <<= (2* numBitsPerBase ); // 1
    ol.ui[0] |= t_[*(ps++)];

    ol.ui[0] <<= (2* numBitsPerBase ); // 2
    ol.ui[0] |= t_[*(ps++)];

    ol.ui[0] <<= (2* numBitsPerBase ); // 3
    ol.ui[0] |= t_[*(ps++)];

    ol.ui[0] <<= (2* numBitsPerBase ); // 4
    ol.ui[0] |= t_[*(ps++)];

    ol.ui[0] <<= (2* numBitsPerBase ); // 5
    ol.ui[0] |= t_[*(ps++)];

    ol.ui[0] <<= (2* numBitsPerBase ); // 6
    ol.ui[0] |= t_[*(ps++)];

    ol.ui[0] <<= (2* numBitsPerBase ); // 7
    ol.ui[0] |= t_[*(ps++)];

    rc.uc[0] = reverseChar[ol.uc[7]];
    rc.uc[1] = reverseChar[ol.uc[6]];
    rc.uc[2] = reverseChar[ol.uc[5]];
    rc.uc[3] = reverseChar[ol.uc[4]];
    rc.uc[4] = reverseChar[ol.uc[3]];
    rc.uc[5] = reverseChar[ol.uc[2]];
    rc.uc[6] = reverseChar[ol.uc[1]];
    rc.uc[7] = reverseChar[ol.uc[0]];

    rc.ui[0] >>= reverseShift1_;
    rc.ui[0] |= (rc.ui[1] << reverseShift2_);
    rc.ui[1] >>= reverseShift1_;

} // ~ExpandedTranslationTable::translate( char* buf, Word ol[2] )


/*****************************************************************************/

FileReader::FileReader(const char* fileName) :
    seqFileName_(fileName), vldFileName_(fileName)
{
    seqFileName_ += (string) ".2bpb";
    vldFileName_ += (string) ".vld";

    int fd;

    char* pChar;

    // memory map sequence file
    fd = open(seqFileName_.c_str(), O_RDONLY, S_IRUSR);
    if (fd == -1)
    {
        cerr << "Error in FileReader: could not open file " << seqFileName_
                << endl;
        exit(1);
    } // ~if
    seqFileSize_ = lseek(fd, 0, SEEK_END);
    cerr << "FileReader: opened file " << seqFileName_ << " of size "
            << seqFileSize_ << " bytes." << endl;

    pMapSeq_ = mmap(0, seqFileSize_, PROT_READ, MAP_SHARED, fd, 0);
    if (pMapSeq_ == MAP_FAILED)
    {
        cerr << "Error in FileReader: could not memory map file "
                << seqFileName_ << ": " << strerror(errno) << endl;
        exit(1);
    } // ~if
    close(fd);

    // memory map validity file
    fd = open(vldFileName_.c_str(), O_RDONLY, S_IRUSR);
    if (fd == -1)
    {
        cerr << "Error in FileReader: could not open file " << vldFileName_
                << endl;
        exit(1);
    } // ~if
    vldFileSize_ = lseek(fd, 0, SEEK_END);
    cerr << "squash: opened file " << vldFileName_ << " of size "
            << vldFileSize_ << " bytes." << endl;

    pMapVld_ = mmap(0, vldFileSize_, PROT_READ, MAP_SHARED, fd, 0);
    if (pMapVld_ == MAP_FAILED)
    {
        cerr << "Error in unsquash: could not memory map file " << vldFileName_
                << ": " << strerror(errno) << endl;
        exit(1);
    } // ~if
    close(fd);

    // set up pointers for Word

    pChar = (char*) pMapSeq_;
    pStart_ = (Word*) pMapSeq_;
    pEnd_ = (Word*) (pChar + seqFileSize_);

    // set up last pointer for valid region

    pChar = (char*) pMapVld_;
    pLastValid_ = (ValidRegion*) (pChar + vldFileSize_);

    cerr << "Comment: ";

    // print out the comment line ...
    while (*pChar != '\n')
    {
        cerr.put(*pChar);
        pChar++;
    } // ~while

    cerr.put('\n');
    // ... and the '\n' afterwards
    //  cout.put(*(pChar++));
    pChar++;
    //  cout << endl; // %%%%%

    // set up
    pValid_ = (ValidRegion*) pChar;
    //  assert( pValid_!=pLastValid_ );
    //  getNextValidRegion();

    cerr << "finished making FileReader" << endl; // %%%%
} // ~FileReader::FileReader

FileReader::~FileReader(void)
{
    cerr << "FileReader: unmapping " << seqFileSize_ << " bytes of memory"
            << endl;
    munmap(pMapSeq_, seqFileSize_);
    cerr << "FileReader: unmapping " << vldFileSize_ << " bytes of memory"
            << endl;
    munmap(pMapVld_, vldFileSize_);

} // ~FileReader::~FileReader


/*****************************************************************************/
// getOligoSource: given the name of a file of oligo data, returns a pointer
// to an instance of the appropriate subclass of OligoSource.
OligoSource* OligoSource::getOligoSource(const char* fileName)
{
    DIR* pDir;
    if ((pDir = opendir(fileName)) != NULL)
    {
        closedir(pDir);
        cerr << "getOligoSource: assuming " << fileName
                << " is a directory of oligo files" << endl;
        return new OligoSourceDirectory(fileName);
    } // ~if
    else
    {
        // makes more sense to delegate rest to getOligoSource until (if) we want
        // to get oligos from non-file sources such as database clients
        return OligoSourceFile::getOligoSourceFile(fileName);
    } // ~else
} // ~getOligoSource

OligoSourceFasta::OligoSourceFasta( const char* oligoFileName ) : OligoSourceFile(oligoFileName)
{
    curSeq_ = 1;
    skippedSequences_ = 0;
}
//    OligoSourceFile(oligoFileName) { curSeq_ = 1; mask_.resize(134217728,true); skippedSequences_ = 0; }

/*****************************************************************************/
// getOligoSource: given the name of a file of oligo data, returns a pointer
// to an instance of the appropriate subclass of OligoSourceFile.
// Examines first character of file:
// If a '>', assumes a fasta file and returns an OligoSourceFasta
// If valid sequence character, assumes raw sequence and return OligoSourceRaw
// else exit with an error
OligoSourceFile* OligoSourceFile::getOligoSourceFile(const char* fileName)
{
    //  ifstream oligoFile(fileName);
    // if (oligoFile.fail())
    //  {
    //   cerr << "Error in getOligoSourceFile: could not open file "
    //	 << fileName << endl;
    //   exit (1);
    //  } // ~if
    //  char firstChar(oligoFile.get());
    //  oligoFile.close();
    FILE* pFile = fopen(fileName, "r");
    if (pFile == NULL)
    {
        cerr << "Error in getOligoSourceFile: could not open file " << fileName
                << endl;
        exit(1);
    } // ~if
    int outInt;
    if ((outInt = fgetc(pFile)) == EOF)
    {
        cerr
                << "Warning in getOligoSourceFile: could not read first char from file '"
                << fileName << "'. Assuming it is Goat formatted." << endl;
        return new OligoSourceGoat(fileName);
        // exit(0); // Allow make pipeline to carry on... Zero exit code.
    } // ~if

    fclose(pFile);
    char firstChar((char) outInt);

    std::ifstream istrm(fileName);
    casava::common::Sequence sequence;

    istrm.clear();
    // DAP-362: prevent the input operator from reading the whole file
    const unsigned int BUFFER_SIZE = 4096;
    std::string buffer;
    buffer.reserve(BUFFER_SIZE);
    char c = 0;
    while (('\n' != c) && istrm && istrm.get(c))
    {
        buffer.push_back(c);
    }
    std::istringstream bufferIstream(buffer);
    bufferIstream >> sequence;

    if (!bufferIstream.fail())
    {
//        cerr << "getOligoSourceFile: detected Goat format oligo file" << endl;
        return new OligoSourceGoat(fileName);
    }
    if (firstChar == '>')
    {
        cerr << "getOligoSourceFile: detected fasta format oligo file" << endl;
        return new OligoSourceFasta(fileName);
    } // ~if
    else if ((whichBase[(uint) firstChar] != nv) || isBlank(firstChar))
    {
        cerr << "getOligoSourceFile: detected raw sequence format oligo file"
                << endl;
        return new OligoSourceRaw(fileName);
    } // ~else if
    else if ((firstChar == ' ') || (firstChar == '-'))
    {
        cerr << "getOligoSourceFile: detected quality value format oligo file"
                << endl;
        return new OligoSourceScore(fileName);
    }
    else
    {
        cerr << "Error in getOligoSourceFile: "
                << "could not recognize format of file " << fileName << endl;
        exit(1);
    }

    return NULL; // compiler too stupid to realize can never get here
} // ~OligoSourceFile* getOligoSourceFile( const char* fileName )

/*****************************************************************************/
// **
// ** Function definitions for OligoSourceGoat
// **
/*****************************************************************************/

OligoSourceGoat::OligoSourceGoat(const char* oligoFileName) :
    OligoSourceFile(oligoFileName)
{
    qseq_file_.open(oligoFileName);

    if (!qseq_file_.is_open())
    {
        cerr << "Error in OligoSourceGoat: could not open file "
                << oligoFileName << endl;
        exit(1);
    }
}

/*****************************************************************************/

OligoSourceGoat::~OligoSourceGoat()
{
    qseq_file_.close();
}

/*****************************************************************************/
// Returns reference to next Sequence (supersedes getNextOligo).
// isValid will be false if there are no sequences left.
const casava::common::Sequence& OligoSourceGoat::getNextSequenceSelect(bool& isValid,
                                                                       const bool,
                                                                       const bool)
{
    if (qseq_file_.eof())
    {
        isValid = sequenceIsValid_ = false;
        return sequence_;
    }

    // Force discovery of EOF if not already noticed.
    qseq_file_.peek();

    if (qseq_file_.eof())
    {
        isValid = sequenceIsValid_ = false;
        return sequence_;
    }

    qseq_file_.clear();
    qseq_file_ >> sequence_;

    if (qseq_file_.fail())
    {
        cerr << "Error: Failed to read qseq line." << endl;
        exit(1);
    }

    isValid = sequenceIsValid_ = true;

    sprintf(nameBuf_, ">%d-%d-%d-%d",
            sequence_.getSpot().getTile().getLaneNumber(),
            sequence_.getSpot().getTile().getTileNumber(),
            sequence_.getSpot().getX(), sequence_.getSpot().getY());

    return sequence_;
}

/*****************************************************************************/
// Returns pointer to ASCII name of last oligo read
const char* OligoSourceGoat::getLastName(void)
{
    return nameBuf_;
}

/*****************************************************************************/

void OligoSourceGoat::rewind(void)
{
    OligoSourceFile::rewind();
    qseq_file_.seekg(0);
}

/*****************************************************************************/
// **
// ** Function definitions for OligoSourceScore
// **

// Returns reference to next Sequence (supersedes getNextOligo).
// isValid will be false if there are no sequences left.
const casava::common::Sequence& OligoSourceScore::getNextSequenceSelect(bool& isValid,
                                                                        const bool,
                                                                        const bool)
{
    isValid = sequenceIsValid_ = (fgets(oligoBuf_, maxLineLength, pFile_)
            != NULL);

    if (!sequenceIsValid_)
    {
        return sequence_;
    }

    const char* p(oligoBuf_);
    int i(0);
    int temp[numDifferentBases];
    // Fill in score values
    while (1)
    {
        if (sscanf(p, "%d%d%d%d\t", &temp[0], &temp[1], &temp[2], &temp[3])
                == EOF)
            break;
        for (int j(0); j < numDifferentBases; j++)
            scoreTable_[j][i] = (BaseScore) temp[j];
        // TBD check for overflow in above loop
        //      cout << i << ": " << qualityStore_[0][i] << " "<< qualityStore_[1][i] << " "<< qualityStore_[2][i] << " "<< qualityStore_[3][i] << endl;
        i++;

        if ((p = strchr(p, '\t')) == NULL)
            break;
        // scan past the deblock cycle stuff
        // no ... don't!!
        //      if ((p=strchr(++p,'\t'))==NULL) break;
        ++p;
        //      if (p!=NULL) { p++; p=strchr(p,'\t'); cout << "scan" << endl;}
        //  if (p!=NULL) p++;
    } // ~while

    // Call consensus based on highest score - TBD handle "dead heats"
    BaseScore maxQual, thisQual;
    int maxPos;
    for (int j(0); j < i; j++)
    {
        maxQual = scoreTable_[0][j];
        maxPos = 0;
        for (int k(1); k < numDifferentBases; k++)
        {
            thisQual = scoreTable_[k][j];
            if (thisQual > maxQual)
            {
                maxPos = k;
                maxQual = thisQual;
            }
            else if (thisQual == maxQual)
            {
                maxPos = -1;
            }
        } // ~for k
        oligoBuf_[j] = ((maxPos == -1) ? '.' : baseNames[maxPos]);
    } // ~for j
#ifdef OLD
    for (int j(0);j<i;j++)
    {
        maxQual=scoreTable_[0][j];
        maxPos=0;
        for (int k(1); k<numDifferentBases; k++)
        {
            if (scoreTable_[k][j]>maxQual)
            {
                maxQual=scoreTable_[k][j];
                maxPos=k;
            } // ~if
        } // ~for k
        oligoBuf_[j]=baseNames[maxPos];
    } // ~for j
#endif
    oligoBuf_[i] = '\0';
    sequence_.setData(oligoBuf_);

    return sequence_;
} // ~OligoSourceScore::getNextOligo( void )


/*****************************************************************************/
// **
// ** Function definitions for OligoSourceDirectory
// **

OligoSourceDirectory::OligoSourceDirectory(const char* dirName) :
    pSource_(NULL)
{
    DIR* pDir;
    if (!(pDir = opendir(dirName)))
    {
        cerr << "Error: failed to open directory " << dirName << endl;
        exit(1);
    } // ~if

    dirent* dirEntry;
    while ((dirEntry = readdir(pDir)))
    {
        if (dirEntry->d_name[0] != '.')
        {
            fileNames_.push_back((string) dirName + (string) "/"
                    + (string) dirEntry->d_name);
        } // ~if
    } // ~while
    closedir(pDir);
    if (fileNames_.size() == 0)
    {
        cerr << "Error: directory " << dirName << " is empty" << endl;
        exit(1);
    } // ~if

    pName_ = fileNames_.begin();
    pSource_ = getOligoSource(pName_->c_str());
} // ~ctor

// Returns reference to next Sequence (supersedes getNextOligo).
// isValid will be false if there are no sequences left.
const casava::common::Sequence&
OligoSourceDirectory::getNextSequenceSelect(bool& isValid,
                                            const bool isProvideHeaders,
                                            const bool isProvideQualities)
{
    if (!(isValid = (pName_ != fileNames_.end())))
    {
        return dummySequence_;
    }

    const casava::common::Sequence& sequence(pSource_->getNextSequenceSelect(isValid,
                                                                             isProvideHeaders,
                                                                             isProvideQualities));

    if (isValid)
    {
        return sequence;
    }
    else
    { // No more in that file - more files?
        if (++pName_ == fileNames_.end())
        {
            return dummySequence_;
        }

        delete pSource_;
        pSource_ = getOligoSource(pName_->c_str());
        // just calling pSource_->getNextOligo here mucks things up
        // if for any reason the source is empty
        return getNextSequenceSelect(isValid,
                                     isProvideHeaders,
                                     isProvideQualities);
    } // ~else
}

// Returns reference to last Sequence fetched (supersedes getLastOligo).
// isValid will be false if no sequences have been successfully read.
const casava::common::Sequence&
OligoSourceDirectory::getLastSequence(bool& isValid) const
{
    return ((isValid = (pName_ != fileNames_.end())) ? pSource_->getLastSequence(
            isValid)
            : dummySequence_);
}

// Returns pointer to ASCII sequence of next oligo, or null if at end
const char* OligoSourceDirectory::getNextOligo(void)
{
    bool isValid(false);
    const casava::common::Sequence& sequence(getNextSequence(isValid));
    return (isValid ? sequence.getData().c_str() : NULL);
} // ~OligoSourceDirectory::getNextOligo( void )

// Returns pointer to ASCII name of last oligo read
const char* OligoSourceDirectory::getLastName(void)
{
    if (pName_ == fileNames_.end())
        return NULL;
    return pSource_->getLastName();
} // ~OligoSourceDirectory::getLastName( void )

// Rewind - next oligo read will be first in list
void OligoSourceDirectory::rewind(void)
{
    delete pSource_;
    pName_ = fileNames_.begin();
    pSource_ = getOligoSource(pName_->c_str());
} // ~OligoSourceDirectory::rewind( void )


/*****************************************************************************/
// ScoreSource subclass function definitions

BaseScore OligoSourceScore::getScore(uint baseNum, uint cycle) const
{
    assert(baseNum < numDifferentBases);
    assert(cycle < maxSeqSize);
    return scoreTable_[baseNum][cycle];
} // ~OligoSourceScore::getScore( uint baseNum, uint cycle ) const


BaseScore OligoSourceScore::getScore(char base, uint cycle) const
{
    return getScore(baseCodes[(uint) base], cycle);
} // ~OligoSourceScore::getScore( char base, uint cycle ) const


/*****************************************************************************/

BaseScore ScoreSourceBasic::getScore(char base, uint cycle) const
{
    assert(cycle < maxSeqSize);
    const char b(oligos_.getLastOligo()[cycle]);
    if (isBlank(b))
        return scoreBlank_;
    else if (tolower(base) == tolower(b))
        return scoreMatch_;
    else
        return scoreMismatch_;
} // ~getScore( char base, uint cycle ) const

BaseScore ScoreSourceBasic::getScore(uint baseNum, uint cycle) const
{
    assert(baseNum < numDifferentBases);
    return ScoreSourceBasic::getScore(baseChars[baseNum], cycle);
} // ~getScore( uint baseNum, uint cycle ) const

/*****************************************************************************/

ScoreSourceCycle::ScoreSourceCycle(const OligoSource& oligos,
        const char* scoreFile) :
    oligos_(oligos)
{

    char lineBuf[maxLineLength];
    int scores[numDifferentBases];
    const BaseScore scoreNotSet(-20000);
    char thisChar;
    cerr << "Reading score information from file " << scoreFile << endl;

    for (uint i(0); i < maxSeqSize; i++)
        for (uint j(0); j < numDifferentBases + 1; j++)
            for (uint k(0); k < numDifferentBases; k++)
                scoreTable_[i][j][k] = scoreNotSet;

    FILE* pScore;
    int i;
    if ((pScore = fopen(scoreFile, "r")) == NULL)
    {
        cerr << "Error: could not open file " << scoreFile << endl;
        exit(1);
    } // ~if

    while (fgets(lineBuf, maxLineLength, pScore) != NULL)
    {
        if (lineBuf[0] != '>')
            continue;
        //      sscanf( lineBuf, ">%i\t%i\t%i",
        //      &i, &match, &mismatch);
        sscanf(lineBuf, ">%i\t%c\t%i\t%i\t%i\t%i", &i, &thisChar, scores,
                scores + 1, scores + 2, scores + 3);
        assert((i >= 1) && (i <= maxSeqSize));
        assert(baseCodes[(uint) thisChar] != nc);

        for (int j(0); j < numDifferentBases; j++)
            scoreTable_[i - 1][baseCodes[(uint) thisChar]][j] = (int) scores[j];
        //      scoreMatch[i-1]=match;
        //   scoreMismatch[i-1]=mismatch;
    } // ~while

    // check we have found all score table values we need
    //    for (uint i(0); i<seqSize; i++)
    //  for (uint j(0); j<numReadableBases; j++)
    //	for (uint k(0); k<numDifferentBases; k++)
    //	  assert(scoreTable[i][j][k]!=scoreNotSet);

    fclose(pScore);

} //~ScoreSourceCycle::ctor


BaseScore ScoreSourceCycle::getScore(uint baseNum, uint cycle) const
{
    assert(baseNum < numDifferentBases);
    assert(cycle < maxSeqSize);
    return scoreTable_[cycle][baseCodes[(uint) oligos_.getLastOligo()[cycle]]][baseNum];
} // ~ScoreSourceCycle::getScore( uint baseNum, uint cycle ) const


BaseScore ScoreSourceCycle::getScore(char base, uint cycle) const
{
    return getScore(baseCodes[(uint) base], cycle);
} // ~ScoreSourceCycle::score( char base, int cycle ) const


const vector<bool> expandUseBases(const std::string &useBases)
{
    std::string str(useBases);
    size_t found=0;
    while ( found < (str.length() - 1) &&
           (found = str.find_first_of("0123456789",found+1)) != std::string::npos )
    {
        std::string s = str.substr(found);
        size_t len = s.find_first_not_of("0123456789");
        unsigned int r = boost::lexical_cast<unsigned int>(s.substr(0,len));
        if (r)
        {
            if (len != std::string::npos)  len++;
            str.replace(found-1,len,r,str[found-1]);
        } else {
            cerr << (boost::format("Wrong format in Use Bases: Cannot repeat '%c' zero times") % str[found-1]).str() << endl;
            exit(1);
        }
    }
    vector<bool> vb;
    for (string::const_iterator it=str.begin(); it != str.end(); ++it)
    {
        unsigned char c = static_cast<unsigned char>(*it);
        vb.push_back(c == 'Y' || c == 'y');
    }
    return vb;
}
// End of file GlobalUtilities.cpp


