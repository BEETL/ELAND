#include <cstdio>
#include <vector>

using namespace std;

enum ExtendedFileFieldName
{
  Machine=0,
  Read, 
  MatchCounter,
  Matches,
  NumberOfEntries // not a field, but if you stick it 
  // on the end of the list of fields it automatically gets right value
};

class ExtendedFileReaderImp;
// ExtendedFileReader
class ExtendedFileReader
{
 public:
  ExtendedFileReader(ExtendedFileReaderImp* p): p_(p) {} 

  ~ExtendedFileReader() {}

  // Returns pointer to ASCII sequence of next oligo, or null if at end
  bool getNextEntry( void ) ;

  // Rewind - next oligo read will be first in list
  void rewind ( void );

  const char* getMachine( void ) const; 
  const char* getRead ( void ) const;
  int getMatchCounter ( void ) const;
  const char* getXYZ ( void ) const;
  const char* getMatches ( void ) const;
 protected:
  ExtendedFileReaderImp* p_;
}; // ~class ExtendedFile

class ExtendedFileReaderImp
{
 public:

  virtual ~ExtendedFileReaderImp() {}

  // Returns pointer to ASCII sequence of next oligo, or null if at end
  virtual bool getNextEntry( void ) =0;

  // Rewind - next oligo read will be first in list
  virtual void rewind ( void ) =0;

  virtual const char* getMachine( void ) const=0; 
  virtual const char* getRead ( void ) const=0;
  virtual int getMatchCounter ( void ) const=0;
  virtual const char* getXYZ ( void ) const=0;
  virtual const char* getMatches ( void ) const=0;
 protected:
}; // ~class ExtendedFileReaderImp

class ExtendedFileReaderActual : public ExtendedFileReaderImp
{
 public:
  ExtendedFileReaderActual( const char* exportFileName );

  //  ExtendedFileReaderActual( const ExtendedFileReaderActual&  e )
  //    {
  //   pFile_ = e.pFile_;
  // }


  virtual ~ExtendedFileReaderActual();

  // Returns pointer to ASCII sequence of next oligo, or null if at end
  virtual bool getNextEntry( void );

  // Rewind - next oligo read will be first in list
  virtual void rewind ( void );

  virtual const char* getMachine( void ) const; 
  virtual const char* getRead ( void ) const;
  virtual int getMatchCounter ( void ) const;
  virtual const char* getXYZ ( void ) const;
  virtual const char* getMatches ( void ) const;

 protected:
  char buf_[120000];
  FILE* pFile_;
  vector<const char*> entry_;


}; // ~class ExtendedFileReaderActual


