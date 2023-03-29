/**
** MG, Feb 9 2023
**/
#ifndef TBFILE_DEFINED
#define TBFILE_DEFINED
#include <iostream>
#include <string>
#include <time.h>
#include<TNamed.h>
#include<TTimeStamp.h>
using namespace std;
// class to store file info

class TBFile : public TNamed
{
public:
  TBFile();
  TBFile(TString name);
  void update(TString name);
  
  //~TBFile();
  string created;
  string modified;
  TTimeStamp tsm; // modified time in seconds
  TTimeStamp tsc; // created time in seconds

  void
  print()
  {
    cout << " TBFile  " << this->GetTitle() << endl << " created " << created << " modified " << modified << endl;
  }

  ClassDef(TBFile, 2)
};
#endif
