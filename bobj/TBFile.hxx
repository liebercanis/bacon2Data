/**
** MG, Feb 9 2023
**/
#ifndef TBFILE_DEFINED
#define TBFILE_DEFINED
#include <iostream>
#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include<TNamed.h>
using namespace std;
// class to store file info

class TBFile : public TNamed
{
public:
  TBFile(TString name);
  //~TBFile();
  struct stat fileInfo;
  string created;
  string modified;
  
  void print()
  {
    cout << " TBFile  " << this->GetName() << " created " << created << " modified " << modified << endl;
  }

  ClassDef(TBFile, 1)
};
#endif
