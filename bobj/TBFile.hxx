/**
** MG, Feb 9 2023
**/
#ifndef TBFILE_DEFINED
#define TBFILE_DEFINED
#include <iostream>
#include <string>
#include<TNamed.h>
using namespace std;
// class to store file info

class TBFile : public TNamed
{
public:
  TBFile();
  TBFile(TString name);
  //~TBFile();
  string created;
  string modified;
  
  void print()
  {
    cout << " TBFile  " << this->GetTitle() << endl << " created " << created << " modified " << modified << endl;
  }

  ClassDef(TBFile, 1)
};
#endif
