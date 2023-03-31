/**
** MG, Feb 9 2023
**/
#ifndef TBEVENTDATA_DEFINED
#define TBEVENTDATA_DEFINED
#include <iostream>
#include <string>
#include <time.h>
#include<TNamed.h>
#include<TTimeStamp.h>
using namespace std;
// class to store file info

class TBEventData : public TNamed
{
public:
  TBEventData();
  TBEventData(TString name);
  void update(TString name);
  
  //~TBEventData();
  string created;
  string modified;
  TTimeStamp tsm; // modified time in seconds
  TTimeStamp tsc; // created time in seconds

  void
  print()
  {
    cout << " TBEventData  " << this->GetTitle() << endl << " created " << created << " modified " << modified << endl;
  }

  ClassDef(TBEventData, 2)
};
#endif
