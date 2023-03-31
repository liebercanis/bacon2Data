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
  TTimeStamp evtime; // eventtime in seconds
  void update(TTimeStamp tstamp){
    evtime = tstamp;
  }

  ClassDef(TBEventData, 1)
};
#endif
