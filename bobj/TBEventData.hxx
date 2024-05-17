/**
** MG, Feb 9 2023
**/
#ifndef TBEVENTDATA_DEFINED
#define TBEVENTDATA_DEFINED
#include <iostream>
#include <string>
#include <time.h>
#include<TNamed.h>
using namespace std;
// class to store file info

class TBEventData : public TNamed
{
public:
  TBEventData();
  void update(time_t time);

  time_t evtime; // event time;
  int sec;
  int min;
  int hour;
  int day;
  int mon;
  int year;
  int isdst;

  void print()
  {
    printf("TBEventData: sec %i min %i hour %i day %i mon %i year %i isdst %i\n", sec, min, hour, day, mon, year, isdst);
  }
  ClassDef(TBEventData, 4)
};

#endif
