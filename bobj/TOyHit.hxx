/**
** MG, March 3 2025
**/
#ifndef TOYHIT_DEFINED
#define TOYHIT_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>

using namespace std;

// class to store pmt hit

class TOyHit : public TNamed
{
public:
  TOyHit();
  virtual ~TOyHit()
  {
  }

  void clear();
  // data elements
  int event;
  double time;
  double val;

  ClassDef(TOyHit, 1)
};
#endif
