/**
** MG, July 2020 
**/
#ifndef TBEVENT_DEFINED
#define TBEVENT_DEFINED
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <TNamed.h>
#include <TTree.h>

using namespace std;

// class to store info for the run

class TBEvent: public TNamed {
  public:
    TBEvent(TString runName = "run0");
    //~TBEvent();

    void clear();
    Long64_t event;
    double trigTime;
    Double_t energy;
    
   
    void print() {
      printf(" %s  event %lld  trigtime %f\n",this->GetName(), event, trigTime);
    }

    ClassDef(TBEvent,1)
};
#endif

