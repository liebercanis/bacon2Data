/**
** MG, July 2020
**/
#ifndef TBRUN_DEFINED
#define TBRUN_DEFINED
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <TNamed.h>
#include <TTree.h>
#include "TBEvent.hxx"
#include "TDet.hxx"

using namespace std;

// class to store info for the run

class TBRun : public TNamed
{
public:
  TBRun(TString tag = TString("none"));
  //TBRun(TTree *btreeInit, TString runName = "run0");
  ~TBRun(){
    detList.clear();
  }
  TTree *btree;
  //TBEvent *bevent;

  void clear();
  // data elements
  vector<TDet *> detList;

//void dumpEvent(ofstream &dumpFile);

Long64_t getNevents() {
  return btree->GetEntriesFast();
}

  bool addDet(int ichan, int ilevel = -1)
  {
    // check that it is not already in list
    for (unsigned i = 0; i < detList.size(); ++ i) {
      if (ichan == detList[i]->channel)
      return false;
    }
    TDet *det = new TDet(ichan, ilevel);
    detList.push_back(det);
    btree->Branch(det->GetName(), det);
    return true;
  }

  TDet *getDet(int ichan)
  {
    TDet  *rev = NULL;
    for (unsigned i = 0; i < detList.size(); ++i)
    {
      if (detList[i]->channel == ichan)
      {
        rev = detList[i];
        break;
      }
    }
    return rev;
  }

  Int_t fill()
  {
    return btree->Fill();
  }

  void print()
  {
    //bevent->print();
    printf(" TBRun has %llu events %lu dets \n", btree->GetEntriesFast(), detList.size());
    for (unsigned i = 0; i < detList.size(); ++i)
      printf(" %s %lu \n", detList[i]->GetName(), detList[i]->hits.size());
    btree->GetListOfBranches()->ls();
  }

  ClassDef(TBRun, 3)
};
#endif
