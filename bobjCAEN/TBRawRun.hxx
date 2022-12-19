/**
** MG, July 2020
**/
#ifndef TBRAWRUN_DEFINED
#define TBRAWRUN_DEFINED
#include <iostream>
#include <string>
#include <map>
#include <TNamed.h>
#include <TTree.h>
#include "TBRawEvent.hxx"

using namespace std;

// class to store info for the run

class TBRawRun : public TNamed
{
public:
  TBRawRun();
  TBRawRun(TString runName = "run0");
  //~TBRawRun();

  void clear();
  TTree *btree;
  // TBEvent *bevent;

  vector<TBRawEvent *> detList;
  void detListClear();
  Int_t fill()
  {
    return btree->Fill();
  }

  void addDet(int ichan)
  {
    TBRawEvent *det = new TBRawEvent(ichan);
    detList.push_back(det);
    btree->Branch(det->GetName(), det);
  }

  TBRawEvent *getDet(int ichan)
  {
    TBRawEvent *rev = NULL;
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

  void print()
  {
    printf(" %s  entries %lld \n", this->GetName(), btree->GetEntries());
    for (unsigned i = 0; i < detList.size(); ++i)
      printf(" %s ev %i %lu \n", detList[i]->GetName(), detList[i]->event, detList[i]->rdigi.size());
    btree->GetListOfBranches()->ls();
  }

  ClassDef(TBRawRun, 1)
};
#endif
