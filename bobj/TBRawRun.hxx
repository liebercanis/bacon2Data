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
#include "TBEventData.hxx"

using namespace std;

// class to store info for the run

class TBRawRun : public TNamed
{
public:
  TBRawRun();
  TBRawRun(TString runName = "run0");
  virtual ~TBRawRun(){};

  void clear();
  TTree *btree;
  TBEventData *eventData;

  vector<TBRawEvent *> detList;
  void detListClear();
  Int_t fill()
  {
    return btree->Fill();
  }

  void updateTime(time_t evtime) {
    eventData->update(evtime);
  }

  TBRawEvent *addDet(unsigned ichan)
  {
    TBRawEvent *det = new TBRawEvent(ichan);
    detList.push_back(det);
    btree->Branch(det->GetName(), det);
    return det;
  }

  TBRawEvent *getDet(unsigned ichan)
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
    if(rev==NULL)
      rev = this->addDet(ichan);
    return rev;
  }

  void print()
  {
    printf(" %s  entries %lld \n", this->GetName(), btree->GetEntries());
    for (unsigned i = 0; i < detList.size(); ++i)
      printf(" %s time %llu  %lu \n", detList[i]->GetName(), detList[i]->time, detList[i]->rdigi.size());
    btree->GetListOfBranches()->ls();
  }

  ClassDef(TBRawRun, 1)
};
#endif
