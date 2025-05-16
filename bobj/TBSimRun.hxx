/**
** MG, July 2020
**/
#ifndef TBSIMRUN_DEFINED
#define TBSIMRUN_DEFINED
#include <iostream>
#include <string>
#include <map>
#include <TNamed.h>
#include <TTree.h>
#include "TDet.hxx"

using namespace std;

// class to store info for the run

class TBSimRun : public TNamed
{
public:
  TBSimRun();
  TBSimRun(TString runName = "simRun0");
  virtual ~TBSimRun() {};

  /* members */
  TTree *btree;
  vector<TDet *> detList;

  /* inlines */
  void clear()
  {
    detListClear();
  }

  void detListClear()
  {
    for (unsigned i = 0; i < detList.size(); ++i)
      detList[i]->clear();
  }

  Int_t fill()
  {
    return btree->Fill();
  }

  TDet *addDet(unsigned ichan)
  {
    TDet *det = new TDet(ichan);
    detList.push_back(det);
    btree->Branch(det->GetName(), det);
    return det;
  }

  TDet *getDet(unsigned ichan)
  {
    TDet *rev = NULL;
    for (unsigned i = 0; i < detList.size(); ++i)
    {
      if (detList[i]->channel == ichan)
      {
        rev = detList[i];
        break;
      }
    }
    if (rev == NULL)
      rev = this->addDet(ichan);
    return rev;
  }

  void print()
  {
    printf(" %s  entries %lld \n", this->GetName(), btree->GetEntries());
    for (unsigned i = 0; i < detList.size(); ++i)
      printf(" channel %i name %s number of photons %lu \n", i, detList[i]->GetName(), detList[i]->hits.size());
    btree->GetListOfBranches()->ls();
  }

  ClassDef(TBSimRun, 1)
};
#endif
