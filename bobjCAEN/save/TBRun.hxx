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
  
  //  data elements
  Long64_t nevents;
  TTree *btree;

//void dumpEvent(ofstream &dumpFile);

TString getDetName(int ichan){
  // check that it is not already in list
  TString detName;
  detName.Form("chan%i", ichan);
  return detName;
}

 bool addDet(int ichan, int ilevel=-1)
  {
    // check that it is not already in list
    TString detName = getDetName(ichan);
    if(btree->FindBranch(detName))
      return false;
    TDet *det = new TDet(ichan, ilevel);
    btree->Branch(det->GetName(), det);
    return true;
  }

  TDet *getDet(int ichan)
  {
    TString detName = getDetName(ichan);
    return (TDet*) btree->FindBranch(detName);
  }

  Int_t fill()
  {
    return btree->Fill();
  }

  void print()
  {
    //bevent->print();
    printf(" TBRun has %llu events %i dets \n", nevents,btree->GetListOfBranches()->GetEntries());
    btree->GetListOfBranches()->ls();
  }

  ClassDef(TBRun, 3)
};
#endif
