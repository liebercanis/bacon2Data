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

class TBRun: public TNamed {
  public:
    TBRun(TString runName = "run0");
    ~TBRun();

    void clear();
    TTree *btree;
    TBEvent *bevent;


    vector<TDet*> detList;
    void detListClear(); 

    void dumpEvent(ofstream& dumpFile);
 
    void fill() {
      btree->Fill();
    }

    void addDet(TString name=TString("UNKNOWN"), const TString title=TString("")) {
      TDet* det = new TDet(name,title);
      detList.push_back(det);
      btree->Branch(name, det);
    }

    void print() {
      bevent->print();
      for(unsigned i=0; i<detList.size(); ++i) printf(" %s %lu \n", detList[i]->GetName(), detList[i]->hits.size() );
      btree->GetListOfBranches()->ls();
    }

    ClassDef(TBRun,1)
};
#endif

