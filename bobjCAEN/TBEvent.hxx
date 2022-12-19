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
#include "TDetHit.hxx"

using namespace std;

// class to store info for the run

class TBEvent : public TNamed
{
public:
  TBEvent(TString runName = "run0");
  //~TBEvent();

  Long64_t event;
  int channel;
  double time;
  Double_t energy;
  Double_t ave;
  Double_t sigma;
  Double_t derSigma;
  Int_t nspe;
  Double_t qPrompt;
  Double_t qSum;
  Double_t hitPrompt;
  Double_t hitSum;
  Int_t nhits;

  std::vector<TDetHit> hits;

  unsigned numhits() { return hits.size(); }

  void clear()
  {
    event = 0;
    channel = 0;
    time = 0;
    energy = 0;
    ave = 0;
    sigma = 0;
    derSigma = 0;
    nspe = 0;
    qPrompt = 0;
    qSum = 0;
    hitPrompt = 0;
    hitSum = 0;
    nhits = 0;
    hits.clear();
  }

  void print()
  {
    printf(" %s  event %lld  trigtime %f %lu \n", this->GetName(), event, time, hits.size());
  }

  ClassDef(TBEvent, 1)
};
#endif
