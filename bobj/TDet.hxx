/**
** MG, March 2020
**/
#ifndef TDET_DEFINED
#define TDET_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>
#include "TDetHit.hxx"

using namespace std;

// class to store info for the det

class TDet : public TNamed
{
public:
  TDet();
  TDet(int ichan, int ilevel=-1);
  // TDet::~TDet(){}

  // data elements
  Long64_t event;
  Long64_t trigger;
  Long64_t buffer;
  Int_t channel;
  Int_t level;
  Double_t ave;
  Double_t sigma;
  Double_t skew;
  Double_t base;
  Double_t peak;
  Double_t sum2;
  Double_t sum;
  Int_t crossings;
  Int_t thresholds;
  Int_t nspe;
  Double_t qPrompt;
  Double_t qSum;
  Double_t hitPrompt;
  Double_t hitSum;
  bool pass;
  std::vector<TDetHit> hits;

  unsigned nhits() { return hits.size(); }

  void clear()
  {
    event = 0;
    trigger = 0;
    buffer = 0;
    ave = 0;
    sigma = 0;
    skew = 0;
    base = 0;
    peak = 0;
    sum2 = 0;
    sum = 0;
    crossings = 0;
    thresholds = 0;
    nspe = 0;
    qPrompt = 0;
    qSum = 0;
    hitPrompt = 0;
    hitSum = 0;
    pass = 0;
    hits.clear();
  }

  ClassDef(TDet, 7)
};
#endif
