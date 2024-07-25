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
  TDet(int ichan, int ilevel = -1);
  virtual ~TDet() {}

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
  Double_t mode;
  Double_t peak;
  Double_t totSum;  // waveform ADC sum
  Double_t preSum;  // pre trigger waveform ADC sum
  Double_t trigSum; // trigger window waveform ADC sum
  Double_t lateSum; // post trigger waveform ADC sum
  // sums of hit peaks
  Double_t totPeakSum;  // waveform ADC sum
  Double_t prePeakSum;  // pre trigger waveform ADC sum
  Double_t trigPeakSum; // trigger window waveform ADC sum
  Double_t latePeakSum; // post trigger waveform ADC sum

  Int_t thresholds;
  Int_t nspe;
  Double_t hitPrompt;
  Double_t qarea;
  Double_t qpeak;
  int pass;
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
    mode = 0;
    peak = 0;
    totSum = 0;
    preSum = 0;
    trigSum = 0;
    lateSum = 0;
    totPeakSum = 0;
    prePeakSum = 0;
    trigPeakSum = 0;
    latePeakSum = 0;
    thresholds = 0;
    nspe = 0;
    hitPrompt = 0;
    qarea = 0;
    qpeak = 0;
    pass = 0;
    hits.clear();
  }

  ClassDef(TDet, 11)
};
#endif
