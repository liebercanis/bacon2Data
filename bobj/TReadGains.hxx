/**
** MG, August 18 2025
**/
#ifndef TREADGAINS_DEFINED
#define TREADGAINS_DEFINED
#include <iostream>
#include <string>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TNamed.h>
#include <vector>
#include <cmath>

using namespace std;

// class to store pmt hit

class TReadGains : public TNamed
{
public:
  enum
  {
    NUMCHANNELS = 13
  };

  TReadGains(bool useFile = true);
  virtual ~TReadGains()
  {
  }

  bool openFile();
  bool readPeakGains();
  bool readSumGains();
  void printGains();
  void clear();
  double getNominalPeak(unsigned ich);
  double getNominalSum(unsigned ich);

  // data elements
  bool readFromFile;
  TString gainFileName;
  TGraphErrors *gPeakGraph;
  TGraphErrors *gSumGraph;
  std::vector<double> sipmPeakGain;
  std::vector<double> sipmPeakGainError;
  std::vector<double> sipmSumGain;
  std::vector<double> sipmSumGainError;
  std::vector<double> relativeEff;
  double nominalGain;     // 170.;     // was 160.0; set Jue 13 2025
  double nominalTrigGain; //
  double nominalQsumGain;
  double nominalQsumTrigGain;
  double nominalPmtGain;
  double nominalQsumPmtGain;

  ClassDef(TReadGains, 1)
};
#endif
