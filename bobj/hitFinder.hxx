////////////////////////////////////////////////////////
//  M.Gold June 2022
// class to make hits from vector data
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <complex> //includes std::pair, std::make_pair
#include <valarray>
//
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TDirectory.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm> // std::sort
#include "TSpectrum.h"
#include "TRandom3.h"
//
#include "TBWave.hxx"
#include "TDetHit.hxx"
#include "TBRun.hxx"

typedef std::vector<std::pair<unsigned, unsigned>> peakType;
typedef std::vector<std::pair<unsigned, unsigned>>::iterator peakTypeIter;
typedef std::map<Double_t, TDetHit, std::less<Double_t>> hitMap;
typedef std::map<Double_t, TDetHit, std::less<Double_t>>::iterator hitMapIter;

const Double_t qnorm = 1.0;

class hitFinder
{
public:
  enum
  {
    UPCROSS,
    DOWNCROSS,
    DOUBLEUPCROSS,
    DOUBLEDOWNCROSS
  };

  TFile *fout;
  TBRun *tbrun;
  TString tag;
  TDirectory *fftDir;
  hitFinder(TFile *theFile, TBRun *brun, TString theTag, int nSamples, vector<int> vchan );
  virtual ~hitFinder() { chanMap.clear(); }
  int nsamples;
  void plotWave(int idet, Long64_t jentry);
  void plotEvent( unsigned idet,Long64_t ievent);

  void findHits(int idet, Long64_t jentry);
  void differentiate(unsigned diffStep);
  std::map<int, int> chanMap;

  std::vector<double> rdigi;
  std::vector<double> digi;
  std::vector<double> ddigi;
  std::vector<double> hdigi;
  std::vector< Double_t > fdigi; // filtered
  hitMap detHits;
  peakType peakList;
  std::vector<Int_t> peakKind;
  double timeUnit;
  double microSec;
  void event(int idet, Long64_t ievent, vector<double> rdigi);
  void derivativePeaks(Int_t idet, Double_t rms);
  hitMap makeHits(int idet, Double_t &triggerTime, Double_t &firstCharge);
  void trimPeaks(int idet, std::vector<Double_t> v);
  void extendPeaks(int idet, std::vector<Double_t> v);
  bool getTransforms();
  // hist
  TH1D *hPeakCount;
  TH1I *hHitLength;
  TH1I *hPeakNWidth;

  /// The fft class to take the fourier transform.
  TVirtualFFT *fFFT;
  /// The fft class to take the inverse fourier transform.
  TVirtualFFT *fInverseFFT;

  std::vector<std::complex<double>> FFT(Int_t idet, std::vector<double> rdigi, bool first = true);
  std::vector<Double_t> inverseFFT(Int_t idet, std::vector<std::complex<double>> VectorComplex, std::vector<double> rdigi);
  bool gotTransforms;

  std::vector<TGraph *> gTransform;
  std::vector<TH1D *> hFFT;
  std::vector<TH1D *> hInvFFT;
  std::vector<TH1D *> hFFTFilt;
  std::vector<TH1D *> hEvWave;
  std::vector<TH1D *> hEvHitWave;
  std::vector<TH1D *> hEvDerWave;
  std::vector<TH1D *> hEvFiltWave;
  std::vector<TH1D *> hDigiVal;
};
