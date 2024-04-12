/*  
  M.Gold June 2022
  class to make hits from vector data
  P. Zugec et al. Pulse processing routines for neutron time-of-flight data. Nucl. Instrum. Meth., A812:134–144, 2016.
*/
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

using namespace std;

typedef std::vector<std::pair<unsigned, unsigned>> peakType;
typedef std::vector<std::pair<unsigned, unsigned>>::iterator peakTypeIter;
typedef std::map<Double_t, TDetHit, std::less<Double_t>> hitMap;
typedef std::map<Double_t, TDetHit, std::less<Double_t>>::iterator hitMapIter;

//const Double_t qnorm = 1.0;

class hitFinder
{
public:
  enum // crossing types
  {
    PUP,
    UPDOWN,
    NUP,
    DOWNUP,
    PDOWN,
    NDOWN
  };
  enum
  {
    CAENLENGTH = 7500
  };
  TFile *fout;
  TBRun *tbrun;
  TString tag;
  TDirectory *fftDir;
  TDirectory *finderDir;
  TDirectory *splitDir;
  TDirectory *sumWaveDir;
  Long64_t theEvent;
  bool isCAEN;
  bool verbose;
  bool splitVerbose;
  bool smoothing;
  const Double_t qnorm = 1.0;
  hitFinder(TFile *theFile, TBRun *brun, TString theTag, int nSamples, vector<int> vchan, vector<double> sigmaValue);
  virtual ~hitFinder() { chanMap.clear(); }
  int nsamples;
  unsigned diffStep;
  unsigned thresholdStepSize;
  double hitThreshold;
  double derivativeThreshold;
  double peakThreshold;
  unsigned maxPeakLength;
  double QPEPeak;
  std::map<int, int> chanMap;
  vector<int> vChannel;
  TNtuple *ntFinder;
  TNtuple *ntSplit;
  std::vector<double> QPEnominal; // nominal QPE
  std::vector<double> rdigi;      //
  std::vector<double> digi;       // baseline subtracted
  std::vector<double> ddigi;      // derivative
  std::vector<double> sdigi;      // smoothed
  std::vector<double> hdigi;      // hits
  std::vector<double> pdigi;      // pulse derivative 
  std::vector<Double_t> fdigi;    // filtered
  std::vector<Double_t> SPEdigi;  // filtered
  std::vector<unsigned> crossings;
  std::vector<unsigned> crossingBin;
  std::vector<double> crossingTime;
  std::vector<unsigned> peakCrossings;
  std::vector<unsigned> peakCrossingBin;
  std::vector<double> peakCrossingTime;
  hitMap detHits;
  peakType peakList;
  std::vector<Int_t> peakKind;
  std::vector<unsigned> splitCount;
  vector<double> channelSigmaValue;
  vector<double> wfilter;
  double timeUnit;
  double microSec;
  void event(int idet, Long64_t ievent, vector<double> rdigi, double theDerivativeThreshold, double theHitThreshold, unsigned step = 3);
  void differentiate();
  vector<double> differentiate(int step, vector<double> pdigi);
  void findThresholdCrossings(Int_t idet, double thresh);
  void findDerivativeCrossings(Int_t idet);
  void findPeakCrossings(Int_t idet, unsigned peakStart, unsigned peakEnd);
  void makePeaks(int idet, std::vector<Double_t> v);
  hitMap makeHits(int idet, Double_t &triggerTime, Double_t &firstCharge);
  void plotWave(int idet, Long64_t jentry);
  void plotEvent(TDirectory* dir, unsigned ichan, Long64_t ievent);
  // careful with indicies ichan and idet!
  void printPeakList();
  void trimPeaks(int idet, std::vector<Double_t> v);
  void splitPeaks(int idet);
  bool gotTemplate;
  bool getTemplate(int ichan);
  //void templateFFT(std::vector<double> rdigi);
  TString templateFileName;
  void fillWFilter(int ichan);
  // hist
    TH1D *htemplate;
    TH1D *htemplateFFT;
    TH1D *hPeakCount;
    TH1D *hPeakValue;
    TH1I *hHitLength;
    TH1I *hPeakNWidth;
    TH1D *hPeakCrossingBin;
    TH1D *hPeakCrossingRatio;

    /// The fft class to take the fourier transform.
    TVirtualFFT *fFFT;
    /// The fft class to take the inverse fourier transform.
    TVirtualFFT *fInverseFFT;

    std::vector<std::complex<double>> forwardFFT(std::vector<double> rdigi);
    std::vector<Double_t> backwardFFT(std::vector<std::complex<double>> VectorComplex);
    std::vector<std::complex<double>> templateTransform;

    TH1D *hWFilter;
    TH1D *hDeriv8;
    std::vector<TH1D *> hCrossingBinA;
    std::vector<TH1D *> hCrossingBinB;
    std::vector<TH1D *> hCrossingBinC;
    std::vector<TH1D *> hCrossingMaxBin;
    std::vector<TH1D *> hFFT;
    std::vector<TH1D *> hInvFFT;
    std::vector<TH1D *> hFFTFilt;
    std::vector<TH1D *> hEvWave;
    std::vector<TH1D *> hPeakCut;
    std::vector<TH2D *> hPeakCutAndTime;
    std::vector<TH1D *> hEvHitPeakWave;
    std::vector<TH1D *> hEvSmooth;
    std::vector<TH1D *> hEvCross;
    std::vector<TH1D *> hEvPeakCross;
    std::vector<TH1D *> hEvHitWave;
    std::vector<TH1D *> hEvDerWave;
    std::vector<TH1D *> hEvFiltWave;
    std::vector<TH1D *> hDigiVal;
    std::vector<TH1D *> hDerivativeVal;
    std::vector<TH2D *> hDerivativeValTime;
    std::vector<TH1D *> hHitSum;
    std::vector<TH1D *> hUnFilteredSummedWave;
    std::vector<TH1D *> hFilteredSummedWave;
    TH1D* hEvAllSumWave;
  };
