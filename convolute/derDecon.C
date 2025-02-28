/* uses TBaconRun Class */
//////////////////////////////////////////////////////////
// original: spulse.C M.Gold May 2020
// modified decon.C for SIPM deconvolution Dec. 24 2024
// this version puts trigger at zero dec 27 2024
// add smoothing
// use derivative deconvolution Feb 26 2025
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
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TH1D.h>
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
#include <TSpectrum.h>
#include <TRandom3.h>

// SG filter
#include "../compiled/SGFilter.hh"
#include "../compiled/dConvolute.hh"

typedef std::complex<double> Complex;
using namespace TMath;

// spe landau shape
static double myLandau(Double_t *xx, Double_t *par)
{
  double x = xx[0];
  return par[2] * TMath::Landau(x, par[0], par[1], 1);
}

// const double sratio = 2.089764;
const double tiny = 1.E-44;
const double sratio = 10;
const int nsamples = 7500;
const int triggerOffset = 730 * 2; // ns
const int triggerTime = 0;
const double sign = 1.0;
const double speMPV = double(triggerOffset); // Landau cannot start from zero
const double speSigma = 2. * 14.;            // ns
const double spe = 1.;
int totalTime = nsamples * 2;

class derDecon
{
public:
  // derDecon(double noiseValue = 0.04, int ifilt = 7500);
  derDecon(double noiseValue = 0.04, int nwindowSG = 10, int npoly = 3);
  virtual ~derDecon() { ; }
  bool addNoise;
  int ismooth = 1;
  bool SGfilter = true; // apply SG smoothing
  int averWindow = 10;  // markov param 3,7,10
  double noiseSigma;
  // double filterMean = double(nsamples / 2);
  double filterMean = 0;
  double filterSigma = double(7500); // turn off
  TVirtualFFT *fFFT;
  TVirtualFFT *fInverseFFT;
  TRandom3 *ran;
  TFile *fout;
  TF1 *speLandau;
  // histos
  TH1D *hResponseWave;
  TH1D *hResponseWaveDer;
  TH1D *hHitWave;
  TH1D *hInputWave;
  TH1D *hInputWaveDer;
  TH1D *hOutputWave;

  double maxtime;
  std::vector<std::complex<double>> fftNoiseWave;
  std::vector<std::complex<double>> fftInputWave;
  void makeHistograms();
  std::vector<double> fillVector(TH1D *hist);
  void fillHisto(std::vector<double> v, TH1D *hist);
  std::vector<double> getResponse();
  std::vector<std::complex<double>> gResponse; // store respoinse
  void getHits(int nSinglet = 1, int nTriplet = 1);
  void hIntegrate(TH1D *hin, TH1D *hout);
};

std::vector<double> derDecon::fillVector(TH1D *hist)
{
  std::vector<double> vreturn;
  // hist starts with bin 1
  for (int ibin = 1; ibin < hist->GetNbinsX(); ++ibin)
    vreturn.push_back(hist->GetBinContent(ibin));
  return vreturn;
}

void derDecon::fillHisto(std::vector<double> v, TH1D *hist)
{
  // hist starts with bin 1
  for (int ibin = 0; ibin < v.size(); ++ibin)
    hist->SetBinContent(ibin + 1, v[ibin]);
}

// histogram based integration
void derDecon::hIntegrate(TH1D *hin, TH1D *hout)
{
  hout->SetBinContent(1, 0.);
  for (int ibin = 2; ibin < hin->GetNbinsX(); ++ibin)
    hout->SetBinContent(ibin, hin->GetBinContent(ibin) + hin->GetBinContent(ibin - 1));
}

void derDecon::makeHistograms()
{

  // output file and histograms with time in ns
  fout = new TFile("derDecon.root", "RECREATE");
  // set trigger time to zero
  hHitWave = new TH1D("HitWave", "HitWave", nsamples, -triggerOffset, totalTime - triggerOffset);
  hHitWave->GetXaxis()->SetTitle("time [ns]");

  hResponseWave = new TH1D("ResponseWave", "ResponseWave", nsamples, -triggerOffset, totalTime - triggerOffset);
  hResponseWave->GetXaxis()->SetTitle("time [ns]");

  hResponseWaveDer = new TH1D("ResponseWaveDer", "ResponseWaveDer", nsamples, -triggerOffset, totalTime - triggerOffset);
  hResponseWaveDer->GetXaxis()->SetTitle("time [ns]");

  hInputWave = new TH1D("InputWave", "InputWave", nsamples, -triggerOffset, totalTime - triggerOffset);
  hInputWave->GetXaxis()->SetTitle("time [ns]");

  hInputWaveDer = new TH1D("InputWaveDer", "InputWaveDer", nsamples, -triggerOffset, totalTime - triggerOffset);
  hInputWaveDer->GetXaxis()->SetTitle("time [ns]");
}

// make SPE pulse in time relative to trigger=0 and number of SPE numSPE
std::vector<double> derDecon::getResponse()
{
  std::vector<double> response;
  double numSPE = 1.;
  // landau response function
  speLandau = new TF1("myLandau", myLandau, 0, double(nsamples), 3);
  // set SPE response parameters
  speLandau->SetParName(1, "MPV");
  speLandau->SetParameter(0, speMPV);
  speLandau->SetParName(1, "sigma");
  speLandau->SetParameter(1, speSigma);
  speLandau->SetParName(2, "norm");
  speLandau->SetParameter(2, numSPE); // single SPE
  printf("line113 getResponse with  MPV %f sigma %f norm %f nbins %i \n", speLandau->GetParameter(0), speLandau->GetParameter(1), speLandau->GetParameter(2), nsamples);

  for (int i = 0; i < nsamples; ++i)
  {
    double val = TMath::Max(tiny, speLandau->Eval(double(i) + 0.5));
    response.push_back(val);
    hResponseWave->SetBinContent(i, val);
  }
  return response;
}

// make hInputHitWave
void derDecon::getHits(int nSinglet = 1, int nTriplet = 1)
{
  double sTau = 7.;
  double tTau = 1600.;
  for (int ir = 0; ir < nSinglet; ++ir)
  {
    double val = sTau * ran->Rndm() + double(triggerTime);
    int ibin = hHitWave->FindBin(val);
    hHitWave->SetBinContent(ibin, hHitWave->GetBinContent(ibin) + spe);
  }
  for (int ir = 0; ir < nTriplet; ++ir)
  {
    double val = tTau * ran->Rndm() + double(triggerTime);
    int ibin = hHitWave->FindBin(val);
    hHitWave->SetBinContent(ibin, hHitWave->GetBinContent(ibin) + spe);
  }
}
derDecon::derDecon(double noiseValue, int nwindowSG, int npoly) // 30/750
{
  ran = new TRandom3();
  // setup FFT in root
  int nfft = nsamples;
  fFFT = TVirtualFFT::FFT(1, &nfft, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nfft, "C2R M K");

  // init SG filter
  SavitzkyGolay *sgfilt = new SavitzkyGolay();
  // init deconvolution
  dConvolute *dConv = new dConvolute(fFFT, fInverseFFT);

  std::vector<std::complex<double>> gNoResponse; // zero length vector will do no convoluton
  std::vector<double> source;
  std::vector<double> rsource;
  std::vector<double> sgout; // vector for output of SG filter

  noiseSigma = noiseValue;
  if (SGfilter)
    noiseSigma *= 0.35;
  // filterSigma = double(ifilt); use default
  printf(" derDecon with noise sigma %E \n", noiseSigma);

  makeHistograms();

  // make response template
  std::vector<std::complex<double>> noResponseFFT;
  std::vector<double> response = getResponse();
  std::vector<std::complex<double>> responseFFT = dConv->FFT(response, noResponseFFT);
  // differentiate
  std::vector<double> dResponse = dConv->differentiate(response, 1);
  fillHisto(dResponse, hResponseWaveDer);
  // FFT without response
  std::vector<std::complex<double>> dResponseFFT = dConv->FFT(dResponse, noResponseFFT);

  // make noise vector
  std::vector<double> vnoise;
  for (int ibin = 0; ibin < nsamples; ++ibin) //
  {
    double gspe = ran->Gaus(0.0, noiseSigma);
    vnoise.push_back(gspe);
  }

  std::vector<double> dNoise = dConv->differentiate(vnoise);
  // noise FFT
  std::vector<std::complex<double>> dNoiseFFT = dConv->FFT(dNoise, noResponseFFT);

  // get hits and make hitWave
  int nsinglet = 7;
  int ntriplet = 7;
  printf("getHits %i %i \n", nsinglet, ntriplet);
  getHits(nsinglet, ntriplet);

  std::vector vhits = fillVector(hHitWave);

  // convolve with response
  std::vector<std::complex<double>> inputFFT = dConv->FFT(vhits, responseFFT);
  // deconvolute without Wiener
  std::vector<double> vinput = dConv->inverseFFT(inputFFT, dNoiseFFT, noResponseFFT, false, 0);

  // add nnoise
  for (unsigned i; i < vinput.size(); ++i)
    vinput[i] = vinput[i] + vnoise[i];

  fillHisto(vinput, hInputWave);
  // input derivative
  std::vector<double> dInput = dConv->differentiate(vinput);
  // quick fix
  cout << " dInput fix " << 1 << " " << dInput[1] << endl;
  dInput[1] = 0;

  fillHisto(dInput, hInputWaveDer);

  return;

  TString cname;
  cname.Form("output-%.2E-sigma-%.0f", noiseSigma, filterSigma);
  TCanvas *canOutput = new TCanvas(cname, cname);
  hOutputWave->SetLineColor(kGreen);
  hOutputWave->Draw();
  hHitWave->Draw("same");
  canOutput->Print(".pdf");

  fout->Write();
  return;
}
