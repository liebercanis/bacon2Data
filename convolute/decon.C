/* uses TBaconRun Class */
//////////////////////////////////////////////////////////
// original: spulse.C M.Gold May 2020
// modified decon.C for SIPM deconvolution Dec. 24 2024
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

typedef std::complex<double> Complex;
using namespace TMath;

// spe landau shape
static double myLandau(Double_t *xx, Double_t *par)
{
  double x = xx[0];
  return par[2] * TMath::Landau(x, par[0], par[1], 1);
}

// const double sratio = 2.089764;
const double sratio = 10;
const double xtau = 0.00483;
const int nsamples = 7500;
const int ntrigger = 730;
const double sign = 1.0;
const double speMPV = double(ntrigger);
const double speSigma = 14.;
const double snoise = 1.E-3; // in SPE units

class decon
{
public:
  decon();
  virtual ~decon() { ; }
  TVirtualFFT *fFFT;
  TVirtualFFT *fInverseFFT;
  TRandom3 *ran;
  TFile *fout;
  TF1 *speLandau;
  TH1D *hResponseWave;
  TH1D *hInputWave;
  TH1D *hOutputWave;
  TH1D *hInputFFT;
  TH1D *hOutputFFT;
  double maxtime;
  void getResponse(double timeOffset, double numSPE);
  void getPulse(double timeOffset = 0, double numSPE = 1);
  std::vector<std::complex<double>> FFT(std::vector<double> vin);
  std::vector<std::complex<double>> inverseFFT(std::vector<std::complex<double>> complexVector);
  std::vector<double> getRealFFTWave(TH1D *h, std::vector<std::complex<double>> complexWave);
  TGraph *makeFFTGraph(std::vector<std::complex<double>> complexVector);
  std::vector<double> xval;
  std::vector<double> yval;
};

// make fft graph
TGraph *decon::makeFFTGraph(std::vector<std::complex<double>> complexVector)
{
  // get real, imaginary parts
  std::vector<double> freal;
  std::vector<double> fimag;
  std::vector<double> freq;
  for (unsigned i = 0; i < complexVector.size(); ++i)
  {
    freal.push_back(complexVector[i].real());
    fimag.push_back(complexVector[i].imag());
    freq.push_back(double(i));
  }
  TGraph *gfft = new TGraph(freq.size(), &freal[0], &fimag[0]);
  return gfft;
}

// make SPE pulse in time relative to trigger timeOffset and number of SPE numSPE
void decon::getResponse(double timeOffset, double numSPE)
{
  // set pulse parameters
  speLandau->SetParName(1, "MPV");
  speLandau->SetParameter(0, speMPV + timeOffset);
  speLandau->SetParName(1, "sigma");
  speLandau->SetParameter(1, speSigma);
  speLandau->SetParName(2, "norm");
  speLandau->SetParameter(2, numSPE);
  printf("getPulse with  MPV %f sigma %f norm %f nbins %i \n", speLandau->GetParameter(0), speLandau->GetParameter(1), speLandau->GetParameter(2), hResponseWave->GetNbinsX());

  // fill histogramW
  for (int i = 0; i < hResponseWave->GetNbinsX(); ++i)
  {
    // add in ramdom noise
    double val = speLandau->Eval(hResponseWave->GetBinCenter(i));
    // printf("getPulse bin %i val %f \n", i, val);
    hResponseWave->SetBinContent(i, val);
    if (i > int(speMPV) && val < numSPE * 1.E-9)
      break;
  }
}

// make SPE pulse in time relative to trigger timeOffset and number of SPE numSPE
void decon::getPulse(double timeOffset, double numSPE)
{
  // set pulse parameters
  speLandau->SetParName(1, "MPV");
  speLandau->SetParameter(0, speMPV + timeOffset);
  speLandau->SetParName(1, "sigma");
  speLandau->SetParameter(1, speSigma);
  speLandau->SetParName(2, "norm");
  speLandau->SetParameter(2, numSPE);
  printf("getPulse with  MPV %f sigma %f norm %f nbins %i \n", speLandau->GetParameter(0), speLandau->GetParameter(1), speLandau->GetParameter(2), hInputWave->GetNbinsX());

  // fill histogram
  for (int i = 0; i < hInputWave->GetNbinsX(); ++i)
  {
    // add in ramdom noise
    double val = speLandau->Eval(hInputWave->GetBinCenter(i)) + snoise * ran->Rndm();
    // printf("getPulse bin %i val %f \n", i, val);
    hInputWave->SetBinContent(i, val);
    if (i > int(speMPV) && val < numSPE * 1.E-9)
      break;
  }
}

// input is real wave to transform output is complex FFT
std::vector<std::complex<double>> decon::FFT(std::vector<double> vin)
{
  std::vector<std::complex<double>> complexVector;
  int nsamples = (int)vin.size();
  for (int is = 0; is < nsamples; ++is)
    fFFT->SetPoint(is, vin[is]);
  fFFT->Transform(); //

  std::vector<Double_t> realVec, imVec;
  for (int i = 0; i < nsamples; ++i)
  {
    double rl, im;
    fFFT->GetPointComplex(i, rl, im);
    std::complex<double> c(rl, im);
    c = c / sqrt(double(nsamples)); // normalize
    complexVector.push_back(c);
  }
  return complexVector;
}

// input is complex FFT output is complex FFT
std::vector<std::complex<double>> decon::inverseFFT(std::vector<std::complex<double>> complexInput)
{
  std::vector<Double_t> Signal;
  std::vector<std::complex<double>> complexOutput;
  int nsamples = (int)complexInput.size();
  for (int is = 0; is < nsamples; ++is)
  {
    fInverseFFT->SetPoint(is, complexInput[is].real(), complexInput[is].imag());
  }
  fInverseFFT->Transform();
  /*
  ** FFTW computes an unnormalized transform, in that there is no coefficient in front of the summation in the DFT.
  ** In other words, applying the forward and then the backward transform will multiply the input by n.
  ** */
  for (int i = 0; i < nsamples; ++i)
  {
    double rl, im;
    fInverseFFT->GetPointComplex(i, rl, im);
    std::complex<double> c(rl, im);
    c = c / sqrt(double(nsamples)); //.real or .imag accessors
    complexOutput.push_back(c);
  }
  return complexOutput;
}

// fill real FFT output
std::vector<double> decon::getRealFFTWave(TH1D *h, std::vector<std::complex<double>> complexWave)
{
  std::vector<double> realWave;
  for (int i = 0; i < complexWave.size(); ++i)
  {
    h->SetBinContent(i, complexWave[i].real());
    realWave.push_back(complexWave[i].real());
  }
  return realWave;
}

decon::decon()
{
  ran = new TRandom3();
  // initialize fft
  TString canName;
  int nfft = nsamples;
  fFFT = TVirtualFFT::FFT(1, &nfft, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nfft, "C2R M K");

  // histograms
  TFile *fout = new TFile("decon.root", "RECREATE");
  hResponseWave = new TH1D("ResponseWave", "ResponseWave", nsamples, 0, nsamples);
  hInputWave = new TH1D("InputWave", "InputWave", nsamples, 0, nsamples);
  hOutputWave = new TH1D("OutputWave", "OutoutWave", nsamples, 0, nsamples);
  hInputFFT = new TH1D("InputFFT", "InputFFT", nsamples, 0, nsamples);
  hOutputFFT = new TH1D("OutputFFT", "OutputFFT", nsamples, 0, nsamples);

  // landau pulse
  speLandau = new TF1("myLandau", myLandau, 0, double(nsamples), 3);

  // fill input into hInputWave
  getPulse();

  // setup for FFT
  xval.resize(hInputWave->GetNbinsX());
  yval.resize(hInputWave->GetNbinsX());
  for (int i = 0; i < hInputWave->GetNbinsX(); ++i)
  {
    xval[i] = hInputWave->GetBinCenter(i);
    yval[i] = hInputWave->GetBinContent(i);
  }

  // FFT of input
  std::vector<std::complex<double>> fftInputWave = FFT(yval);
  std::vector<double> realInputFFT = getRealFFTWave(hInputFFT, fftInputWave);

  // inverse FFT
  std::vector<std::complex<double>> fftOutputWave = inverseFFT(fftInputWave);
  std::vector<double> realOutputFFT = getRealFFTWave(hOutputFFT, fftOutputWave);

  // fill output wave
  for (int i = 0; i < hOutputWave->GetNbinsX(); ++i)
  {
    hOutputWave->SetBinContent(i, realOutputFFT[i]);
  }

  TGraph *gfft = makeFFTGraph(fftInputWave);
  gfft->SetName("SPE-InFFT");
  gfft->SetTitle("SPE-InFFT");
  gfft->SetMarkerSize(0.6);
  gfft->SetMarkerStyle(21);
  gfft->SetMarkerColor(kBlue);
  gfft->GetXaxis()->SetTitle("real");
  gfft->GetYaxis()->SetTitle("imag");
  fout->Append(gfft);

  TGraph *gfftOut = makeFFTGraph(fftInputWave);
  gfftOut->SetName("SPE-OutFFT");
  gfftOut->SetTitle("SPE-OutFFT");
  gfftOut->SetMarkerSize(0.8);
  gfftOut->SetMarkerStyle(4);
  gfftOut->SetMarkerColor(kRed);
  gfftOut->GetXaxis()->SetTitle("real");
  gfftOut->GetYaxis()->SetTitle("imag");
  fout->Append(gfftOut);

  // make canvas
  canName.Form("SPE-FFT");
  TCanvas *canFFT = new TCanvas(canName, canName);
  canFFT->SetGrid();
  gfft->Draw("ap");
  gfftOut->Draw("psame");
  canFFT->BuildLegend();

  fout->Write();
}
