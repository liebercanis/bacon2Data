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
const int nsamples = 7500;
const int ntrigger = 730 * 2; // ns
const double sign = 1.0;
const double speMPV = double(ntrigger);
const double speSigma = 2. * 14.; // ns
const double spe = 1.;
const double snoise = 1.E-6; // in SPE units

class decon
{
public:
  decon();
  virtual ~decon() { ; }
  bool addNoise;
  TVirtualFFT *fFFT;
  TVirtualFFT *fInverseFFT;
  TRandom3 *ran;
  TFile *fout;
  TF1 *speLandau;
  TH1D *hResponseWave;
  TH1D *hResponseFFT;
  TH1D *hHitWave;
  TH1D *hInputWave;
  TH1D *hOutputWave;
  TH1D *hInputFFT;
  double maxtime;
  void getResponse(double timeOffset = 0);
  std::vector<std::complex<double>> gResponse; // store respoinse
  void getHits(int nSinglet = 1, int nTriplet = 1);
  void getPulse(double timeOffset = 0, double numSPE = 1);
  std::vector<std::complex<double>> FFT(std::vector<double> vin, std::vector<std::complex<double>> gResponse);
  std::vector<std::complex<double>> inverseFFT(std::vector<std::complex<double>> complexVector, std::vector<std::complex<double>> gResponse);
  std::vector<double> getRealFFTWave(TH1D *h, std::vector<std::complex<double>> complexWave, int offset = 0, double norm = 1.);
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
    double val = speLandau->Eval(hInputWave->GetBinCenter(i)) + snoise * ran->Rndm();
    // printf("getPulse bin %i val %f \n", i, val);
    hInputWave->SetBinContent(i, val);
  }
}

// input is real wave to transform output is complex FFT
std::vector<std::complex<double>> decon::FFT(std::vector<double> vin, std::vector<std::complex<double>> gResponse)
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

  // convolution is H[w]=C[w]*G[w]
  printf("convolute complex size %lu \n", complexVector.size());
  if (gResponse.size() > 0)
    for (int i = 0; i < complexVector.size(); ++i)
      complexVector[i] *= gResponse[i];
  return complexVector;
}

// input is complex FFT output is complex FFT
std::vector<std::complex<double>> decon::inverseFFT(std::vector<std::complex<double>> complexInput, std::vector<std::complex<double>> gResponse)
{
  std::vector<std::complex<double>> complexOutput;
  int nsamples = (int)complexInput.size();
  // deconvolution is C[w]=H[w]/G[w]
  printf("deconvolute complex size %lu response %lu \n", complexInput.size(), gResponse.size());
  if (gResponse.size() > 0)
  {
    for (int i = 0; i < complexInput.size(); ++i)
    {
      complexInput[i] /= gResponse[i];
    }
  }
  for (int is = 0; is < nsamples; ++is)
  {
    fInverseFFT->SetPoint(is, complexInput[is].real(), complexInput[is].imag());
  }
  fInverseFFT->Transform();
  /*
  ** FFTW computes an unnormalized transform, in that there is no coefficient in front of the summation in the DFT.
  ** In other words, applying the forward and then the backward transform will multiply the input by n.
  ** */
  std::vector<Double_t> realVec, imVec;
  for (int i = 0; i < nsamples; ++i)
  {
    double rl, im;
    fInverseFFT->GetPointComplex(i, rl, im);
    std::complex<double> c(rl, im);
    c = c / sqrt(double(nsamples)); // normalize
    // c = c / double(nsamples);
    complexOutput.push_back(c);
  }
  return complexOutput;
}

// fill real FFT output
std::vector<double> decon::getRealFFTWave(TH1D *h, std::vector<std::complex<double>> complexWave, int offset, double norm)
{
  if (norm != 1.)
    printf("getRealFFTwave norm  %f \n", norm);
  std::vector<double> realWave;
  for (int i = 0; i < complexWave.size(); ++i)
  {
    // offset is in sample number = half of time in ns
    h->SetBinContent(i - offset / 2, complexWave[i].real() * norm);
    realWave.push_back(complexWave[i].real() * norm);
  }
  return realWave;
}

// make SPE pulse in time relative to trigger timeOffset and number of SPE numSPE
void decon::getResponse(double timeOffset)
{
  double numSPE = 1.;
  // landau response function
  speLandau = new TF1("myLandau", myLandau, 0, double(nsamples), 3);
  // set SPE response parameters
  speLandau->SetParName(1, "MPV");
  speLandau->SetParameter(0, speMPV + timeOffset);
  speLandau->SetParName(1, "sigma");
  speLandau->SetParameter(1, speSigma);
  speLandau->SetParName(2, "norm");
  speLandau->SetParameter(2, numSPE); // single SPE
  printf("line113 getResponse with  MPV %f sigma %f norm %f nbins %i \n", speLandau->GetParameter(0), speLandau->GetParameter(1), speLandau->GetParameter(2), nsamples);

  // fill histogram
  for (int i = 0; i < hResponseWave->GetNbinsX(); ++i)
  {
    double val = speLandau->Eval(hResponseWave->GetBinCenter(i));
    hResponseWave->SetBinContent(i, val);
    if (i > speMPV && val < numSPE * 1.E-9)
      break;
  }

  // containers for input to FFT
  yval.clear();
  yval.resize(hResponseWave->GetNbinsX());
  // make response FFT setup for FFT
  for (int i = 0; i < hResponseWave->GetNbinsX(); ++i)
  {
    yval[i] = hResponseWave->GetBinContent(i);
  }
  std::vector<std::complex<double>> gInput; // zero length
  // FFT of response function
  gResponse = FFT(yval, gInput);
  std::vector<double> realInputFFT = getRealFFTWave(hResponseFFT, gResponse);
  TGraph *graphResponseFFT = makeFFTGraph(gResponse);
  graphResponseFFT->SetName("responseFFT");
  graphResponseFFT->SetTitle("responseFFT");
  graphResponseFFT->SetMarkerSize(0.6);
  graphResponseFFT->SetMarkerStyle(21);
  graphResponseFFT->SetMarkerColor(kBlue);
  graphResponseFFT->GetXaxis()->SetTitle("real");
  graphResponseFFT->GetYaxis()->SetTitle("imag");
  fout->Add(graphResponseFFT);
}

// make hInputHitWave
void decon::getHits(int nSinglet = 1, int nTriplet = 1)
{
  double sTau = 7.;
  double tTau = 1600.;
  for (int ir = 0; ir < nSinglet; ++ir)
  {
    double val = sTau * ran->Rndm() + double(ntrigger);
    int ibin = hHitWave->FindBin(val);
    printf("xxx bin %i time %f \n", ibin, val);
    hHitWave->SetBinContent(ibin, hHitWave->GetBinContent(ibin) + spe);
  }
  for (int ir = 0; ir < nTriplet; ++ir)
  {
    double val = tTau * ran->Rndm() + double(ntrigger);
    int ibin = hHitWave->FindBin(val);
    printf("xxx bin %i time  %f \n", ibin, val);
    hHitWave->SetBinContent(ibin, hHitWave->GetBinContent(ibin) + spe);
  }
}
decon::decon()
{
  addNoise = true;
  ran = new TRandom3();
  // initialize fft
  TString canName;
  int nfft = nsamples;
  fFFT = TVirtualFFT::FFT(1, &nfft, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nfft, "C2R M K");

  // output file and histograms with time in ns
  fout = new TFile("decon.root", "RECREATE");
  hHitWave = new TH1D("HitWave", "HitWave", nsamples, 0, 2 * nsamples);
  hHitWave->GetXaxis()->SetTitle("time [ns]");
  hResponseWave = new TH1D("ResponseWave", "ResponseWave", nsamples, 0, 2 * nsamples);
  hResponseWave->GetXaxis()->SetTitle("time [ns]");
  hResponseFFT = new TH1D("ResponseFFT", "ResponseFFT", nsamples, 0, 2 * nsamples);
  hResponseFFT->GetXaxis()->SetTitle("time [ns]");
  hInputWave = new TH1D("InputWave", "InputWave", nsamples, 0, 2 * nsamples);
  hInputWave->GetXaxis()->SetTitle("time [ns]");
  hOutputWave = new TH1D("OutputWave", "OutputWave", nsamples, 0, 2 * nsamples);
  hOutputWave->GetXaxis()->SetTitle("time [ns]");
  hInputFFT = new TH1D("InputFFT", "InputFFT", nsamples, 0, 2 * nsamples);
  hInputFFT->GetXaxis()->SetTitle("time [ns]");

  // make response template and FFT of response and fill FFT gResponse
  getResponse();

  // make hitWave
  int nsinglet = 7;
  int ntriplet = 7;
  printf("getHits %i %i \n", nsinglet, ntriplet);
  getHits(nsinglet, ntriplet);

  // setup for FFT
  yval.clear();
  yval.resize(hHitWave->GetNbinsX());
  for (int i = 0; i < hHitWave->GetNbinsX(); ++i)
  {
    yval[i] = hHitWave->GetBinContent(i);
  }

  // FFT of input
  printf("FFT convolute with response and fill hInputFFT, hInputWave\n");
  std::vector<std::complex<double>> fftInputWave = FFT(yval, gResponse);
  std::vector<double> realInputFFT = getRealFFTWave(hInputFFT, fftInputWave);

  // inverse FFT
  int offset = ntrigger; // offset is needed because singnal does not start at time = 0
  printf("inverse FFT without decovolution and fill hInputFFT, hInputWave offset %i\n", offset);
  std::vector<std::complex<double>> gInput; // zero length to skip deconvolution
  std::vector<std::complex<double>> fftOutputWave = inverseFFT(fftInputWave, gInput);
  // offset is trigger time do not understand why I need to normalize hInputWave
  std::vector<double> realOutputFFT = getRealFFTWave(hInputWave, fftOutputWave, offset, double(fftOutputWave.size()));

  // copy of fft input without noise
  std::vector<std::complex<double>> fftInputNoiseWave = fftInputWave;
  // add noise on top of signal
  // From TRandom Gaus (Double_t mean=0, Double_t sigma=1)
  if (addNoise)
  {
    printf("\t add noise hist size %i noise sigma %f \n", hInputWave->GetNbinsX(), snoise);
    for (int ibin = 0; ibin < hInputWave->GetNbinsX(); ++ibin)
    {
      double gspe = ran->Gaus(0.0, spe * snoise);
      hInputWave->SetBinContent(ibin, hInputWave->GetBinContent(ibin) + gspe);
      // if (ibin > ntrigger / 2 && ibin < ntrigger / 2 + 10)
      //   printf("addNoise ibin %i gspe %f content %f \n", ibin, gspe, hInputWave->GetBinContent(ibin));
    }

    // now make FFT of input wave with noise added
    // setup for FFT
    yval.clear();
    yval.resize(hInputWave->GetNbinsX());
    for (int i = 0; i < hInputWave->GetNbinsX(); ++i)
    {
      yval[i] = hInputWave->GetBinContent(i);
    }

    // FFT of input without convolution
    printf("FFT of signal with no convolution and fill hInputFFT, hInputWave\n");
    fftInputNoiseWave = FFT(yval, gInput);
  }

  // inverse FFT
  printf("deconvolute inverse FFT using response and fill hOutputWave, hInputWave offset %i\n", offset);
  std::vector<std::complex<double>> fftDeconvolveWave = inverseFFT(fftInputNoiseWave, gResponse);
  offset = 0;
  std::vector<double> realDeconvolveFFT = getRealFFTWave(hOutputWave, fftDeconvolveWave, offset);

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

  TCanvas *canInput = new TCanvas("Input", "Input");
  hHitWave->SetLineColor(kRed);
  hInputWave->SetLineColor(kBlue);
  hInputWave->Draw();
  hHitWave->Draw("same");
  canInput->Print(".pdf");

  TCanvas *canOutput = new TCanvas("Output", "Output");
  hOutputWave->SetLineColor(kGreen);
  hOutputWave->Draw();
  hHitWave->Draw("same");
  canOutput->Print(".pdf");

  fout->Write();
}
