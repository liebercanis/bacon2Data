/* uses TBaconRun Class */
//////////////////////////////////////////////////////////
// original: spulse.C M.Gold May 2020
// modified decon.C for SIPM deconvolution Dec. 24 2024
// this version puts trigger at zero dec 27 2024
// add smoothing
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
#include "SGFilter.hh"

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

class decon
{
public:
  // decon(double noiseValue = 0.04, int ifilt = 7500);
  decon(double noiseValue = 0.04, int nwindowSG = 10, int npoly = 3);
  virtual ~decon() { ; }
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
  TH1D *hResponseWaveUn;
  TH1D *hResponseWave;
  TH1D *hResponseFFT;
  TH1D *hResponseFFTUn;
  TH1D *hHitWave;
  TH1D *hInputWave;
  TH1D *hInputNoiseWave;
  TH1D *hInpuNoiseWaveSG;
  TH1D *hDataUnNoise;
  TH1D *hDataNoise;
  TH1D *hInputWaveShift;
  TH1D *hInputNoiseWaveShift;
  TH1D *hOutputWave;
  TH1D *hOutputNoiseWave;
  TH1D *hOutputNoiseWaveGold;
  TH1D *hOutputNoiseWaveSmooth;
  TH1D *hInputFFT;
  TH1D *hInputNoiseFFT;
  TH1D *hNoise;
  TH1D *hNoiseFFT;
  TH1D *hSignalWeight;
  TH1D *hResponseWeight;
  TH1D *hFilterNS;
  TH1D *hFilterWeight;
  double maxtime;
  std::vector<std::complex<double>> fftNoiseWave;
  std::vector<std::complex<double>> fftInputWave;
  void makeHistograms();
  void getResponse();
  std::vector<std::complex<double>> gResponse; // store respoinse
  void getHits(int nSinglet = 1, int nTriplet = 1);
  std::vector<std::complex<double>> FFT(TH1D *hin, std::vector<std::complex<double>> gResponse);
  std::vector<std::complex<double>> inverseFFT(std::vector<std::complex<double>> complexVector, std::vector<std::complex<double>> gResponse, bool filter = false, int shift = 0);
  void getRealFFTWave(TH1D *h, std::vector<std::complex<double>> complexWave, double norm = 1.);
  void getRealInvFFTWave(TH1D *h, std::vector<std::complex<double>> complexWave, double norm = 1.);
  TGraph *makeFFTGraph(std::vector<std::complex<double>> complexVector);
  void boxCar(TH1D *hist, TH1D *hsooth, int window = 3);
};

// implimintation of boxcar smoothing
void decon::boxCar(TH1D *hist, TH1D *hsmooth, int window)
{
  std::vector<double> source;
  source.resize(hist->GetNbinsX());
  printf(" boxCar size %lu \n", source.size());
  for (int i = 0; i < hist->GetNbinsX(); i++)
    source[i] = hist->GetBinContent(i + 1);

  double ave;
  for (int i = 0; i < hist->GetNbinsX(); i++)
  {
    int mstart = std::max(0, i - window);
    int mend = std::min(hist->GetNbinsX(), i + window);
    int norm = mend - mstart + 1;
    ave = 0;
    for (int j = mstart; j < mend; j++)
      ave += source[j];
    ave /= double(norm);
    // if (i / 100 * 100 == i)
    //   printf(" smooth bin %i val %f ave %f  norm %i \n ", i, source[i], ave, norm);
    hsmooth->SetBinContent(i + 1, ave);
  }
}

void decon::makeHistograms()
{

  // output file and histograms with time in ns
  fout = new TFile("decon.root", "RECREATE");
  // set trigger time to zero
  hHitWave = new TH1D("HitWave", "HitWave", nsamples, -triggerOffset, totalTime - triggerOffset);
  hHitWave->GetXaxis()->SetTitle("time [ns]");
  hNoise = new TH1D("Noise", "Noise", nsamples, -triggerOffset, totalTime - triggerOffset);
  hNoise->GetXaxis()->SetTitle("time [ns]");

  hResponseWaveUn = new TH1D("ResponseWaveUn", "ResponseWave unfiltered", nsamples, -triggerOffset, totalTime - triggerOffset);
  hResponseWaveUn->GetXaxis()->SetTitle("time [ns]");

  hResponseWave = new TH1D("ResponseWave", "ResponseWave", nsamples, -triggerOffset, totalTime - triggerOffset);
  hResponseWave->GetXaxis()->SetTitle("time [ns]");

  hInputWave = new TH1D("InputWave", "InputWave", nsamples, -triggerOffset, totalTime - triggerOffset);
  hInputWave->GetXaxis()->SetTitle("time [ns]");
  hInputNoiseWave = new TH1D("InputNoiseWave", "InputNoiseWave", nsamples, -triggerOffset, totalTime - triggerOffset);
  hInputNoiseWave->GetXaxis()->SetTitle("time [ns]");

  hInputWaveShift = new TH1D("InputWaveShift", "InputWaveShift", nsamples, -triggerOffset, totalTime - triggerOffset);
  hInputWaveShift->GetXaxis()->SetTitle("time [ns]");
  hInputNoiseWaveShift = new TH1D("InputNoiseWaveShift", "InputNoiseWaveShift", nsamples, -triggerOffset, totalTime - triggerOffset);
  hInputNoiseWaveShift->GetXaxis()->SetTitle("time [ns]");

  hInpuNoiseWaveSG = new TH1D("InputNoiseWaveSG", "InputNoiseWave SG filtered", nsamples, -triggerOffset, totalTime - triggerOffset);
  hInpuNoiseWaveSG->GetXaxis()->SetTitle("time [ns]");

  hDataUnNoise = new TH1D("DataUnNoise", "unfiltered Data noise", 200, -.4, .4);
  hDataUnNoise->GetXaxis()->SetTitle("SPE");

  hDataNoise = new TH1D("DataNoise", "SG filtered Data noise", 200, -.4, .4);
  hDataNoise->GetXaxis()->SetTitle("SPE");

  hOutputWave = new TH1D("OutputWave", Form("OutputWave noise %.2E", noiseSigma), nsamples, -triggerOffset, totalTime - triggerOffset);
  hOutputWave->GetXaxis()->SetTitle("time [ns]");

  hOutputNoiseWave = new TH1D("OutputNoiseWave", Form("OutputNoiseWave noise %.2E averWindow %i", noiseSigma, averWindow), nsamples, -triggerOffset, totalTime - triggerOffset);
  hOutputNoiseWave->GetXaxis()->SetTitle("time [ns]");

  hOutputNoiseWaveGold = new TH1D("OutputNoiseWaveGold", Form("OutputNoiseWave gold deconvolution noise %.2E", noiseSigma), nsamples, -triggerOffset, totalTime - triggerOffset);
  hOutputNoiseWaveGold->GetXaxis()->SetTitle("time [ns]");

  hOutputNoiseWaveSmooth = new TH1D(Form("OutputNoiseWaveSmooth%i", averWindow), Form("OutputNoiseWave smooth %i", averWindow), nsamples, -triggerOffset, totalTime - triggerOffset);
  hOutputNoiseWaveSmooth->GetXaxis()->SetTitle("time [ns]");

  // FFT histos
  int maxFreqency = 500.; // 1/2ns
  hResponseFFTUn = new TH1D("ResponseFFTUn", "ResponseFFT not SG filtered", nsamples / 2, 0, maxFreqency);
  hResponseFFTUn->GetXaxis()->SetTitle("frequency MHz");
  hInputFFT = new TH1D("InputFFT", "InputFFT", nsamples / 2, 0, maxFreqency);
  hInputFFT->GetXaxis()->SetTitle("frequency MHz");

  hResponseFFT = new TH1D("ResponseFFT", "ResponseFFT", nsamples / 2, 0, maxFreqency);
  hResponseFFT->GetXaxis()->SetTitle("frequency MHz");
  hInputFFT = new TH1D("InputFFT", "InputFFT", nsamples / 2, 0, maxFreqency);
  hInputFFT->GetXaxis()->SetTitle("frequency MHz");

  hInputNoiseFFT = new TH1D("InputNoiseFFT", "InputNoiseFFT", nsamples / 2, 0, maxFreqency);
  hInputNoiseFFT->GetXaxis()->SetTitle("frequency MHz");

  hNoiseFFT = new TH1D("NoiseFFT", "NoiseFFT", nsamples / 2, 0, maxFreqency);
  hNoiseFFT->GetXaxis()->SetTitle("frequency MHz");

  hSignalWeight = new TH1D("SignalWeight", "Signal weight", nsamples / 2, 0, maxFreqency);
  hSignalWeight->GetXaxis()->SetTitle("frequency MHz");

  hResponseWeight = new TH1D("ResponseWeight", "Response weight", nsamples / 2, 0, maxFreqency);
  hResponseWeight->GetXaxis()->SetTitle("frequency MHz");

  hFilterWeight = new TH1D("FilterWeight", Form("Filter weight noise=%.2E sigma %.0f", noiseSigma, filterSigma), nsamples / 2, 0, maxFreqency);
  hFilterWeight->GetXaxis()->SetTitle("frequency MHz");

  hFilterNS = new TH1D("FilterNS", Form("Filter noise/signal ratio noise=%.2E sigma %.0f", noiseSigma, filterSigma), nsamples / 2, 0, maxFreqency);
  hFilterNS->GetXaxis()->SetTitle("frequency MHz");
}

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

// input is real wave to transform output is complex FFT
std::vector<std::complex<double>> decon::FFT(TH1D *hin, std::vector<std::complex<double>> gResponse)
{
  std::vector<std::complex<double>> complexVector;
  int nsamples = hin->GetNbinsX();
  for (int is = 0; is < nsamples; ++is)
    fFFT->SetPoint(is, hin->GetBinContent(is));
  fFFT->Transform(); //

  std::vector<Double_t> realVec, imVec;
  for (int i = 0; i < nsamples; ++i)
  {
    double rl, im;
    fFFT->GetPointComplex(i, rl, im);
    std::complex<double> c(rl, im);
    // c = c / sqrt(double(nsamples)); // normalize
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
std::vector<std::complex<double>> decon::inverseFFT(std::vector<std::complex<double>> complexInput, std::vector<std::complex<double>> gResponse, bool filter, int shift = 0)
{
  std::vector<std::complex<double>> complexOutput;
  int nsamples = (int)complexInput.size();
  // deconvolution is C[w]=H[w]/G[w] with Wiener from wikipedia https://en.wikipedia.org/wiki/Wiener_deconvolution
  printf("deconvolute complex size %lu response %lu  filter %i \n", complexInput.size(), gResponse.size(), filter);
  if (gResponse.size() > 0)
  {

    double responseNorm = 0;
    double weigntNorm = 0;
    // loop over frequency bins
    for (int i = 0; i < complexInput.size(); ++i)
    {
      std::complex<double> W;
      if (filter)
      {
        std::complex<double> signalPower = std::norm(complexInput[i]);
        std::complex<double> noisePower = std::norm(fftNoiseWave[i]);
        std::complex<double> zero = 0.0;
        std::complex<double> inputPower = std::norm(fftInputWave[i]);
        double noiseSignalRatio = double(std::norm(fftNoiseWave[i])) / double(std::norm(signalPower));
        // modified with cut-off
        // Gaus args (x, filterMean, sigma , false = unnormalized = default)
        std::complex<double> noiseWeight = std::max(tiny, noiseSignalRatio * TMath::Gaus(double(i), filterMean, filterSigma));

        hFilterNS->SetBinContent(i + 1, noiseSignalRatio);
        hSignalWeight->SetBinContent(i + 1, std::abs(complexInput[i]));
        hResponseWeight->SetBinContent(i + 1, std::abs(gResponse[i]));
        hFilterWeight->SetBinContent(i + 1, std::abs(noiseWeight));

        responseNorm += std::norm(gResponse[i]);
        weigntNorm += std::norm(noiseWeight);

        // noiseSignalRatio = std::min(1.E-5, noiseSignalRatio);
        /* Vivek's suggestion
          std::complex<double> signalMinusNoise = signalPower - noisePower;

        if (std::norm(signalMinusNoise) < 0)
          signalMinusNoise = 0;
          */

        // printf("f=%i N/S %E \n", i, noiseSignalRatio);
        //  cut off
        //  if (noiseSignalRatio > 1)
        //   noiseSignalRatio = 1;
        std::complex denom = std::norm(gResponse[i]) + noiseWeight;
        W = std::conj(gResponse[i]) / (std::norm(gResponse[i]) + noiseWeight);
        std::complex<double> z2 = 1. / gResponse[i];
        // norm
        // W /= std::abs(W);
        // bunch of cross checks on complex
        std::complex<double> z1 = gResponse[i];
        std::complex<double> z3 = std::conj(gResponse[i]) / std::norm(gResponse[i]);
        // print some values
        /*
        if (i > -1)
        {
          printf("filter int %i ratio %E mag weight %E abs inv g %E abs denom %E abs W %E W = (%f,%f) inv g = (%f,%f)=(%f,%f) \n",
                 i, noiseSignalRatio, std::abs(noiseWeight), std::abs(z2), std::abs(denom), std::abs(W), W.real(), W.imag(), std::abs(z2), std::arg(z2), std::abs(z3), std::arg(z3));
        }
        */
      }
      else
      {
        W = 1. / gResponse[i];
        /*
        if (i > -1)
          printf("no filter int %i abs W %E \n", i, std::abs(W));
          */
        // apply transform
      }
      complexInput[i] *= W;
    }
  }
  // if applying time shift
  if (shift != 0)
    for (int i = 0; i < complexInput.size(); ++i)
    {
      // add a phase shift 1i is sqrt(-1) frequency i/nsamples
      std::complex<double> phase = std::exp(1i * double(shift) * double(i) / double(nsamples) * 2. * TMath::Pi());
      // cout << shift << " " << i << " phase " << phase << std::endl;
      complexInput[i] *= phase;
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
    // normalize
    c /= double(nsamples);
    complexOutput.push_back(c);
  }
  return complexOutput;
}

// fill real FFT output
void decon::getRealFFTWave(TH1D *h, std::vector<std::complex<double>> complexWave, double norm)
{
  // The DFT results are stored in-order in the array out, with the zero-frequency (DC) component in out[0]
  for (int i = 0; i < complexWave.size(); ++i)
    h->SetBinContent(i + 1, std::abs(complexWave[i].real()) * norm);
}

// fill real invFFT output
void decon::getRealInvFFTWave(TH1D *h, std::vector<std::complex<double>> complexWave, double norm)
{
  for (int i = 0; i < complexWave.size(); ++i)
    h->SetBinContent(i + 1, complexWave[i].real() * norm); // time bin
}

// make SPE pulse in time relative to trigger=0 and number of SPE numSPE
void decon::getResponse()
{
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

  // fill unfiltered response histogram
  for (int i = 0; i < hResponseWaveUn->GetNbinsX(); ++i)
  {
    double val = TMath::Max(tiny, speLandau->Eval(hResponseWaveUn->GetBinCenter(i)));
    hResponseWaveUn->SetBinContent(i, val);
  }
}

// make hInputHitWave
void decon::getHits(int nSinglet = 1, int nTriplet = 1)
{
  double sTau = 7.;
  double tTau = 1600.;
  for (int ir = 0; ir < nSinglet; ++ir)
  {
    double val = sTau * ran->Rndm() + double(triggerTime);
    int ibin = hHitWave->FindBin(val);
    printf("xxx bin %i time %f \n", ibin, val);
    hHitWave->SetBinContent(ibin, hHitWave->GetBinContent(ibin) + spe);
  }
  for (int ir = 0; ir < nTriplet; ++ir)
  {
    double val = tTau * ran->Rndm() + double(triggerTime);
    int ibin = hHitWave->FindBin(val);
    printf("xxx bin %i time  %f \n", ibin, val);
    hHitWave->SetBinContent(ibin, hHitWave->GetBinContent(ibin) + spe);
  }
}
decon::decon(double noiseValue, int nwindowSG, int npoly) // 30/750
{
  // init SG filter
  SavitzkyGolay *sgfilt = new SavitzkyGolay();
  TSpectrum *spect = new TSpectrum();
  std::vector<std::complex<double>> gNoResponse; // zero length vector will do no convoluton
  std::vector<double> source;
  std::vector<double> rsource;
  std::vector<double> sgout; // vector for output of SG filter

  addNoise = true;
  noiseSigma = noiseValue;
  if (SGfilter)
    noiseSigma *= 0.35;
  // filterSigma = double(ifilt); use default
  printf(" decon with noise sigma %E \n", noiseSigma);

  ran = new TRandom3();
  // initialize fft
  TString canName;
  int nfft = nsamples;
  fFFT = TVirtualFFT::FFT(1, &nfft, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nfft, "C2R M K");

  makeHistograms();

  // make response template and FFT of response and fill FFT gResponse
  getResponse();

  // FFT before SG
  // FFT of response function
  std::vector<std::complex<double>> gResponseUn = FFT(hResponseWaveUn, gNoResponse);
  getRealFFTWave(hResponseFFTUn, gResponseUn);

  // SG filter response
  sgout = sgfilt->SavGolFilter(hResponseWaveUn, nwindowSG, npoly);
  for (int i = 0; i < hResponseWave->GetNbinsX(); i++)
  {
    hResponseWave->SetBinContent(i + 1, sgout[i]);
  }

  // FFT of response function
  gResponse = FFT(hResponseWave, gNoResponse);
  getRealFFTWave(hResponseFFT, gResponse);
  TGraph *graphResponseFFT = makeFFTGraph(gResponse);
  graphResponseFFT->SetName("responseFFT");
  graphResponseFFT->SetTitle("responseFFT");
  graphResponseFFT->SetMarkerSize(0.6);
  graphResponseFFT->SetMarkerStyle(21);
  graphResponseFFT->SetMarkerColor(kBlue);
  graphResponseFFT->GetXaxis()->SetTitle("real");
  graphResponseFFT->GetYaxis()->SetTitle("imag");
  fout->Add(graphResponseFFT);

  // fill vedtor with response wave for Gold deconvolution
  rsource.clear();
  rsource.resize(hResponseWave->GetNbinsX());
  for (int i = 0; i < hResponseWave->GetNbinsX(); i++)
    rsource[i] = hResponseWave->GetBinContent(i + 1);

  // make hitWave
  int nsinglet = 7;
  int ntriplet = 7;
  printf("getHits %i %i \n", nsinglet, ntriplet);
  getHits(nsinglet, ntriplet);
  // FFT of input convolve with response
  printf("\t FFT convolute with response and fill hInputFFT, hInputWave\n");
  fftInputWave = FFT(hHitWave, gResponse);
  getRealFFTWave(hInputFFT, fftInputWave);

  // inverse FFT without deconvolution to get input wave
  int shift = 1500 - 40; // 3000 ns
  printf("\t inverse FFT without decovolution and fill hInputFFT, hInputWave with shift %i \n", shift);
  std::vector<std::complex<double>> fftOutputWave = inverseFFT(fftInputWave, gNoResponse, false);
  //  do not understand why I need to normalize hInputWave
  getRealInvFFTWave(hInputWave, fftOutputWave, sqrt(double(nsamples)));
  // for display only
  std::vector<std::complex<double>> fftOutputWaveShift = inverseFFT(fftInputWave, gNoResponse, false, shift);
  getRealInvFFTWave(hInputWaveShift, fftOutputWaveShift, sqrt(double(nsamples)));

  // add noise to hInputWave
  if (addNoise)
  {
    printf("\t add noise hist size %i noise sigma %f \n", hInputWave->GetNbinsX(), noiseSigma);
    for (int ibin = 0; ibin < hInputWave->GetNbinsX(); ++ibin) //
    {
      double gspe = ran->Gaus(0.0, noiseSigma);
      hNoise->SetBinContent(ibin, gspe);
      hInputNoiseWave->SetBinContent(ibin, hInputWave->GetBinContent(ibin) + gspe);
      hInputNoiseWaveShift->SetBinContent(ibin, hInputWaveShift->GetBinContent(ibin) + gspe);

      // if (ibin > triggerTime / 2 && ibin < triggerTime / 2 + 10)
      //  printf("addNoise hInputWave ibin %i was %f gspe %f content %f \n", ibin, hInputWave->GetBinContent(ibin), gspe, hInputWave->GetBinContent(ibin));
    }
  }

  // SG filter data
  sgout = sgfilt->SavGolFilter(hInputNoiseWaveShift, nwindowSG, npoly);
  for (int i = 0; i < hInpuNoiseWaveSG->GetNbinsX(); i++)
  {
    hInpuNoiseWaveSG->SetBinContent(i + 1, sgout[i]);
  }

  // noise histos
  for (int i = 0; i < hInputNoiseWave->GetNbinsX(); i++)
  {
    hDataUnNoise->Fill(hInputNoiseWave->GetBinContent(i));
    hDataNoise->Fill(hInpuNoiseWaveSG->GetBinContent(i));
  }

  // try gold deconvolution
  source.clear();
  source.resize(hInpuNoiseWaveSG->GetNbinsX());

  for (int i = 0; i < hInpuNoiseWaveSG->GetNbinsX(); i++)
    source[i] = hInpuNoiseWaveSG->GetBinContent(i + 1);

  spect->Deconvolution(&source[0], &rsource[0], int(source.size()), 1000, 1, 1);
  // refill histogram
  for (int i = 0; i < hOutputNoiseWaveGold->GetNbinsX(); i++)
    hOutputNoiseWaveGold->SetBinContent(i + 1, source[i]);
  // end of try gold

  // make FFT of the noise for Wiener filter
  fftNoiseWave = FFT(hNoise, gNoResponse);
  getRealFFTWave(hNoiseFFT, fftNoiseWave);

  // now FFT the input noise SG filtered wave without dconvolution
  std::vector<std::complex<double>> fftInputNoiseWave;
  if (SGfilter)
  {
    printf("\t FFT of hInpuNoiseWaveSG without response  / \n");
    fftInputNoiseWave = FFT(hInpuNoiseWaveSG, gNoResponse);
  }
  else
  {
    printf("\t FFT of hInputNoiseWave  without response  / \n");
    fftInputNoiseWave = FFT(hInputNoiseWave, gNoResponse);
  }

  // inverse FFT of input no filter
  printf("\t deconvolute fftInputWave using response and fill hOutputWave \n");
  std::vector<std::complex<double>> fftDeconvolveWave = inverseFFT(fftInputWave, gResponse);
  getRealInvFFTWave(hOutputWave, fftDeconvolveWave);

  // inverse FFT of signal applying filter
  printf("\t deconvolute fftInputNoiseWave using response X filter and fill hOutputWave \n");
  if (SGfilter)
    shift = 6040;
  else
    shift = 0;
  std::vector<std::complex<double>> fftDeconvolveNoiseWave = inverseFFT(fftInputNoiseWave, gResponse, true, shift);
  // again do not understand this normalization and it is noise dependent!
  getRealInvFFTWave(hOutputNoiseWave, fftDeconvolveNoiseWave);

  // smoothing
  if (ismooth == 1)
  {
    /* try smoothing just around singlet */
    int endBin = hOutputNoiseWave->FindBin(200.);
    int startBin = hOutputNoiseWave->FindBin(-200.);
    int ssize = endBin - startBin;
    // int ssize = pow(2, 10); // not sure if it needs to be power of 2
    // int nbins = hOutputNoiseWave->GetNbinsX();
    printf(" MARKOV smooth from bin %i center %.0f to bin %i center %.0f ssize %i\n", startBin, hOutputNoiseWave->GetBinCenter(startBin), endBin, hOutputNoiseWave->GetBinCenter(endBin), ssize);
    source.clear();
    source.resize(ssize);
    std::fill(source.begin(), source.end(), 0.0);
    // source.resize(hOutputNoiseWave->GetNbinsX());

    // for (int i = 0; i < hOutputNoiseWave->GetNbinsX(); i++)
    for (int i = startBin; i < endBin; i++)
      source[i - startBin] = hOutputNoiseWave->GetBinContent(i + 1);

    printf("smooth size %lu number %i \n", source.size(), averWindow);
    spect->SmoothMarkov(&source[0], ssize, averWindow);

    // refill histogram
    for (int i = 0; i < hOutputNoiseWave->GetNbinsX(); i++)
      hOutputNoiseWaveSmooth->SetBinContent(i + 1, hOutputNoiseWave->GetBinContent(i));

    for (int i = startBin; i < endBin; i++)
      hOutputNoiseWaveSmooth->SetBinContent(i + 1, source[i - startBin]);
  }
  else if (ismooth == 2)
  {
    boxCar(hOutputNoiseWave, hOutputNoiseWaveSmooth, averWindow);
  }
  else
  {
    sgout = sgfilt->SavGolFilter(hOutputNoiseWave, nwindowSG, npoly);
    for (int i = 0; i < hOutputNoiseWave->GetNbinsX(); i++)
    {
      hOutputNoiseWaveSmooth->SetBinContent(i + 1, sgout[i]);
    }
  }

  // why is hOututWave not normalized? and what is correct normalization?
  // for (int ibin = 0; ibin < hOutputWave->GetNbinsX(); ++ibin)
  //  hOutputWave->SetBinContent(ibin, hOutputWave->GetBinContent(ibin) / double(nsamples * nsamples));

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

  TString cname;
  cname.Form("input-%.2E-sigma-%.0f", noiseSigma, filterSigma);
  TCanvas *canInput = new TCanvas(cname, cname);
  hHitWave->SetLineColor(kRed);
  hInputNoiseWaveShift->SetLineColor(kBlue);
  hInputNoiseWaveShift->Draw();
  hHitWave->Draw("same");
  canInput->Print(".pdf");

  // hOutputWave->GetXaxis()->SetRangeUser(-100, 2000);
  // hOutputNoiseWave->GetXaxis()->SetRangeUser(-100, 2000);

  cname.Form("output-%.2E-sigma-%.0f", noiseSigma, filterSigma);
  TCanvas *canOutput = new TCanvas(cname, cname);
  hOutputWave->SetLineColor(kGreen);
  hOutputWave->Draw();
  hHitWave->Draw("same");
  canOutput->Print(".pdf");

  cname.Form("noise-%.2E-sigma-%.0f", noiseSigma, filterSigma);
  TCanvas *canWeight = new TCanvas(cname, cname);
  canWeight->SetLogy();
  hFilterWeight->SetLineColor(kRed);
  hFilterWeight->Draw("");
  hFilterNS->Draw("same");
  canWeight->Print(".pdf");

  if (ismooth == 1)
    cname.Form("InputNoise-%.2E-sigma-%.0f-MarkovSmooth%i", noiseSigma, filterSigma, averWindow);
  else if (ismooth == 2)
    cname.Form("InputNoise-%.2E-sigma-%.0f-BoxSmooth%i", noiseSigma, filterSigma, nwindowSG);
  else
    cname.Form("InputNoise-%.2E-sigma-%.0f-SGSmooth%i", noiseSigma, filterSigma, nwindowSG);
  TCanvas *canInputNoise = new TCanvas(cname, cname);
  hInpuNoiseWaveSG->SetLineColor(kMagenta);
  hInputNoiseWaveShift->SetLineColor(kGreen);
  hInputNoiseWaveShift->Draw("");
  hInpuNoiseWaveSG->Draw("same");
  // hHitWave->Draw("same");
  canInputNoise->Print(".pdf");

  if (ismooth == 1)
    cname.Form("OutputNoise-%.2E-sigma-%.0f-MarkovSmooth%i", noiseSigma, filterSigma, averWindow);
  else if (ismooth == 2)
    cname.Form("OutputNoise-%.2E-sigma-%.0f-BoxSmooth%i", noiseSigma, filterSigma, nwindowSG);
  else
    cname.Form("OutputNoise-%.2E-sigma-%.0f-SGSmooth%i", noiseSigma, filterSigma, nwindowSG);

  TCanvas *canOutputNoise = new TCanvas(cname, cname);
  hOutputNoiseWaveSmooth->SetLineColor(kMagenta);
  hOutputNoiseWave->SetLineColor(kGreen);
  hOutputNoiseWaveSmooth->Draw("");
  hOutputNoiseWave->Draw("same");
  hOutputNoiseWaveSmooth->Draw("same");
  hHitWave->Draw("same");
  canOutputNoise->Print(".pdf");

  /*

  TCanvas *canResponse = new TCanvas("response", "response");
  canResponse->SetLogy();
  hResponse->SetLineColor(kGreen);
  hResponse->Draw("");
  canResponse->Print(".pdf");
  */

  fout->Write();
}
