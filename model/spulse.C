/* uses TBaconRun Class */
//////////////////////////////////////////////////////////
//  M.Gold May 2020 
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <complex>//includes std::pair, std::make_pair
#include <valarray>
//
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TTree.h>
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
#include <algorithm>    // std::sort
#include <TSpectrum.h>
#include <TRandom3.h>

typedef std::complex<double> Complex;
using namespace TMath;

//const double sratio = 2.089764;
const double sratio = 10;
const double xtau = 0.00483;
const int nsamples = 10000;
const double sign = 1.0;

class spulse
{
  public:
    spulse();
    virtual ~spulse() { ; }
    TVirtualFFT *fFFT;
    TVirtualFFT *fInverseFFT;
    TFile *fout;
    double maxtime;
    std::vector<std::complex<double> > FFT(std::vector<double> vin);
    std::vector<Double_t > inverseFFT(std::vector<std::complex<double> > VectorComplex);
    void getPulse();
    std::vector<double>  xval;
    std::vector<double>  yval;


};

void spulse::getPulse()
{
  xval.clear(); yval.clear();
  int istart =0;
  int iwidth = 1;

  for(int i=0;  i<nsamples; ++i) {
    xval.push_back(double(i));
    double v=0;
    if(i>=istart&&i<istart+iwidth) v=sign;
    if(i>=istart+iwidth) v = -1.0*sign/sratio*Exp(-double(i-istart-iwidth)*xtau);
    yval.push_back(v);
  }

  printf(" getPulse %lu %lu \n",xval.size(),yval.size());
  return;

}



std::vector<std::complex<double> > spulse::FFT(std::vector<double> vin)
{
  std::vector<std::complex<double> > VectorComplex;
  int nsamples = (int) vin.size();
  for(int is =0; is<nsamples; ++is) fFFT->SetPoint(is, vin[is]);
  fFFT->Transform();//

  std::vector<Double_t> realVec,imVec;
  for (int i = 0; i<nsamples; ++i) {
    double rl, im;
    fFFT->GetPointComplex(i,rl,im);
    std::complex<double>  c(rl,im);
    c=c/sqrt(double(nsamples)); //.real or .imag accessors
    VectorComplex.push_back(c);
  }
  return VectorComplex;
}

std::vector< double > spulse::inverseFFT(std::vector<std::complex<double> > VectorComplex)
{
  std::vector<Double_t > Signal;
  int nsamples = (int)  VectorComplex.size();
  for(int is =0; is<nsamples; ++is) {
    fInverseFFT->SetPoint(is, VectorComplex[is].real(),VectorComplex[is].imag() );
  }
  fInverseFFT->Transform();
  /*
  ** FFTW computes an unnormalized transform, in that there is no coefficient in front of the summation in the DFT.
  ** In other words, applying the forward and then the backward transform will multiply the input by n.
  ** */
  for (int i = 0; i<nsamples; ++i) {
    double rl = fInverseFFT->GetPointReal(i)/sqrt(double(nsamples));
    Signal.push_back(rl);
  }
  return Signal;
}



spulse::spulse() 
{
   // initialize fft
  TString canName;
  int nfft = nsamples;
  fFFT = TVirtualFFT::FFT(1, &nfft, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nfft, "C2R M K");


  getPulse();

  TGraph *gr = new TGraph(xval.size(),&xval[0],&yval[0]);
  gr->SetName("SPEpulse");
  gr->SetTitle("SPEpulse");
  gr->SetLineWidth(2);
  gr->SetMarkerSize(0.4);
  gr->SetMarkerStyle(21);
  gr->SetMarkerColor(kBlue);

  
  // FFT of filtered wave
  std::vector<std::complex<double> > trans = FFT(yval);
  std::vector< double > vout = inverseFFT(trans);

  vector<double> freal;
  vector<double> fimag;
  vector<double> freq;

  for(unsigned i=0; i< vout.size()/2; ++i) {
    freal.push_back(trans[i].real());
    fimag.push_back(trans[i].imag());
    freq.push_back(double(i));
    //if(i%10000==0) printf(" %i vin %f freq %f rin %f rout %f \n",i,vin[i],freq[rin.size()-1],rin[rin.size()-1], rout[rin.size()-1] );
  }

  printf(" freq %lu \n",freq.size());

  TGraph* gfft  = new TGraph(freq.size(),&freal[0],&fimag[0]);
  gfft->SetTitle("SPE-FFT");
  gfft->SetName("SPEFFT ");
  gfft->SetMarkerSize(0.8);
  gfft->SetMarkerStyle(21);
  gfft->SetMarkerColor(kBlue);


  canName.Form("SPE-FFT");
  TCanvas *canFFT = new TCanvas(canName,canName);
  canFFT->SetGrid();
  gfft->GetXaxis()->SetTitle("real");
  gfft->GetYaxis()->SetTitle("imag");
  gfft->Draw("ap");


  TGraph *gout  = new TGraph(xval.size(),&xval[0],&vout[0]);
  gout->SetTitle("transformed-back");
  gout->SetName("transformed back");

  gout->SetMarkerSize(.2);
  gout->SetMarkerStyle(4);
  gout->SetMarkerColor(kRed);


  canName.Form("SPE-pulse");
  TCanvas *can = new TCanvas(canName,canName);
  can->SetGrid();
  gr->Draw("ap");
  gr->GetXaxis()->SetTitle("time");
  gr->GetYaxis()->SetTitle("response");

  gout->Draw("psame");
  can->BuildLegend();


  // build square wave and convolve
  vector<double> swave;
  unsigned width=1000;
  unsigned start = 1000;
  for(unsigned iw = 0; iw< xval.size(); ++iw) {
    if(iw<start) swave.push_back(0);
    else if(iw>start&&iw<start+width) swave.push_back(-1.0*sign);
    else swave.push_back(0.);
  }

  double norm=0;
  for(unsigned iw = 0; iw< swave.size(); ++iw) norm+= swave[iw];
  norm=1.0;

  std::vector<std::complex<double> > transWave = FFT(swave);
  unsigned maxFrequency = transWave.size();

  std::vector<std::complex<double> > rFFTWave;
  for(unsigned iw =0 ; iw<maxFrequency; ++iw) {
    std::complex<double> d = transWave[iw]*trans[iw];
    rFFTWave.push_back(d) ;
  }

  std::vector< double > rWave = inverseFFT(rFFTWave);

  // make graphs
  TGraph* gWave  = new TGraph(xval.size(),&xval[0],&swave[0]);
  gWave->SetTitle("squareWave");
  gWave->SetName("squareWave");

   // make graphs
  TGraph* gResWave  = new TGraph(xval.size(),&xval[0],&rWave[0]);
  gResWave->SetTitle("responseWave");
  gResWave->SetName("responseWave");
  gWave->SetMarkerSize(.2);
  gWave->SetMarkerStyle(4);
  gWave->SetMarkerColor(kBlue);
  gResWave->SetMarkerSize(.2);
  gResWave->SetMarkerStyle(4);
  gResWave->SetMarkerColor(kRed);


  canName.Form("SquareWave");
  TCanvas *canw = new TCanvas(canName,canName);
  canw->SetGrid();
  gWave->GetXaxis()->SetTitle("time");
  gWave->GetYaxis()->SetTitle("response");
  gResWave->GetXaxis()->SetTitle("time");
  gResWave->GetYaxis()->SetTitle("response");
  gResWave->Draw("ap");
  gWave->Draw("psame");
  canw->BuildLegend();




  TFile *fout = new TFile("spulse-output.root","recreate");
  fout->Append(gr);
  fout->Append(gfft);
  fout->Append(gWave);
  fout->Append(gResWave);
  fout->Write();


}
