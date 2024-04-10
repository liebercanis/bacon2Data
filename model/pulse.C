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


class pulse
{
  public:
    pulse();
    virtual ~pulse() { ; }
    TVirtualFFT *fFFT;
    TVirtualFFT *fInverseFFT;
    TFile *fout;
    int nsamples;
    double maxtime;
    std::vector<std::complex<double> > FFT(std::vector<double> vin);
    std::vector<Double_t > inverseFFT(std::vector<std::complex<double> > VectorComplex);
    std::vector<std::complex<double> > vpulse(std::vector<std::complex<double> >  vin);


};


std::vector<std::complex<double> > pulse::FFT(std::vector<double> vin)
{
  std::vector<std::complex<double> > VectorComplex;
  int nsamples = (int) vin.size();
  for(int is =0; is<nsamples; ++is) fFFT->SetPoint(is, vin[is]);
  fFFT->Transform();//

  std::vector<Double_t> realVec,imVec;
  for (int i = 0; i<nsamples; ++i) {
    double rl, im;
    fFFT->GetPointComplex(i,rl,im);
    std::complex<double>  c(rl,im); //.real or .imag accessors
    VectorComplex.push_back(c);
  }
  return VectorComplex;
}

std::vector< double > pulse::inverseFFT(std::vector<std::complex<double> > VectorComplex)
{
  std::vector<Double_t > Signal;
  int nsamples = (int)  VectorComplex.size();
  for(int is =0; is<nsamples; ++is) {
    fInverseFFT->SetPoint(is, VectorComplex[is].real(),VectorComplex[is].imag() );
  }
  fInverseFFT->Transform();
  for (int i = 0; i<nsamples; ++i) {
    double rl = fInverseFFT->GetPointReal(i);
    Signal.push_back(rl);
  }
  return Signal;
}


std::vector<std::complex<double>>  pulse::vpulse(std::vector<std::complex<double>>  vin)
{
  /*
  ** function vo = threecap(vi,w)
  ** This finds the voltage, as a function of frequency w and signal in, vi, 
  ** for a three-capacitor model of pmt output. 
  ** time is measured in nanoseconds.
  */
  double capu = 1.;//1.E-9;
  double R3 = 300;
  double R1 = 50;
  double C1 = 3*capu;
  double C2 = 30*capu;
  double C3 = 3*capu;
  double b = 0.25;
  std::complex<double>  ZR1(R1,0);
  std::complex<double>  ZR3(R3,0);
  std::vector<std::complex<double>> vout;

  for(unsigned i=1; i< vin.size(); ++i) {
    double w = 2.*TMath::Pi()*double(i)/double(vin.size());
    std::complex<double>  ZC1(0,1./(w*C1));
    std::complex<double>  ZC3(0,1./(w*C3));
    std::complex<double>  Z2(0,1./(w*C2));
    std::complex<double>  Z1 = ZC1 + ZR1;
    std::complex<double>  Z3 = ZR3*ZC3/(ZR3 + ZC3);
    std::complex<double>  i1 = vin[i]*(Z2 + b*Z3)/(Z1*Z2 + Z2*Z3 + Z1*Z3);
    std::complex<double>  vo = vin[i] - i1*ZC1;
    //printf(" vin %f %f current  %f %f vout %f %f  \n",norm(vin[i]),arg(vin[i]),norm(i1),arg(i1),norm(vo),arg(vo));
    vout.push_back(vo);
  }
  return vout;
}


pulse::pulse() 
{
  fout = new TFile("pulse.root","RECREATE");
  fout->cd();
  maxtime = 3000;
  nsamples = int(maxtime);
  // initialize fft
  fFFT = TVirtualFFT::FFT(1, &nsamples, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nsamples, "C2R M K");


  std::vector<double>  vin;
  vector<double> time;

  for(unsigned i=0; i<maxtime; ++i) {
    time.push_back(i);
    if(i<1) vin.push_back(1);
    else vin.push_back(0);
  }
  printf(" %lu %lu \n",time.size(),vin.size());

 
  // FFT of filtered wave
  std::vector<std::complex<double> > vinTrans;
  vinTrans = FFT(vin);


  std::vector<std::complex<double>> fvout=vpulse(vinTrans);
  std::vector< double > signal = inverseFFT(fvout);

  // norm
  double sum=0;
  for(unsigned i=0; i<signal.size(); ++i) {
    signal[i]/= double(nsamples);
    sum += signal[i];
  }

  printf(" %lu %lu %lu %lu \n",vin.size(),vinTrans.size(),fvout.size(),signal.size());
  printf(" sum %f to time %.0f \n",sum,maxtime);



  vector<double> rin;
  vector<double> rout;
  vector<double> freq;

  for(unsigned i=0; i< fvout.size()/2; ++i) {
    rin.push_back(vinTrans[i].real());
    rout.push_back(fvout[i].real());
    freq.push_back(double(i));
    //if(i%10000==0) printf(" %i vin %f freq %f rin %f rout %f \n",i,vin[i],freq[rin.size()-1],rin[rin.size()-1], rout[rin.size()-1] );
  }


  TGraph* gin  = new TGraph(freq.size(),&freq[0],&rin[0]);
  TGraph* gout = new TGraph(freq.size(),&freq[0],&rout[0]);
  gout->SetTitle(" FT vout ");
  gin->SetTitle(" FT vin ");


  TCanvas *can = new TCanvas("volt","volt");
  can->SetGridx(); can->SetGridy();
  can->SetTitle(" pulse ");
  gin->GetXaxis()->SetTitle(" freq");
  gin->GetYaxis()->SetTitle(" volt");
  gout->GetXaxis()->SetTitle(" freq");
  gout->GetYaxis()->SetTitle(" volt");
  gout->SetMarkerColor(kRed);
  gout->SetLineColor(kRed);
  gout->SetLineWidth(2);
  gout->SetMarkerStyle(7);
  gin->SetLineColor(kBlue);
  gin->SetLineWidth(2);
  gin->SetLineStyle(1);
  gout->SetLineStyle(1);
  gout->SetMarkerSize(1);
  gin->SetMarkerStyle(1);
  gin->SetMarkerColor(kBlue);
  gin->SetMarkerSize(.7);

  gout->Draw("ac");
  gin->Draw("csame");


  TGraph* gvin  = new TGraph(vin.size(),&time[0],&vin[0]);
  TGraph* gsignal = new TGraph(vin.size(),&time[0],&signal[0]);
  TCanvas *cani = new TCanvas("signal","signal");
  cani->SetGridx(); cani->SetGridy();
  gvin->SetTitle(" vin ");
  gvin->GetXaxis()->SetTitle(" time");
  gvin->GetYaxis()->SetTitle(" volt");
  gvin->SetLineColor(kBlack);
  gvin->SetLineWidth(3);
  gvin->SetLineStyle(1);
  gvin->SetMarkerStyle(4);
  gvin->SetMarkerSize(.4);
  gvin->SetMarkerColor(kBlack);

  gsignal->SetTitle(" vout ");
  gsignal->GetXaxis()->SetTitle(" time");
  gsignal->GetYaxis()->SetTitle(" volt");
  gsignal->SetLineColor(kRed);
  gsignal->SetLineWidth(3);
  gsignal->SetLineStyle(1);
  gsignal->SetMarkerStyle(22);
  gsignal->SetMarkerSize(.4);
  gsignal->SetMarkerColor(kRed);
  //gsignal->GetXaxis()->SetRangeUser(2,nsamples);
  //gsignal->GetYaxis()->SetRangeUser(-.002,.002);

  gsignal->Draw("apc");
  gvin->Draw("psame");



  fout->ls();
  fout->Write();
}
