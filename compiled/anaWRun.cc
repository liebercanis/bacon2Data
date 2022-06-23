////////////////////////////////////////////////////////
//  M.Gold August 2020
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
#include "TBRawEvent.hxx"
#include "TBRun.hxx"

typedef std::complex<double> Complex;
typedef std::map<Double_t,TDetHit,std::less<Double_t> >  hitMap;
typedef std::map<Double_t,TDetHit,std::less<Double_t> >::iterator  hitMapIter;

typedef std::vector<std::pair<unsigned,unsigned> >  peakType;
typedef std::vector<std::pair<unsigned,unsigned> >::iterator  peakTypeIter;

const Double_t qnorm=1.0;
class anaWRun
{
  public:
    enum {NDET=4};
    enum {UPCROSS,DOWNCROSS,DOUBLEUPCROSS,DOUBLEDOWNCROSS};
    TTree *btree;
    TBRawEvent *detList[NDET];
    ofstream dumpFile;
    TFile *fin;
    TFile *fout;
    TBRun *brun;
    TString runTag;
    enum { minLength=7};
    enum { baseLineHalfWindow = 200 }; // even integer

    anaWRun(Int_t maxEvents = 0, TString tag = "9_25_2020");
    virtual ~anaWRun() { ; }
    bool openFile();
    void analyze(Long64_t maxEvents);
    bool getWiggle();
    std::vector<Double_t> differentiate(std::vector<Double_t> v, unsigned diffStep); 
    void getAverage(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma,unsigned max=0);
    void getAverages(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma,unsigned max=0);
    void plotWave(Long64_t ievent,unsigned idet);
    std::vector<Double_t> getBaselineWeights(unsigned arraySize, peakType peakList, Int_t& maxWeightWidth);
    void getBaselineWMARecursive(int startBin, int stopBin, std::vector<Double_t> sig, std::vector<Double_t> weight, Int_t NWindow,std::vector<Double_t>& baseline);
    peakType  derivativePeaks(std::vector<Double_t> v, Int_t idet , Int_t nwindow, Double_t rms, std::vector<Int_t>& peakKind);
    peakType  simplePeaks(std::vector<Double_t> v, Int_t idet , unsigned minWidth, unsigned maxWidth, Double_t rms, std::vector<Int_t>& peakKind);
    void trimPeaks(int idet, std::vector<Double_t> v, peakType &peakList );
    void extendPeaks(int idet, std::vector<Double_t> v, peakType &peakList );
    bool getTransforms();
    /// The fft class to take the fourier transform.
    TVirtualFFT *fFFT;
    /// The fft class to take the inverse fourier transform.
    TVirtualFFT *fInverseFFT;

    std::vector<std::complex<double> > FFT(Int_t idet, std::vector<double> ddigi,bool first=true);
    std::vector<Double_t > inverseFFT(Int_t idet, std::vector<std::complex<double> > VectorComplex,std::vector<double> ddigi);



    hitMap makeHits(int idet, peakType peakList, std::vector<Int_t> peakKind, std::vector<Double_t> digi, Double_t sigma, Double_t& firstTime, Double_t& firstCharge);

    double timeUnit;
    double derivativeSigma[NDET];
    double vsign[NDET];
    double microSec;
    double maxLife;
    Int_t nSigma;
    Int_t aveWidth;
    unsigned long digiSize;
    unsigned diffStep;
    unsigned baselineBreak;
  

    TGraph *gTransform[NDET];
    TH1D* hSWavePMT;
    TH1D *hEvent[NDET];
    TH1D *hFFT[NDET];
    TH1D *hFFTFilt[NDET];
    TH1D *hInvFFT[NDET];
    TH1D *hNoise[NDET];
    TH1D *hBase[NDET];
    TH1D *hDerivative[NDET];
    TH1D *hEvRawWave[NDET];
    TH1D *hEvFiltWave[NDET];
    TH1D *hEvWave[NDET];
    TH1D *hEvDWave[NDET];
    TH1D *hEvSignal[NDET];
    TH1D *hEvPeaks[NDET];
    TH1D *hEvDPeaks[NDET];
    TH1D *hBaselineWMA[NDET];

    TH1D *hWave[NDET];
    TH1D *hSWaveRaw[NDET];
    TH1D *hSWaveSum[NDET];
    TH1D *hSWaveSumSubtracted[NDET];
    TH1D *hEvWeight[NDET];
    TH1D *hDPeaks[NDET];
    TH1D *hDNoise[NDET]; 
    TH1D *hSNoise[NDET];
    TH1D *hLife[NDET];
    TH1I* hPeakNWidth;
    TH1I* hHitLength;
    TH1D *hPeakCount;
    TNtuple *nt;
    TNtuple *ntDer;
    TNtuple *ntSimplePeaks;
    TNtuple *ntWave;
    TNtuple *ntCal;
    TNtuple *ntHit;
    TH1D* hPulseShape[NDET];
    TH1D* hWiggle;
    int pulseCount[NDET];

    int pmtMode;
    TDirectory *evDir;

};

bool anaWRun::getWiggle()
{
  TH1D* hW=NULL;
  hWiggle=NULL;
  TString fileName = TString("rootData/anaSum")+runTag+TString(".root");
  TFile *f1 = new TFile(fileName,"readonly");
  if(!f1) {
    printf(" no wiggle file for %s \n",runTag.Data());
    return false;
  } 
  f1->GetObject("SWaveSumwave1",hW);
  if(!hW) return false;
  hWiggle = (TH1D*) hW->Clone(Form("SWaveSumwave1-%s",runTag.Data())); 
  double norm = hWiggle->GetEntries()/hWiggle->GetNbinsX();
  printf(" got wiggle %s bins %i norm %.2f for run %s \n",hWiggle->GetName(),hWiggle->GetNbinsX(),norm,runTag.Data());
  for(int i=0; i< hWiggle->GetNbinsX(); ++i) hWiggle->SetBinContent(i, hW->GetBinContent(i)/norm);
  return true;
}

std::vector<std::complex<double> > anaWRun::FFT(Int_t id,std::vector<double> ddigi, bool first)
{
  std::vector<std::complex<double> > VectorComplex;
  int nsamples = (int) ddigi.size();
  for(int is =0; is<nsamples; ++is) fFFT->SetPoint(is, ddigi[is]);
  fFFT->Transform();//
  //cout << " FFT histogram " << hFFT[id]->GetName() << endl;

  std::vector<Double_t> realVec,imVec;
  for (int i = 0; i<nsamples; ++i) {
    double rl, im;
    fFFT->GetPointComplex(i,rl,im);
    std::complex<double>  c(rl,im); //.real or .imag accessors
    VectorComplex.push_back(c);
    // skip first bin which is pedestal
    if(i<nsamples/2) {
      if(first) hFFT[id]->SetBinContent(i,hFFT[id]->GetBinContent(i)+std::abs(c));
      else hFFTFilt[id]->SetBinContent(i,hFFTFilt[id]->GetBinContent(i)+std::abs(c));
    }
    realVec.push_back( VectorComplex[i].real());
    imVec.push_back(VectorComplex[i].imag() );
  }
  return VectorComplex;
}

std::vector< Double_t > anaWRun::inverseFFT(Int_t id,std::vector<std::complex<double> > VectorComplex,std::vector<double> ddigi)
{
  std::vector<Double_t > Signal;
  int nsamples = (int) detList[id]->digi.size();
  double sum=0;
  for(int is =0; is<nsamples; ++is) {
    fInverseFFT->SetPoint(is, VectorComplex[is].real(),VectorComplex[is].imag() );
    sum+= ddigi[is];
  }
  fInverseFFT->Transform();
  //cout << " FFT histogram " << hInvFFT[id]->GetName() << " ddigi sum " << sum << endl;

  Double_t norm = 0;
  for (int i = 0; i<nsamples; ++i) {
    double rl = fInverseFFT->GetPointReal(i);
    norm += rl;
    Signal.push_back(rl);
  }
  for(int i = 0; i < nsamples; i++){
    Signal[i] = Signal[i]*sum/norm;
    hInvFFT[id]->SetBinContent(i+1,hInvFFT[id]->GetBinContent(i+1) + Signal[i]);
  }
  return Signal;
}


bool anaWRun::getTransforms()
{
  TGraph* gw=NULL;
  TString fileName = TString("WienerTransforms.root");
  TFile *f1 = new TFile(fileName,"readonly");
  if(!f1) {
    printf(" no  file for %s \n",runTag.Data());
    return false;
  } 
  bool val=true;
  for(int i=0; i<NDET; ++i ) {
    gTransform[i]=NULL;
    f1->GetObject(Form("WienerTransDet%i",i),gw);
    if(!gw) {
      val=false;
      continue;
    }
    gTransform[i] = (TGraph*) gw->Clone(Form("WeinerForDet%i-%s",i,runTag.Data())); 
    printf(" got transform %s points %i for run %s \n",gTransform[i]->GetName(),gTransform[i]->GetN(),runTag.Data());
    fout->Append(gTransform[i]);
  }
  return val;
}


////////////////////////////////////////////////////////
anaWRun::anaWRun(Int_t maxEvents, TString tag)
{
  runTag = tag;
  hWiggle=NULL;
  pmtMode=1;
  fin=NULL;
  timeUnit = 2.0 ; // ns per count
  microSec=1.0E-3;
  diffStep=7;
  nSigma=5;
  aveWidth=20;
  baselineBreak=3000;

  if(!openFile()) { cout << " file with tag not found " << endl; return; }
  cout <<  " opened file with tag " << endl;

  btree->GetEntry(0);
  for(unsigned id=0; id< NDET; ++id) {
    pulseCount[id]=0;
    derivativeSigma[id]=3.0;
    vsign[id]=1.0;
    TString dname = detList[id]->GetName();
    if(dname.Contains( TString("1")))  {
      derivativeSigma[id]=5.0;
      vsign[id]=-1.0;
    }
    printf(" setting  derivative Sigma for %s to %.3f \n",dname.Data(),derivativeSigma[id]);
  }
  double baseOffset = 8100;
  int nsamples = (int) detList[0]->digi.size();
  maxLife = double(nsamples)*timeUnit*microSec; // microseconds
  cout << tag << " samples " << nsamples << endl;
  for(unsigned id=0; id<NDET; ++id) { 
    hEvRawWave[id] = new TH1D(Form("EvRawWave%s",detList[id]->GetName()),Form("RawWave%s",detList[id]->GetName()), nsamples,0, nsamples);
    hEvFiltWave[id] = new TH1D(Form("EvFiltWave%s",detList[id]->GetName()),Form("RawWave%s",detList[id]->GetName()), nsamples,0, nsamples);
    hEvWave[id] = new TH1D(Form("EvWave%s",detList[id]->GetName()),Form("Wave%s",detList[id]->GetName()), nsamples,0, nsamples);
    //hFFT[id] = new TH1D(Form("FFTDET%i",id),Form("FFT Channel %i ",id),nsamples/2,0,nsamples/2);
    //hInvFFT[id] = new TH1D(Form("InvFFTDET%i",id),Form("Inverse FFT Channel %i ",id),nsamples,0,nsamples);
    hEvDWave[id] = new TH1D(Form("EvDWave%s",detList[id]->GetName()),Form("DWave%s",detList[id]->GetName()), nsamples,0, nsamples);
    hEvSignal[id] = new TH1D(Form("EvSignal%s", detList[id]->GetName()), Form("DetSignal%s", detList[id]->GetName()), nsamples, 0, nsamples);
    hEvPeaks[id] = new TH1D(Form("EvPeaks%s", detList[id]->GetName()), Form("DetPeaks%s", detList[id]->GetName()), nsamples, 0, nsamples);
    hEvDPeaks[id] = new TH1D(Form("EvDPeaks%s", detList[id]->GetName()), Form("DerivDetPeaks%s", detList[id]->GetName()), nsamples, 0, nsamples);
    hEvWeight[id] = new TH1D(Form("EvWeight%s",detList[id]->GetName()),Form("baseline weight det %s",detList[id]->GetName()), nsamples,0, nsamples);
    hBaselineWMA[id]=new TH1D(Form("BaselineWMA%s", detList[id]->GetName()), Form("BaselineWMA%s", detList[id]->GetName()), nsamples, 0, nsamples);
  }

  if(pmtMode==0) fout = new TFile( Form("ana-Filt-%s.root",tag.Data()),"recreate");
  else  if(pmtMode==1) fout = new TFile( Form("ana-%s.root",tag.Data()),"recreate");
  else fout = new TFile( Form("ana-NoFilt-%s.root",tag.Data()),"recreate");
  
  printf(" \t **** pmtMode %i output file name is %s **** \n",pmtMode,fout->GetName());

  if(!getTransforms()) {
    printf(" !!!! no transforms !!!! \n");
    return;
  }

  if(pmtMode==1) {
    if(!getWiggle()) {
      printf(" !!!! no wigge histograms !!!! \n");
      return;
      if(hWiggle) fout->Append(hWiggle);
    }
  }

  evDir = fout->mkdir("EvDir");
  fout->cd();

  cout << tag << " samples " << nsamples << endl;
  for(unsigned id=0; id<NDET; ++id) { 
    hFFT[id] = new TH1D(Form("FFTDET%i",id),Form("FFT Channel %i ",id),nsamples/2,0,nsamples/2);
    hInvFFT[id] = new TH1D(Form("InvFFTDET%i",id),Form("Inverse FFT Channel %i ",id),nsamples,0,nsamples);
    hFFTFilt[id] = new TH1D(Form("FFTFiltDET%i",id),Form("filtered FFT Channel %i ",id),nsamples/2,0,nsamples/2);
  }


  brun = new TBRun(btree,tag);

  printf(" ***** assigned det names *****  %lu \n", brun->detList.size() );
  for(unsigned id=0; id<brun->detList.size() ; ++id) printf(" %i %s %s \n",id,brun->detList[id]->GetName() , detList[id]->GetName() );

  // open dump file
  dumpFile.open (Form("dump-%s.txt",brun->bevent->GetName()),ios::out);
  dumpFile << "C****** open  dump for run " <<  brun->bevent->GetName() << endl;
  cout << " ****** open  dump for run " <<  brun->bevent->GetName() << endl;


  hSWavePMT = new TH1D("SWavePMT","Wave PMT", nsamples,0, nsamples);
  for(unsigned id=0; id<NDET; ++id) { 
    hPulseShape[id] = new TH1D(Form("PulseShape%s",detList[id]->GetName()),Form("PulseShape%s",detList[id]->GetName()), nsamples,0, nsamples);
    hSWaveRaw[id] = new TH1D(Form("SWaveRaw%s",detList[id]->GetName()),Form("Wave%s",detList[id]->GetName()), nsamples,0, nsamples);
    hSWaveSum[id] = new TH1D(Form("SWaveSum%s",detList[id]->GetName()),Form("Wave%s",detList[id]->GetName()), nsamples,0, nsamples);
    hSWaveSumSubtracted[id] = new TH1D(Form("SWaveSumSubtracted%s",detList[id]->GetName()),Form("Wave%s",detList[id]->GetName()), nsamples,0, nsamples);
    hEvent[id] = new TH1D(Form("Event%s",detList[id]->GetName()),Form("Sum%s",detList[id]->GetName()), nsamples,0, nsamples);
    hBase[id] = new TH1D(Form("RawBase%s",detList[id]->GetName()),Form("Base det %s",detList[id]->GetName()), 400,baseOffset-200,baseOffset+200);
    hNoise[id] = new TH1D(Form("RawNoise%s",detList[id]->GetName()),Form("Noise det %s",detList[id]->GetName()), 50,0, 50);
    hNoise[id]->GetXaxis()->SetTitle(" ADC counts ");
    hNoise[id]->GetYaxis()->SetTitle(" frequency ");
    hDNoise[id] = new TH1D(Form("DNoise%s",detList[id]->GetName()),Form("DNoise%s",detList[id]->GetName()), 1000,-200,200);
    hSNoise[id] = new TH1D(Form("SNoise%s", detList[id]->GetName()), Form("SNoise det %s", detList[id]->GetName()), 50, -50, 50);
    hLife[id] = new TH1D(Form("Life%i", id), Form(" lifetime DET %i ", id), 4000, 0, maxLife);
    hLife[id]->Sumw2();
  }
  hPeakCount = new TH1D("PeakCount"," peaks by det ",NDET,0,NDET);
  hHitLength = new TH1I("HitLength", " hit length", 100, 0, 100);
  hPeakNWidth = new TH1I("PeakNWidth", "PeakNWidth", 100, 0, 100);
  nt = new TNtuple("ntBase", Form("ntBase%s", tag.Data()), "id:rawBase:rawNoise:rawAve:rawSigma");
  ntDer =  new TNtuple("ntDer"," deriviative ","t:sigma:d0:kover:type");
  ntSimplePeaks =  new TNtuple("ntSimplePeaks"," simple peaks ","idet:rms:v:t:width");
  ntWave = new TNtuple("ntWave"," wave ","event:v:d");
  ntCal =  new TNtuple("ntCal","ntuple Cal","iev:ipmt:base:sigma:dbase:dsigma:sbase:ssigma");
  ntHit = new TNtuple("ntHit", "ntuple Hit", "det:hits:order:istart:time:q:width:peak:good:kind");


  // initialize fft
  fFFT = TVirtualFFT::FFT(1, &nsamples, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nsamples, "C2R M K");
  nsamples = (int) detList[0]->digi.size();

  analyze(maxEvents);

  // normalize pulse shapes
  for(unsigned id=0; id<NDET; ++id) {
    printf(" pulse shapes averaged over %.0E for det %i \n",float(pulseCount[id]),id);
    if(pulseCount[id]>0) {
      hPulseShape[id]->Scale(1./double(pulseCount[id]));
    }
    hPulseShape[id]->SetTitle(Form(" pulse shape det %i n=%i",id,pulseCount[id]));
  }

  fin->Close();
  fout->Write();
  if (dumpFile.is_open()) { 
    cout << " ****** close dump for run " <<  brun->bevent->GetName() << endl;
    dumpFile << "C****** close dump for run " <<  brun->bevent->GetName() << endl;
    dumpFile.close();
  }


//fout->Close()

}

void anaWRun::analyze(Long64_t nmax)
{
  peakType peakLateList;
  peakType peakList;
  std::vector<Double_t> baselineDigi;
  std::vector<Int_t> peakKind;
  std::vector<double> vderiv;
  std::vector<double> ddigi;
  std::vector<double> sdigi;
  hitMap  detHits;
      
  unsigned prePulse=800;


  // reset wave histos
  for (unsigned id=0; id<NDET; ++id) {
    hEvWave[id]->Reset();
    hEvRawWave[id]->Reset();
    hEvFiltWave[id]->Reset();
    //hFFT[id]->Reset();
    hInvFFT[id]->Reset();
    hEvDWave[id]->Reset();
    hEvSignal[id]->Reset();
    hEvPeaks[id]->Reset();
    hEvDPeaks[id]->Reset();
    hEvWeight[id]->Reset();
    hBaselineWMA[id]->Reset();
  }
  //
  Long64_t nentries = btree->GetEntries();
  if(nmax>0) nentries=nmax;
  cout << " \t pmtMode " << pmtMode << " analyze " <<  nentries << endl;
  double trigTime=0;
  for(Long64_t ievent=0; ievent< nentries; ++ievent) {
    if (ievent%100==0) printf(" ... %lld \n", ievent);
    brun->clear();
    btree->GetEntry(ievent);
    // fill raw ntuple... these are not used
    double rawBaseline[NDET];
    double noise[NDET];
    for(unsigned id=0; id<NDET; ++id) {
      for(unsigned long i=0; i< detList[id]->digi.size() ; ++i) hEvent[id]->SetBinContent(i+1, hEvent[id]->GetBinContent(1+i)+detList[id]->digi[i]);
      vector<double> adigi;
      for (unsigned is=0; is<prePulse; ++is) adigi.push_back(detList[id]->digi[is]);
      std::sort(adigi.begin(), adigi.end());
      double nsamples = (double) detList[id]->digi.size();
      rawBaseline[id] = adigi[0.5*prePulse];
      noise[id] = std::abs( adigi[0.68*prePulse] - rawBaseline[id]);
    }

    fout->cd();
    bool isPMT=false;
    for(int id = 0; id < NDET; id++){
      double firstPulseTime=0;
      TString dname = detList[id]->GetName();
      if(dname.Contains( TString("1"))) isPMT=true;
      else isPMT=false;
      Double_t derAve, derSigma;
      Double_t rawAve, rawSigma;
      getAverages(detList[id]->digi, rawAve, rawSigma, prePulse); // first bins
      hBase[id]->Fill(rawAve);
      hNoise[id]->Fill(rawSigma);


      // first subtract rawBaseline
      nt->Fill(id,rawBaseline[id],noise[id],rawBaseline[id],noise[id],rawAve,rawSigma);
      //printf(" det %i raw ave %f base %f noise %f %f \n",id,rawAve,rawBaseline[id] ,rawSigma, noise[id] );
      ddigi.clear();
      for (unsigned i = 0; i < detList[id]->digi.size(); i++) ddigi.push_back(detList[id]->digi[i]-rawAve); 
      digiSize =  ddigi.size();

      for(unsigned isample = 0; isample < ddigi.size(); isample++) hEvRawWave[id]->SetBinContent(isample+1,ddigi[isample]);
      for (unsigned i = 0; i < detList[id]->digi.size(); i++) hSWaveRaw[id]->SetBinContent(i+1,hSWaveRaw[id]->GetBinContent(i+1)+ddigi[i]);
     

      // make the derivative waveform 
      vderiv.clear();
      vderiv = differentiate(ddigi,diffStep);
      getAverages(vderiv,derAve,derSigma);
      for(unsigned isample = 0; isample < vderiv.size(); isample++){
        hEvWave[id]->SetBinContent(isample+1,ddigi[isample]);
        hEvDWave[id]->SetBinContent(isample+1,vderiv[isample]);
        ntWave->Fill(id,ddigi[isample],vderiv[isample]);
      }

      //FFT and filter
      std::vector<std::complex<double> > complexTransform;
      complexTransform = FFT(id,ddigi);
      unsigned maxFrequency = complexTransform.size();
      unsigned maxN = gTransform[id]->GetN();
      double xpoint,ypoint;
      for(unsigned iw =1 ; iw<maxFrequency; ++iw) {
        int ipoint=iw-1;
        if(iw> maxFrequency/2) ipoint = maxFrequency-iw-1;
        gTransform[id]->GetPoint(ipoint,xpoint,ypoint);
        complexTransform[iw] *= ypoint;
      }
      std::vector< Double_t > fdigi = inverseFFT(id,complexTransform,ddigi);
      for(unsigned isample = 0; isample < fdigi.size(); isample++) {
        hEvFiltWave[id]->SetBinContent(isample+1,fdigi[isample]);
      }

      // use filtered wave
      if(pmtMode==0) ddigi=fdigi;

      // FFT of filtered wave
      std::vector<std::complex<double> > complexTransformFilt;
      complexTransformFilt = FFT(id,ddigi,false);



      // derAVe, derSigma  best to just fit Gaussian
      hDNoise[id]->Reset();
      for(unsigned isample = 0; isample < vderiv.size(); isample++) hDNoise[id]->Fill(vderiv[isample]);
      hDNoise[id]->Fit("gaus","0Q");
      TF1 *fit1 = (TF1*)hDNoise[id]->GetFunction("gaus");
      derAve = fit1->GetParameter(1);
      derSigma = fit1->GetParameter(2);

      // find  peaks 
      // for derivativePeaks, window in time is timeUnit*windowSize (ns) . timeUnit = 2
      Int_t windowSize=10;              
      // min, max width in time bins for simple peaks
      unsigned maxWidth = 100000; 
      unsigned minWidth = 10;
      peakList.clear();
      peakKind.clear();
      if(isPMT) peakList= derivativePeaks(vderiv, id, windowSize, derSigma, peakKind);
      else peakList= simplePeaks(ddigi, id, minWidth, maxWidth, rawSigma, peakKind);
      hPeakCount->Fill(id,peakList.size());

      

      //printf("\t det %i %s has %lu peaks \n",id,dname.Data(),peakList.size() );
      //for (int ip=0; ip<int(peakList.size()); ++ip) printf(" \t ip %i start %i end %i \n",ip, std::get<0>(peakList[ip]) ,std::get<1>(peakList[ip]) );

      // get WMARecursive baseline weigthts 
    

      Double_t filtAve, filtSigma;
      getAverages(ddigi, filtAve, filtSigma, prePulse); // first bins
      std::vector<Double_t> baselineDigi(ddigi.size(),filtAve);


      int maxWeightWidth=0;
      std::vector<Double_t> weight = getBaselineWeights(ddigi.size(), peakList,maxWeightWidth);
      // set pre-trigger weights to zero
      unsigned wbin = 1000;
      for (unsigned i = 0; i < wbin; i++) weight[i]=0; 
      for (unsigned i = 0; i < weight.size(); i++) hEvWeight[id]->SetBinContent(i+1, weight[i]); 

      // get WMARecursive baseline 
      /*
      if(isPMT) {
        int NWindow= 1000;  // size is 8200
        getBaselineWMARecursive(0,ddigi.size(),ddigi, weight, NWindow,baselineDigi);
      }
      // zero baseline until after first pulse
      if(isPMT&&peakList.size()>0) for(unsigned i=0; i<= std::get<1>(peakList[0]) ; ++i ) baselineDigi[i]=0;
      */

      
      // save wiggle 
      for (unsigned i = 0; i < detList[id]->digi.size(); i++) hSWaveSum[id]->SetBinContent(i+1,hSWaveSum[id]->GetBinContent(i+1)+ddigi[i] - baselineDigi[i]);

      // remove PMT wiggle
      if(pmtMode==1&&peakList.size()>0&&isPMT&&hWiggle) 
        for (unsigned i =  std::get<1>(peakList[0])+1; i < baselineBreak; i++) baselineDigi[i] += hWiggle->GetBinContent(i+1);

      /* get WMARecursive baseline for all */
      int NWindow=0;  // size is 8200
      if(isPMT)  NWindow= 1000;
      else  NWindow= 8200;

      bool doWMA =true;
      if(doWMA&&isPMT) {
        getBaselineWMARecursive(0,ddigi.size(),ddigi, weight, NWindow,baselineDigi);
        // zero baseline until after first pulse
        if(peakList.size()>0) for(unsigned i=0; i<= std::get<1>(peakList[0]) ; ++i ) baselineDigi[i]=0;
      }
      // subract base 
      sdigi.clear();
      double qDetSum=0;
      double qDetPrompt=0;
      for (unsigned i = 0; i < detList[id]->digi.size(); i++) {
        hBaselineWMA[id]->SetBinContent(i, baselineDigi[i]); 
        double sval = ddigi[i] - baselineDigi[i];
        qDetSum += vsign[id]*sval;
        if(i>720&&i<770) qDetPrompt+= vsign[id]*sval;
        sdigi.push_back(sval);
        hEvSignal[id]->SetBinContent(i+1,sval);
        hSWaveSumSubtracted[id]->SetBinContent(i+1,hSWaveSumSubtracted[id]->GetBinContent(i+1)+sval);
      }

      
      // redo peak finding for PMT if wiggle subtracted
      if(isPMT&&pmtMode==1) {
        vderiv.clear();
        vderiv = differentiate(sdigi,diffStep);
        unsigned maxWidth = 100000; 
        unsigned minWidth = 10;
        peakList.clear();
        peakKind.clear();
        peakList= derivativePeaks(vderiv, id, windowSize, derSigma, peakKind);
        //hPeakCount->Fill(id,peakList.size())
      }
      

      // extend peaks to baseline subtracted wave 
      if(isPMT&&pmtMode==0) trimPeaks(id,sdigi,peakList);
      else extendPeaks(id,sdigi,peakList);
      /*** peak finding */
      for(unsigned isample = 0; isample < sdigi.size(); isample++) hSNoise[id]->Fill(sdigi[isample]);
      hSNoise[id]->Fit("gaus","0Q");
      TF1 *fit2 = (TF1*)hSNoise[id]->GetFunction("gaus");
      double sAve = 0;
      double sSigma = 0;
      if(fit2) {
        sAve=fit2->GetParameter(1);
        sSigma=fit2->GetParameter(2);
      } else 
        printf(" no SNoise fit event %lld det %i ",ievent,id);

      ntCal->Fill(ievent,id,rawAve,rawSigma,derAve,derSigma,sAve,sSigma);

      /***********  make hits *********/
      double triggerTime=0;
      double firstCharge=0;
      int hitCount=0;
      hitMap  detHits = makeHits(id,peakList,peakKind,sdigi,sSigma,triggerTime,firstCharge);

      //printf("after makeHits detHits det %i name %s has %lu peaks \n",id,dname.Data(),detHits.size() );

      int ihit=0;
      double hitSum=0;
      double hitPrompt=0;
      /*** loop over hits ***/
      for (hitMapIter hitIter=detHits.begin(); hitIter!=detHits.end(); ++hitIter) {
        TDetHit  hiti = hitIter->second;
        Double_t hitQ = hiti.qsum;
        hitSum += hitQ;
        if(hiti.startTime>720&&hiti.startTime>770) hitPrompt += hitQ;
        brun->detList[id]->hits.push_back(hiti);
        //Double_t phitQErr = phiti.qerr*timeUnit*1E9;
        Int_t width = hiti.lastBin - hiti.firstBin +1;
        //printf(" \t ihit %i time %f start bin  %i end bin %i kind %i charge %f  \n",ihit++,hitIter->first, hiti.firstBin, hiti.lastBin, hiti.kind,hiti.qsum );

        for (int ibin=hiti.firstBin; ibin <= hiti.lastBin; ++ibin) {
          hEvPeaks[id]->SetBinContent(ibin+1,sdigi[ibin]);
          hEvDPeaks[id]->SetBinContent(ibin+1,vderiv[ibin]);
        }
        /*** PulseShape sum ***/
        bool doPulse = true;
        if(hiti.qsum<0) doPulse=false;
        // skip early PMT pulses for sum
        if(id==1&&hiti.startTime*timeUnit<3000.0) doPulse=false;
        if(doPulse) {
          int startPulse = TMath::Max(0,hiti.peakBin-500);
          int endPulse  = TMath::Min(int(sdigi.size()-1),hiti.peakBin+500);
          if(firstPulseTime==0) firstPulseTime=hiti.peakBin;
          if(firstPulseTime==hiti.peakBin||hiti.peakBin>firstPulseTime+1000)  {
            ++pulseCount[id];
            for (int ibin=startPulse; ibin <= endPulse; ++ibin) {
              int pbin = ibin+1-hiti.peakBin+1000;
              hPulseShape[id]->SetBinContent(pbin,hPulseShape[id]->GetBinContent(pbin)+sdigi[ibin]);
            }
            firstPulseTime=hiti.peakBin;
          }
        }
        /******/
        Double_t hitTime =  hiti.startTime*timeUnit*microSec;
        Int_t istartBin =  hLife[id]->FindBin(hitTime);
        ntHit->Fill(id,detHits.size(), ++hitCount, istartBin, hiti.startTime*timeUnit*microSec, hitQ, width, hiti.qpeak, hiti.good, hiti.kind);
        hLife[id]->SetBinContent(istartBin, hLife[id]->GetBinContent(istartBin)+hitQ);
        hLife[id]->SetBinError(istartBin, sqrt(pow(hLife[id]->GetBinError(istartBin), 2)+pow(hiti.qerr, 2)));
      }

      if (ievent<10) plotWave(ievent, id);

      // fill in det header variables
      brun->detList[id]->event =ievent;
      brun->detList[id]->nspe = int(peakList.size());// needs to be fixed
      brun->detList[id]->qSum=qDetSum;  
      brun->detList[id]->qPrompt=qDetPrompt;  
      brun->detList[id]->hitSum=hitSum;  
      brun->detList[id]->hitPrompt=hitPrompt; // same 

    } // end of det loop
    // fill BRun
    brun->bevent->event =ievent;
    brun->bevent->trigTime =trigTime;
    brun->fill();

    if(dumpFile.is_open() ) brun->dumpEvent(dumpFile);
  } // loop over events
}

bool  anaWRun::openFile() 
{
  // open ouput file and make some histograms
  TString fileName; fileName.Form("rootData/%s.root",runTag.Data());
  printf(" looking for file %s\n",fileName.Data());
  fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return false;
  }

  fin->GetObject("BTree",btree);
  TObjArray *brList = btree->GetListOfBranches();

  TIter next(brList);
  TBranch *aBranch=NULL;
  int idet=0;

  while( ( aBranch = (TBranch *) next() ) ) {
    detList[idet] = new TBRawEvent(aBranch->GetName());
    btree->SetBranchAddress(aBranch->GetName() , &detList[idet]);
    detList[idet]->SetName(aBranch->GetName());
    ++idet;
  }

  for(unsigned id=0; id<NDET; ++id) hEvent[id]=NULL;
  cout << " number of entries in BTree is " <<  btree->GetEntries() << endl;
  btree->GetListOfBranches()->ls();
  return true;
}

// WMARecursive weights
std::vector<Double_t> anaWRun::getBaselineWeights(unsigned arraySize, peakType peakList, Int_t& maxWeightWidth)
{
    std::vector<Double_t> weight(digiSize, 0);
    if (peakList.size() == 0) return weight;

    // construct the weights
    int alpha;
    int beta;
    maxWeightWidth=0;
    for (int ip=0; ip<int(peakList.size()); ++ip) {
        if (ip==0) beta = -1;
        else beta = std::get<1>(peakList[ip-1]);
        if (ip==peakList.size()) alpha = weight.size();
        else alpha = std::get<0>(peakList[ip]);
        int width = std::get<1>(peakList[ip])-std::get<0>(peakList[ip]);
        //printf(" \t weight %lu peak %i at %i width %i\n",weight.size(),ip,std::get<0>(peakList[ip]),width);
        if (width>maxWeightWidth) maxWeightWidth=width;
        for (int jp=beta+1; jp<alpha; ++jp)  weight[jp]=alpha-beta-1;
    }
    return weight;
}
//Weighted Moving Average baseline from Zugec et al https://arxiv.org/pdf/1601.04512.pdf 
void  anaWRun::getBaselineWMARecursive(int startBin, int stopBin,std::vector<Double_t> sig, std::vector<Double_t> weight,Int_t NWindow,std::vector<Double_t>& baseline)
{

  // starting values 
  Double_t KW=0;
  Double_t CW=0;
  Double_t SW=0;
  Double_t KS=0;
  Double_t CS=0;
  Double_t SS=0;

  // recursive sums
  Double_t kcos = TMath::Cos( TMath::Pi()/Double_t(NWindow));
  Double_t ksin = TMath::Sin( TMath::Pi()/Double_t(NWindow));
  Int_t resetBin = 100;

  for(int iw=int(startBin); iw<int(stopBin); ++iw) {
    // reset
    int low =  TMath::Max(0,iw-NWindow);
    int high = TMath::Min(iw+NWindow,stopBin-1);
    if(iw%resetBin ==  0) {
      KW=0;CW=0;SW=0;KS=0; CS=0;SS=0;
      for(int jp =low; jp < high; ++jp) {
        Double_t cosj = TMath::Cos(  Double_t(jp-iw)*TMath::Pi()/Double_t(NWindow) );
        Double_t sinj = TMath::Sin(  Double_t(jp-iw)*TMath::Pi()/Double_t(NWindow) );
        KW+= weight[jp];
        KS+= sig[jp]*weight[jp];
        CW+= weight[jp]*cosj;
        CS+= sig[jp]*weight[jp]*cosj;
        SW+= weight[jp]*sinj;
        SS+= sig[jp]*weight[jp]*sinj;
      }
    } else {
      //save previous iteraton
      int i1 = iw+NWindow;
      int i2 = iw-1-NWindow;
      Double_t CWlast=CW;
      Double_t CSlast=CS;

      // case sums 
      CW = kcos*CWlast+ksin*SW;
      CS = kcos*CSlast+ksin*SS;
      SW = kcos*SW-ksin*CWlast;
      SS = kcos*SS-ksin*CSlast;

      if(iw  <= NWindow && iw + NWindow <= stopBin - 1){
        KW += weight[i1];
        KS += sig[i1]*weight[i1];
        CW += -1.*weight[i1];
        CS += -1.*sig[i1]*weight[i1];
      } else if(iw > NWindow && iw + NWindow <= stopBin - 1){
        KW += -1.*weight[i2] + weight[i1];
        KS += -1.* sig[i2]*weight[i2] +  sig[i1]*weight[i1];
        CW +=  kcos*weight[i2] - weight[i1];
        CS +=  kcos*sig[i2]*weight[i2] - sig[i1]*weight[i1];
        SW += -1.*ksin*weight[i2] ;
        SS += -1.*ksin*sig[i2]*weight[i2] ;
      } else if(iw  > NWindow && iw + NWindow > stopBin - 1){
        KW += -1.*weight[i2];
        KS +=  -1.*sig[i2]*weight[i2];
        CW +=  kcos*weight[i2];
        CS +=  kcos*sig[i2]*weight[i2];
        SW += -1.*ksin*weight[i2];
        SS += -1.*ksin*sig[i2]*weight[i2];
      }  // sum baseline 
    }
    if( TMath::Abs(KW+CW)>1E-12) baseline[iw]=(KS+CS)/(KW+CW);
  }
  return;
}

peakType anaWRun::simplePeaks(std::vector<Double_t> v, Int_t idet, unsigned minWidth, unsigned maxWidth, Double_t rms,std::vector<Int_t>& peakKind ) 
{
  peakType peakList;
  peakKind.clear();
  std::vector<unsigned> crossings;
  std::vector<unsigned> crossingBin;
  std::vector<double> crossingTime;
  unsigned vsize = v.size();
  Double_t cut = derivativeSigma[idet]*rms;
   //cout << " >>>> det "<< idet << " rms " << rms << " cut " << cut << endl;

  /* loop over bins finding v>cut.  continue until v<cut */
  for( unsigned ibin=1; ibin< vsize-1; ++ibin ) {
    Double_t vi = vsign[idet]*v[ibin];
    Double_t vmax=0;
    Double_t vj=0;
    unsigned jbin=ibin+1;
    if(vi>cut) {
      //if(idet==2) printf(" at  %i ibin %u vi %f\n",idet,ibin,vi);
      for(jbin=ibin+1; jbin< vsize; ++jbin ) {
        vj = vsign[idet]*v[jbin];
        if(vj>vmax) vmax=vj;
        //if(idet==2) printf(" \t   jbin %u vj %f \n",jbin,vj);
        if(vj<cut) break;
      }
      if(jbin-ibin>minWidth&&jbin-ibin<maxWidth) {
        peakList.push_back( std::make_pair(ibin,jbin-1) );
        peakKind.push_back(1);
        ntSimplePeaks->Fill(idet,rms,vmax,float(ibin),float(jbin-ibin));
        //printf(" found peak det %i vi %f vj %f start bin %u end  %u \n",idet,vi,vsign[idet]*v[jbin-1],ibin,jbin-ibin);
        ibin = jbin+1;
      }
    }
  }
  // return list
  return peakList;
}

peakType anaWRun::derivativePeaks(std::vector<Double_t> v,  Int_t idet, Int_t nsum, Double_t rms,std::vector<Int_t>& peakKind ) 
{
  peakType peakList;
  peakKind.clear();
  std::vector<unsigned> crossings;
  std::vector<unsigned> crossingBin;
  std::vector<double> crossingTime;
  unsigned vsize = v.size();
  Double_t cut = derivativeSigma[idet]*rms;
  //cout << " for det " << idet << " in derivative peaks >>>> rms " << rms << " cut " << cut << endl;
  Double_t ncut = -derivativeSigma[idet]*rms;
  // find all crossings
  for( unsigned ibin=1; ibin< vsize; ++ibin ) {
    Double_t u = double(ibin)*timeUnit*microSec;
    Double_t vi = vsign[idet]*v[ibin];
    Double_t vj = vsign[idet]*v[ibin-1];
    unsigned ctype=10;
    if(vj>cut&&vi<ncut) {
      crossings.push_back(DOUBLEUPCROSS);
      ctype=DOUBLEUPCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    } else if(vj<ncut&&vi>cut) {
      crossings.push_back(DOUBLEDOWNCROSS);
      ctype=DOUBLEDOWNCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    } else if( vi>cut && vj <cut ) {
      crossings.push_back(UPCROSS);
      ctype=UPCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    } else if( vi<cut && vj > cut ) {
      crossings.push_back(UPCROSS);
      ctype=UPCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    } else if( vi<ncut && vj > ncut ) {
      crossings.push_back(DOWNCROSS);
      ctype=DOWNCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    } else if( vi >ncut && vj < ncut ) {
      crossings.push_back(DOWNCROSS);
      ctype=DOWNCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    //if (idet==1&&ctype<10)  printf("....... %u vj %f vi %f cut %f cross type %u \n", ibin, vj, vi, cut, ctype );
    //if (idet==1&&ibin>2350&&ibin<2450)  printf("\t %u vj %f vi %f ctype %u  \n", ibin, vj, vi, ctype );
  }

  if(crossings.size()<4) return peakList;

  // label found crossings, intially all false
  std::vector<bool> crossingFound;
  crossingFound.resize(crossings.size());
  for(unsigned jc=0; jc<crossings.size(); ++jc)  {
    //printf(" det %i crossing %i bin %i time %f type %i \n",idet,jc,crossingBin[jc],crossingTime[jc],crossings[jc]);
    crossingFound[jc]=false;
  }


  // parse crossings to make pairs 
  /* first find sequence of 0,2,x x>0 crossings */
  unsigned ip=0;
  while ( ip<= crossings.size() -3 ) {
    int ibin = crossingBin[ip];
     if(crossings[ip]==UPCROSS&&crossings[ip+1]==DOUBLEUPCROSS&&crossings[ip+2]>UPCROSS) {
        //if (idet==1)  printf("\t peak %lu ibin %i type (%i,%i,%i) \n", peakList.size() ,ibin, crossings[ip],crossings[ip+1],crossings[ip+2] );
        peakList.push_back( std::make_pair(crossingBin[ip],crossingBin[ip+2]) );
        peakKind.push_back(0);
        //printf(" det %i make peak  (%i,%i) kind %i  \n",idet,crossingBin[ip],crossingBin[ip+2],peakKind[peakKind.size()-1]);
        ntDer->Fill(rms,v[crossingBin[ip]],double(crossingBin[ip+2]-crossingBin[ip]),double(0));//sigma:d0:step:dstep
        crossingFound[ip]=true;
        crossingFound[ip+1]=true;
        crossingFound[ip+2]=true;
        ip=ip+3;
     } else {
       ++ip;
     }
  }


  
  /* 
  ** next find sequence of 0,0,1,1 crossings */
  ip=0;
  while ( ip<= crossings.size() -4 ) {
    int skip=0;
    if(crossingFound[ip])   ++skip;
    if(crossingFound[ip+1]) ++skip;
    if(crossingFound[ip+2]) ++skip;
    if(crossingFound[ip+3]) ++skip;
    if(skip>0) {
      ip+=skip;
      continue;
    }
    if(crossings[ip]==UPCROSS&&crossings[ip+1]==UPCROSS&&crossings[ip+2]==DOWNCROSS&&crossings[ip+3]==DOWNCROSS) {
      int ibin = crossingBin[ip];
      //if (idet==1)  printf("\t peak %lu ibin %i type (%i,%i,%i,%i) \n", peakList.size() ,ibin, crossings[ip],crossings[ip+1],crossings[ip+2],crossings[ip+3]);
      peakList.push_back( std::make_pair(crossingBin[ip],crossingBin[ip+2]) );
      peakKind.push_back(1);
      //printf(" det %i make peak  (%i,%i) kind %i  \n",idet,crossingBin[ip],crossingBin[ip+3],peakKind[peakKind.size()-1]);
      ntDer->Fill(rms,v[crossingBin[ip]],double(crossingBin[ip+2]-crossingBin[ip]),double(1));//sigma:d0:step:dstep
      crossingFound[ip]=true;
      crossingFound[ip+1]=true;
      crossingFound[ip+2]=true;
      crossingFound[ip+3]=true;
      ip=ip+4;
    } else {
      ++ip;
    }
  }



  // return list
  return peakList;
}
void anaWRun::getAverages(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma, unsigned max) 
{
  if(max==0) max = digi.size();
  std::vector<double> vsort;
  for (unsigned is=0; is< max; ++is) vsort.push_back(digi[is]); 
  std::sort(vsort.begin(), vsort.end());
  ave = vsort[unsigned(0.5*double(vsort.size()))];
  sigma = TMath::Abs(vsort[0.659*vsort.size()]-ave);
}

void anaWRun::getAverage(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma, unsigned max) 
{
  
  if(max==0) max = digi.size();
  double sum=0;
  double sum2=0;
  for (unsigned is=0; is< max; ++is) {
    sum += digi[is];
    sum2 += pow(digi[is],2.);
  }
  ave = sum/double(max);
  sigma = sqrt( sum2/double(max) - pow(ave,2.));
}
//implimented as in zugec et al, arXiv:1601.04512v1
std::vector<Double_t> anaWRun::differentiate(std::vector<Double_t> v, unsigned diffStep)
{
  std::vector<Double_t> d;
  unsigned nsamples = v.size();
  Double_t sump=0;
  Double_t summ=0;
  d.push_back(0); // first entry is zero
  for(unsigned i=1; i<nsamples; ++i) {
    unsigned i2 = 2*i;
    unsigned max = TMath::Min( diffStep, i);
    max = TMath::Min(max, nsamples - 1 -i);
    // beginning, middle, end cases
    if(i<=diffStep && i2 <= nsamples -1 ) {
      sump = sump - v[i] + v[i2-1] + v[i2];
      summ = summ + v[i-1];
    }
    else if(i>diffStep && i+diffStep <= nsamples -1){
      sump = sump - v[i] + v[i+diffStep];
      summ = summ + v[i-1]-v[i-1-diffStep];
    }
    else if(i+diffStep> nsamples-1 && i2 > nsamples-1) {
      sump = sump - v[i];
      summ = summ + v[i-1] - v[i2-nsamples-1] - v[i2-nsamples];
    }
    d.push_back(sump-summ);
  }
  return d;
}

hitMap anaWRun::makeHits(int idet, peakType peakList, std::vector<Int_t> peakKind, std::vector<Double_t> ddigi, Double_t sigma, Double_t& triggerTime, Double_t& firstCharge)
{
    triggerTime=1E9;
    firstCharge=0;
    hitMap detHits;
    if (peakList.size()<1) return detHits;
    Double_t qmax=0;

    unsigned minLength=5;
    for (unsigned ip=0; ip<peakList.size(); ++ip) {
        unsigned klow  = std::get<0>(peakList[ip]);
        unsigned khigh = std::get<1>(peakList[ip]);
        //printf(" hit  %u (%u,%u) kind %i length %u ",ip,klow,khigh,peakKind[ip],khigh-klow+1);
        hHitLength->Fill(khigh-klow+1);
        if (khigh-klow+1<minLength) {
            continue;
        }
        Double_t qhit=0;
        UInt_t peakt=0;
        Double_t qpeak=0;
        Double_t qsum = 0;
        for (unsigned k=klow; k<khigh; ++k) {
            double qdigik = vsign[idet]*ddigi[k];
            qsum += qdigik;
            if (qdigik>qpeak) {
                peakt=k;
                qpeak = qdigik;
            }
        }

        TDetHit dhit;
        dhit.peakBin=Int_t(peakt);
        dhit.qsum=qsum;
        dhit.qpeak=qpeak;
        dhit.firstBin = klow;
        dhit.lastBin = khigh;
        dhit.peakMaxTime=peakt;
        dhit.peakt=peakt;
        dhit.startTime=klow;
        dhit.peakWidth=khigh-klow;
        // this is N= q/qnorm and delta q = root(n)*qnorm;
        dhit.qerr = sqrt(pow(sigma*Double_t(dhit.peakWidth), 2)+ qnorm*qsum);
        dhit.kind = peakKind[ip];

        // just use the biggest pulse 
        if (qsum>qmax) {
            qmax=qsum;
            triggerTime=dhit.startTime*timeUnit*microSec;
            firstCharge = qsum;
        }
        Double_t hitTime = dhit.startTime*timeUnit*microSec;
        //printf("  insert hit  %lu time %f (%u,%u) kind %i length %u  \n",detHits.size(),hitTime,dhit.firstBin,dhit.lastBin, peakKind[ip],khigh-klow+1 );
        //for (unsigned k=klow; k<khigh; ++k) printf(" \t %u %f ; ", k, ddigi[k]);
        //cout << endl;
        detHits.insert(std::pair<Double_t, TDetHit>(hitTime, dhit));
        hPeakNWidth->Fill(dhit.lastBin-dhit.firstBin+1);
    }

    // first time, charge from map
    /*
    hitMapIter hitIter;
    hitIter=detHits.begin();
    TDetHit dhit0 = hitIter->second;
    triggerTime = dhit0.startTime*microSec;
    firstCharge = dhit0.qsum;
    */
    //printf(" return from makeHits with %lu made \n",detHits.size()); 
    return  detHits;
}

void anaWRun::trimPeaks(int idet, std::vector<Double_t> v, peakType &peakList ) {

  if(peakList.size()<1) return;
  for(unsigned ip=0; ip<peakList.size(); ++ip)  {
    // trim first peak 
    unsigned peakStart = std::get<0>(peakList[ip]);
    unsigned peakEnd = std::get<1>(peakList[ip]);
    for(unsigned kp= peakEnd; kp > peakStart ; --kp ) {
      double vp = vsign[idet]*v[kp];
      if( vp >0) break;
      std::get<1>(peakList[ip])=kp;
    }

    for(unsigned kp= peakStart; kp < peakEnd ; ++kp ) {
      double vp = vsign[idet]*v[kp];
      if( vp >0 ) break;
      std::get<0>(peakList[ip])=kp;
    }
  }
}

// extend peaks to zero of waveform  
void anaWRun::extendPeaks(int idet, std::vector<Double_t> v, peakType &peakList ) {
  if(peakList.size()<1) return;
  for(unsigned ip=0; ip<peakList.size(); ++ip)  {
    // high direction
    unsigned high = std::get<1>(peakList[ip]);
    unsigned next = v.size();
    if(ip<peakList.size()-1) next = std::get<1>(peakList[ip+1]);
    for(unsigned kp= high; kp < next ; ++kp ) {
      double vp = vsign[idet]*v[kp];
      if( vp <0) break;
      std::get<1>(peakList[ip])=kp;
    }

    // low direction
    unsigned low = std::get<0>(peakList[ip]);
    unsigned prev = 1;
    if(ip>0) prev = TMath::Max(prev,std::get<0>(peakList[ip-1]));
    for(unsigned kp= low; kp > prev ; --kp ) {
      double vp = vsign[idet]*v[kp];
      if( vp <0 ) break;
      std::get<0>(peakList[ip])=kp;
    }
    //printf("\t  extend peak %u from (%u,%u) to  (%u,%u)  \n",ip,low,high,std::get<0>(peakList[ip]),std::get<1>(peakList[ip])    );
  }
}



void anaWRun::plotWave(Long64_t ievent, unsigned idet) {
    evDir->cd();
    TString histName;
    TString detName = detList[idet]->GetName();

    
    histName.Form("EvRawWave%lli_DET%1i_%s", ievent,idet, detName.Data());
    TH1D* hrawwave = (TH1D*)hEvRawWave[idet]->Clone(histName);

    histName.Form("EvFiltWave%lli_DET%1i_%s", ievent,idet, detName.Data());
    TH1D* hfiltwave = (TH1D*)hEvFiltWave[idet]->Clone(histName);


    histName.Form("EvWave%lli_DET%1i_%s", ievent,idet, detName.Data());
    TH1D* hwave = (TH1D*)hEvWave[idet]->Clone(histName);

    histName.Form("EvFFT%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hfft = (TH1D*)hFFT[idet]->Clone(histName);

    histName.Form("EvFFTFilt%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hfftfilt = (TH1D*)hFFTFilt[idet]->Clone(histName);


    histName.Form("EvInvFFT%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hinvfft = (TH1D*)hInvFFT[idet]->Clone(histName);

    histName.Form("EvBase%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hbase = (TH1D*)hBaselineWMA[idet]->Clone(histName);

    histName.Form("EvDWave%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hdwave = (TH1D*)hEvDWave[idet]->Clone(histName);

    histName.Form("EvSignal%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hsignal = (TH1D*)hEvSignal[idet]->Clone(histName);

    histName.Form("EvPeaks%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hpeaks = (TH1D*)hEvPeaks[idet]->Clone(histName);

    histName.Form("EvDPeaks%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hdpeaks = (TH1D*)hEvDPeaks[idet]->Clone(histName);

    histName.Form("EvWeight%lli_DET%1i_%s", ievent,idet,detName.Data());
    TH1D* hweight = (TH1D*)hEvWeight[idet]->Clone(histName);

    fout->cd();
}
