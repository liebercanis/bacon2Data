////////////////////////////////////////////////////////
//  M.Gold May 2022
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
#include "TBRun.hxx"

typedef std::vector<std::pair<unsigned, unsigned>> peakType;
typedef std::vector<std::pair<unsigned, unsigned>>::iterator peakTypeIter;
typedef std::map<Double_t, TDetHit, std::less<Double_t>> hitMap;
typedef std::map<Double_t, TDetHit, std::less<Double_t>>::iterator hitMapIter;

const Double_t qnorm = 1.0;

class anaRun
{
  public:
    enum
    {
      UPCROSS,
      DOWNCROSS,
      DOUBLEUPCROSS,
      DOUBLEDOWNCROSS
    };
    enum
    {
      NDET = 4
    };

    double vsign[NDET] = {1.,1.,1.,1.};
    ofstream dumpFile;
    TFile *fout;
    TBRun *brun;
    TString runTag;
    vector<TBWave*> waveList;
    TString tag;
    TBRun* tbrun;
    anaRun(const char* theTag="run");
    virtual ~anaRun() { ; }
    void plotRawWave(int ifile,  Long64_t jentry);
    void findHits(int ifile, Long64_t jentry);
    void differentiate( unsigned diffStep);
    std::vector<double> digi;
    std::vector<double> ddigi;
    peakType peakList;
    std::vector<Int_t> peakKind;
    double timeUnit;
    double microSec;
    void derivativePeaks(Int_t idet, Double_t rms);
    hitMap makeHits(int idet, Double_t &triggerTime, Double_t &firstCharge);
    void trimPeaks(int idet, std::vector<Double_t> v);
    void extendPeaks(int idet, std::vector<Double_t> v);
    // hist
    TH1D *hPeakCount;
    TH1I *hHitLength;
    TH1I *hPeakNWidth;
    };

void anaRun::plotRawWave(int ifile, Long64_t jentry) 
{
  TString hname;
  hname.Form("det-%i-event-%lli",ifile,jentry);
  Long64_t nsamples =  waveList[ifile]->Samples->GetSize();
  TH1S* hist = new TH1S(hname,hname,nsamples,0,nsamples);
  for(int i=0; i<nsamples; ++i) hist->SetBinContent(i+1,hist->GetBinContent(i+1) + waveList[ifile]->Samples->GetAt(i));
  TCanvas *can = new TCanvas(hname,hname);
  hist->Draw();
  can->Print(".gif");
}

anaRun::anaRun(const char* theTag)
{
  microSec=1.0E-3;
  timeUnit = 1.0 ; // ns per count

  fout = new TFile("anaRun.root","recreate");
  hPeakCount = new TH1D("PeakCount", " peaks by det ", NDET, 0, NDET);
  hHitLength = new TH1I("HitLength", " hit length", 100, 0, 100);
  hPeakNWidth = new TH1I("PeakNWidth", "PeakNWidth", 100, 0, 100);

  TString dirName;
  tag = TString(theTag);
  dirName.Form("data/XenonDoped1_5PPMSiPM4/DAQ/%s/RAW",tag.Data());
  //dirName.Form("data/SiPM_1_inside/DAQ/%s/RAW",tag.Data());
  cout << tag  << "  dirName " << dirName << endl;
  TSystemDirectory dir("RAW",dirName);
  TList *files = dir.GetListOfFiles();
  tbrun = new TBRun(tag);
  cout << " create TBRun  " << tbrun->GetName() << endl;

  TIter next(files);
  TSystemFile *file;
  unsigned idet = 0; // detectors sequentially numbered
  while ((file = (TSystemFile *)next()))
  {
    string name = string(file->GetName());
    cout << name << endl;
    string exten  = name.substr( name.find_last_of(".")+1 );
    if(exten!=string("root")) continue;
    string tag = name.substr( 0, name.find_last_of(".") );
    //string chan = tag.substr( tag.find_first_of("CH"), tag.find_first_of("@") - tag.find_first_of("CH"));
    string chan = tag.substr( tag.find_first_of("det"), tag.find_first_of(".root") - tag.find_first_of("det"));
    //cout << " tag " << tag << "  " << tag.find_first_of("CH") <<  "   " << " chan  " << chan << endl;
    string fullName =  string( dirName.Data())  + string("/")+name;
    cout << " tag " << tag << " open " << fullName << endl;
    TFile *f = new TFile(fullName.c_str(),"readonly");
    if(!f) continue;
    TTree *dtree=NULL;
    f->GetObject("Data_R",dtree);
    if(!dtree) continue;
    TBranch *b_Samples = dtree->FindBranch("Samples");
    if(!b_Samples) continue;
    if(b_Samples->GetEntries()<1) {
      cout << " skipping  " << tag << endl;
      continue;
    }
    cout << f->GetName()  << "  " << b_Samples->GetEntries() << endl;
    TBWave* wave = new TBWave(tag.c_str());
    TTree *tree = NULL;
    f->GetObject("Data_R",tree);
    wave->Init(tree);
    waveList.push_back(wave);
    cout << " file " << wave->GetName() << " open " << waveList.size() << endl;
    tbrun->addDet(Form("DET%i",idet++), wave->GetName()) ;
  }

  cout << " got " << waveList.size() << " files " <<  endl;

  if( waveList.size()<1) {   
    fout->Close();
    return;
  }
 
  for (unsigned ifile = 0; ifile < waveList.size(); ++ifile)
  {
    Long64_t nentries = waveList[ifile]->fChain->GetEntries();
    cout << ifile << " for wave  " << waveList[ifile]->GetName() << " det entries  " << nentries  << endl;
  }

  for (unsigned ifile = 0; ifile < waveList.size(); ++ifile)
    {
      cout << ifile << " for wave  " << waveList[ifile]->GetName() << " det " << tbrun->detList[ifile]->GetName() << endl;
      Long64_t nbytes = 0, nb = 0;
      Long64_t ientry = waveList[ifile]->LoadTree(0);
      nb = waveList[ifile]->fChain->GetEntry(0);
      Long64_t nsamples = waveList[ifile]->Samples->GetSize();

      fout->cd();

      Long64_t nentries = waveList[ifile]->fChain->GetEntries();
      cout << " entries " << nentries << endl;
      nentries = 1000;
      for (Long64_t jentry = 0; jentry < nentries; jentry++)
      {
        if(jentry/1000*1000==jentry)
          printf("\t ... %lld \n", jentry);
        ientry = waveList[ifile]->LoadTree(jentry);
        if (ientry < 0)
          break;
        nb = waveList[ifile]->GetEntry(jentry);
        nbytes += nb;
        nsamples = waveList[ifile]->Samples->GetSize();
        findHits(ifile, jentry);
        double triggerTime = 0;
        double firstCharge = 0;
        hitMap detHits = makeHits(ifile, triggerTime, firstCharge);
        // fill hits
        cout << " << detHits size " << detHits.size() << endl;

        for (hitMapIter hitIter = detHits.begin(); hitIter != detHits.end(); ++hitIter)
        {
          TDetHit hiti = hitIter->second;
          tbrun->detList[ifile]->hitSum += hiti.qsum;
          if (hiti.startTime > 720 && hiti.startTime > 770)
            tbrun->detList[ifile]->hitPrompt += hiti.qsum;
          //
          tbrun->detList[ifile]->hits.push_back(hiti);
          // Double_t phitQErr = phiti.qerr*timeUnit*1E9;
        }
        //
        tbrun->fill();
        if (jentry < 4)
            plotRawWave(ifile, jentry);
      }
    }
    for (int i = 0; i < tbrun->detList.size(); ++i)
      printf(" det %i %s \n", i, tbrun->detList[i]->GetName());
    // fout->ls();
    hPeakCount->Print("all");
    fout->Write();
    fout->Close();
  }

void anaRun::findHits(int ifile, Long64_t jentry) {

  // fill vector
  Long64_t nsamples =  waveList[ifile]->Samples->GetSize();
  digi.clear();
  for(int i=0; i<nsamples; ++i) digi.push_back( double(waveList[ifile]->Samples->GetAt(i)) );

  // get averages
  std::vector<double>vsort;
  for (unsigned is=0; is< nsamples; ++is) vsort.push_back(digi[is]); 
  std::sort(vsort.begin(), vsort.end());
  double ave = vsort[unsigned(0.5*double(vsort.size()))];
  double sigma = TMath::Abs(vsort[0.659*vsort.size()]-ave);
  tbrun->detList[ifile]->event=jentry;
  tbrun->detList[ifile]->ave=ave;
  tbrun->detList[ifile]->sigma=sigma;

  unsigned diffStep = 7;
  differentiate(diffStep);

  // find peaks
          // for derivativePeaks, window in time is timeUnit*windowSize (ns) . timeUnit = 2
  // min, max width in time bins for simple peaks
  Int_t windowSize = 10;
  unsigned maxWidth = 100000;
  unsigned minWidth = 10;
  derivativePeaks(ifile, 5.0);
  hPeakCount->Fill(ifile, peakList.size());

}

/*
   void anaRun::getAverage(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma, unsigned max) 
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
*/
void anaRun::differentiate(unsigned diffStep)
{
  ddigi.clear();
  unsigned nsamples = digi.size();
  Double_t sump = 0;
  Double_t summ = 0;
  ddigi.push_back(0); // first entry is zero
  for (unsigned i = 1; i < nsamples; ++i)
  {
    unsigned i2 = 2 * i;
    unsigned max = TMath::Min(diffStep, i);
    max = TMath::Min(max, nsamples - 1 - i);
    // beginning, middle, end cases
    if (i <= diffStep && i2 <= nsamples - 1)
    {
      sump = sump - digi[i] + digi[i2 - 1] + digi[i2];
      summ = summ + digi[i - 1];
    }
    else if (i > diffStep && i + diffStep <= nsamples - 1)
    {
      sump = sump - digi[i] + digi[i + diffStep];
      summ = summ + digi[i - 1] - digi[i - 1 - diffStep];
    }
    else if (i + diffStep > nsamples - 1 && i2 > nsamples - 1)
    {
      sump = sump - digi[i];
      summ = summ + digi[i - 1] - digi[i2 - nsamples - 1] - digi[i2 - nsamples];
    }
    ddigi.push_back(sump - summ);
  }
}

void anaRun::derivativePeaks(Int_t idet, Double_t rms)
{
  peakList.clear();
  peakKind.clear();
  std::vector<unsigned> crossings;
  std::vector<unsigned> crossingBin;
  std::vector<double> crossingTime;
  unsigned vsize = ddigi.size();
  Double_t cut = tbrun->detList[idet]->sigma * rms;
  // cout << " for det " << idet << " in derivative peaks >>>> rms " << rms << " cut " << cut << endl;
  Double_t ncut = -cut;
  // find all crossings
  for (unsigned ibin = 1; ibin < vsize; ++ibin)
  {
    Double_t u = double(ibin) * timeUnit * microSec;
    Double_t vi = vsign[idet] * ddigi[ibin];
    Double_t vj = vsign[idet] * ddigi[ibin - 1];
    unsigned ctype = 10;
    if (vj > cut && vi < ncut)
    {
      crossings.push_back(DOUBLEUPCROSS);
      ctype = DOUBLEUPCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    else if (vj < ncut && vi > cut)
    {
      crossings.push_back(DOUBLEDOWNCROSS);
      ctype = DOUBLEDOWNCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    else if (vi > cut && vj < cut)
    {
      crossings.push_back(UPCROSS);
      ctype = UPCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    else if (vi < cut && vj > cut)
    {
      crossings.push_back(UPCROSS);
      ctype = UPCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    else if (vi < ncut && vj > ncut)
    {
      crossings.push_back(DOWNCROSS);
      ctype = DOWNCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    else if (vi > ncut && vj < ncut)
    {
      crossings.push_back(DOWNCROSS);
      ctype = DOWNCROSS;
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
    // if (idet==1&&ctype<10)  printf("....... %u vj %f vi %f cut %f cross type %u \n", ibin, vj, vi, cut, ctype );
    // if (idet==1&&ibin>2350&&ibin<2450)  printf("\t %u vj %f vi %f ctype %u  \n", ibin, vj, vi, ctype );
  }

  if (crossings.size() < 4)
    return;
  ;

  // label found crossings, intially all false
  std::vector<bool> crossingFound;
  crossingFound.resize(crossings.size());
  for (unsigned jc = 0; jc < crossings.size(); ++jc)
  {
    // printf(" det %i crossing %i bin %i time %f type %i \n",idet,jc,crossingBin[jc],crossingTime[jc],crossings[jc]);
    crossingFound[jc] = false;
  }

  // parse crossings to make pairs
  /* first find sequence of 0,2,x x>0 crossings */
  unsigned ip = 0;
  while (ip <= crossings.size() - 3)
  {
    int ibin = crossingBin[ip];
    if (crossings[ip] == UPCROSS && crossings[ip + 1] == DOUBLEUPCROSS && crossings[ip + 2] > UPCROSS)
    {
      // if (idet==1)  printf("\t peak %lu ibin %i type (%i,%i,%i) \n", peakList.size() ,ibin, crossings[ip],crossings[ip+1],crossings[ip+2] );
      peakList.push_back(std::make_pair(crossingBin[ip], crossingBin[ip + 2]));
      peakKind.push_back(0);
      /*
      std::pair<unsigned,unsigned>  pp = std::make_pair(crossingBin[ip], crossingBin[ip + 2]);
      unsigned ip = peakList.size() - 1;
      printf(" derivativePeaks ip = %lu \n", peakList.size() );
      unsigned klow = pp.first;
      unsigned khigh = pp.second;
      printf(" derivativePeaks   (%u,%u) \n ", klow, khigh);
      */
      // printf(" det %i make peak  (%i,%i) kind %i  \n",idet,crossingBin[ip],crossingBin[ip+2],peakKind[peakKind.size()-1]);
      //ntDer->Fill(rms, v[crossingBin[ip]], double(crossingBin[ip + 2] - crossingBin[ip]), double(0)); // sigma:d0:step:dstep
      crossingFound[ip] = true;
      crossingFound[ip + 1] = true;
      crossingFound[ip + 2] = true;
      ip = ip + 3;
    }
    else
    {
      ++ip;
    }
  }
  return;
}

hitMap anaRun::makeHits(int idet,Double_t &triggerTime, Double_t &firstCharge)
{
  double sigma = tbrun->detList[idet]->sigma;
  triggerTime = 1E9;
  firstCharge = 0;
  hitMap detHits;
  if (peakList.size() < 1)
    return detHits;
  Double_t qmax = 0;

  unsigned minLength = 5;
  cout << " makeHits peak size " << peakList.size() << " peakKind  " << peakKind.size() << endl;
  if (peakList.size()<1)
    return detHits;
  for (unsigned ip = 0; ip < peakList.size(); ++ip)
  {
    cout << "\t....  " << ip  << endl;
    unsigned klow = std::get<0>(peakList[ip]);
    unsigned khigh = std::get<1>(peakList[ip]);
    printf(" hit  %u (%u,%u) kind %i length %u ",ip,klow,khigh,peakKind[ip],khigh-klow+1);
    hHitLength->Fill(khigh - klow + 1);
    if (khigh - klow + 1 < minLength)
    {
      continue;
    }
    Double_t qhit = 0;
    UInt_t peakt = 0;
    Double_t qpeak = 0;
    Double_t qsum = 0;
    for (unsigned k = klow; k < khigh; ++k)
    {
      double qdigik = vsign[idet] * ddigi[k];
      qsum += qdigik;
      if (qdigik > qpeak)
      {
        peakt = k;
        qpeak = qdigik;
      }
    }

    TDetHit dhit;
    dhit.peakBin = Int_t(peakt);
    dhit.qsum = qsum;
    dhit.qpeak = qpeak;
    dhit.firstBin = klow;
    dhit.lastBin = khigh;
    dhit.peakMaxTime = peakt;
    dhit.peakt = peakt;
    dhit.startTime = klow;
    dhit.peakWidth = khigh - klow;
    // this is N= q/qnorm and delta q = root(n)*qnorm;
    dhit.qerr = sqrt(pow(sigma * Double_t(dhit.peakWidth), 2) + qnorm * qsum);
    dhit.kind = peakKind[ip];

    // just use the biggest pulse
    if (qsum > qmax)
    {
      qmax = qsum;
      triggerTime = dhit.startTime * timeUnit * microSec;
      firstCharge = qsum;
    }
    Double_t hitTime = dhit.startTime * timeUnit * microSec;
    // printf("  insert hit  %lu time %f (%u,%u) kind %i length %u  \n",detHits.size(),hitTime,dhit.firstBin,dhit.lastBin, peakKind[ip],khigh-klow+1 );
    // for (unsigned k=klow; k<khigh; ++k) printf(" \t %u %f ; ", k, ddigi[k]);
    // cout << endl;
    detHits.insert(std::pair<Double_t, TDetHit>(hitTime, dhit));
    hPeakNWidth->Fill(dhit.lastBin - dhit.firstBin + 1);
    }

  // first time, charge from map
  /*
  hitMapIter hitIter;
  hitIter=detHits.begin();
  TDetHit dhit0 = hitIter->second;
  triggerTime = dhit0.startTime*microSec;
  firstCharge = dhit0.qsum;
  */
  // printf(" return from makeHits with %lu made \n",detHits.size());
  return detHits;
}

void anaRun::trimPeaks(int idet, std::vector<Double_t> v)
{

  if (peakList.size() < 1)
    return;
  for (unsigned ip = 0; ip < peakList.size(); ++ip)
  {
    // trim first peak
    unsigned peakStart = std::get<0>(peakList[ip]);
    unsigned peakEnd = std::get<1>(peakList[ip]);
    for (unsigned kp = peakEnd; kp > peakStart; --kp)
    {
      double vp = vsign[idet] * v[kp];
      if (vp > 0)
        break;
      std::get<1>(peakList[ip]) = kp;
    }

    for (unsigned kp = peakStart; kp < peakEnd; ++kp)
    {
      double vp = vsign[idet] * v[kp];
      if (vp > 0)
        break;
      std::get<0>(peakList[ip]) = kp;
    }
  }
}

// extend peaks to zero of waveform
void anaRun::extendPeaks(int idet, std::vector<Double_t> v)
{
  if (peakList.size() < 1)
    return;
  for (unsigned ip = 0; ip < peakList.size(); ++ip)
  {
    // high direction
    unsigned high = std::get<1>(peakList[ip]);
    unsigned next = v.size();
    if (ip < peakList.size() - 1)
      next = std::get<1>(peakList[ip + 1]);
    for (unsigned kp = high; kp < next; ++kp)
    {
      double vp = vsign[idet] * v[kp];
      if (vp < 0)
        break;
      std::get<1>(peakList[ip]) = kp;
    }

    // low direction
    unsigned low = std::get<0>(peakList[ip]);
    unsigned prev = 1;
    if (ip > 0)
      prev = TMath::Max(prev, std::get<0>(peakList[ip - 1]));
    for (unsigned kp = low; kp > prev; --kp)
    {
      double vp = vsign[idet] * v[kp];
      if (vp < 0)
        break;
      std::get<0>(peakList[ip]) = kp;
    }
    // printf("\t  extend peak %u from (%u,%u) to  (%u,%u)  \n",ip,low,high,std::get<0>(peakList[ip]),std::get<1>(peakList[ip])    );
  }
}