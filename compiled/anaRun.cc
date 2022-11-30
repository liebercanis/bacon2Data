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
#include <numeric>
#include <algorithm> // std::sort

#include "TBRawEvent.hxx"
#include "TBRun.hxx"
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TBranchElement.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
//
#include "TBRun.hxx"
#include "hitFinder.hxx"

class anaRun
{
public:
  TFile *fout;
  TFile *fin;
  TTree *rtree;
  TBRun *tbrun;
  TObjArray *brList;
  vector<TBRawEvent *> detList;
  vector<TH1D *> noiseHist;
  vector<TH1D *> skewHist;
  vector<TH1D *> sumWave;
  vector<TH1D *> valHist;
  TH1D *hEvRawWave;
  TH1D *hEvGaus;
  vector<int> vchan;
  ofstream dumpFile;
  vector<TBWave *> waveList;
  hitFinder *finder;
  TString tag;
  anaRun(const char *theTag = "run", Long64_t maxEntries = 0);
  virtual ~anaRun() { ; }
  bool openFile();
  long getEvent(Long64_t entry);
  void anaEvent(Long64_t entry);
  Long64_t nentries;
  int nfill3;
};

/* get detList */
long anaRun::getEvent(Long64_t entry)
{
  rtree->GetEntry(entry);
  // find branches and cast as TBRawEvent. each is a channel
  // cout << " getEvent " << entry << endl;
  // brList->ls();
  TIter next(brList);
  TBranchElement *aBranch = NULL;
  detList.clear();
  while ((aBranch = (TBranchElement *)next()))
  {
    TBRawEvent *det = (TBRawEvent *)aBranch->GetObject();
    detList.push_back(det);
  }
  return detList.size();
}

/* analyze detList */
void anaRun::anaEvent(Long64_t entry)
{
  // loop over channels
  for (int ichan = 0; ichan < detList.size(); ++ichan)
  {
    int nbins = detList[ichan]->rdigi.size();
    if (nbins < 1)
      continue;

    // simple baseline
    //double base = std::accumulate(detList[ichan]->rdigi.begin(), detList[ichan]->rdigi.end(), 0) / detList[ichan]->rdigi.size();
    unsigned  baseLength = 100;

    double base = 0;
    for (unsigned j = 0; j < baseLength; ++j){
      base += (double) detList[ichan]->rdigi[j];
    }
    base /= double(baseLength);

    TString hname;
    hEvRawWave = NULL;

    // baseline correction from fitted Gaussian
    hEvGaus->Reset();
    for (unsigned j = 0; j < detList[ichan]->rdigi.size(); ++j)
      hEvGaus->Fill((double)detList[ichan]->rdigi[j] - base);

    hEvGaus->Fit("gaus", "QO");
    TF1 *gfit = (TF1 *)hEvGaus->GetListOfFunctions()->FindObject("gaus");
    double ave = gfit->GetParameter(1);
    double sigma = gfit->GetParameter(2);
    double skew = hEvGaus->GetSkewness();
    // fitptr->Print();
    // cout << hname << " " <<   " ave " << ave << " sigma " << sigma << " skew " << skew <<  endl;
    noiseHist[ichan]->Fill(sigma);
    skewHist[ichan]->Fill(skew);
    double sign = TMath::Sign(1., skew);
    // skip events with no hits
    if (abs(skew) < 1)
      continue;

    if(ichan==3)
      ++nfill3;

    // plot some single events
    if (entry < 100)
    {
      hname.Form("EvRawWaveEv%ich%i", int(entry), ichan);
      hEvRawWave = new TH1D(hname, hname, nbins, 0, nbins);
      for (unsigned j = 0; j < detList[ichan]->rdigi.size(); ++j)
      {
        //double val = sign * ((double)detList[ichan]->rdigi[j] - base - ave - valHist[ichan]->GetMean());
        double val = sign * ((double)detList[ichan]->rdigi[j] - base);
        hEvRawWave->SetBinContent(j + 1, val);
      }
    }

    // summed plots
    for (unsigned j = 0; j < detList[ichan]->rdigi.size(); ++j)
    {
      double val = sign * ((double)detList[ichan]->rdigi[j] - base);
      //double val = sign * ((double)detList[ichan]->rdigi[j] - base - ave) - valHist[ichan]->GetMean();
      sumWave[ichan]->SetBinContent(j + 1, sumWave[ichan]->GetBinContent(j + 1) + val);
      valHist[ichan]->Fill(val);
    }
  }
  /*
  printf("xxxxx %i  bins %i rdigi %lu sumwave %f val %f \n", nfill3, sumWave[3]->GetNbinsX(), detList[3]->rdigi.size(),
         double(sumWave[3]->GetEntries()) / double(sumWave[3]->GetNbinsX()),
          valHist[3]->GetEntries() / double(detList[3]->rdigi.size()));
          */
}

bool anaRun::openFile()
{
  // open input file and make some histograms
  TString fileName;
  fileName.Form("data/rootData/%s.root", tag.Data());
  printf(" looking for file %s\n", fileName.Data());
  fin = new TFile(fileName, "readonly");
  if (fin->IsZombie())
   {
    printf(" couldnt open file %s\n", fileName.Data());
     return false;
    }

  rtree = NULL;
  brList = NULL;
  fin->GetObject("RawTree", rtree);
  if (!rtree)
  {
    printf("!! RawTree not found \n");
    return false;
  }

  printf(" RawTree entries = %llu\n", rtree->GetEntries());
  brList = rtree->GetListOfBranches();
  brList->ls();

  TIter next(brList);
  TBranchElement *aBranch = NULL;
  detList.clear();

  // set branches
  while ((aBranch = (TBranchElement *)next()))
  {
    TString brname = TString(aBranch->GetName());
    int chan = TString(brname(brname.First('n') + 1, brname.Length() - brname.First('n'))).Atoi();
     vchan.push_back(chan);
     aBranch->SetAddress(0);
 }
 printf("total channels  =  %lu \n", vchan.size());

  // make histograms on output file
  fout->cd();
  for (unsigned ichan = 0; ichan < vchan.size(); ++ichan)
  {
    noiseHist.push_back(new TH1D(Form("NoiseChan%i", ichan), Form("NoiseChan%i", ichan), 100, 0, 10));
    skewHist.push_back(new TH1D(Form("SkewChan%i", ichan), Form("SkewChan%i", ichan), 100, -10, 10));
    sumWave.push_back(new TH1D(Form("sumWave%i", ichan), Form("sumWave%i", ichan), 4100, 0, 4100));
    valHist.push_back(new TH1D(Form("ValChan%i", ichan), Form("ValChan%i", ichan), 1100, -20, 200));
   }
  hEvGaus = new TH1D("baseline", "baseline", 200, -100, 100);

  return true;
}

anaRun::anaRun(const char *theTag, Long64_t maxEntries)
{
  nfill3 = 0;
  tag = TString(theTag);
  cout << " starting anaRun entries = " << maxEntries << " tag =  " << tag << endl;
  // do not create tbrun TTree on output file
  tbrun = new TBRun(tag);
  fout = new TFile(Form("anaRun-%s-%llu.root", tag.Data(), maxEntries), "recreate");
  cout << " opened output file " << fout->GetName() << endl;
  fout->Append(tbrun);
  if (!openFile())
  {
    printf("cannot open file!\n");
    return;
  }
  Long64_t nentries = rtree->GetEntries();
  if (maxEntries > 0)
    nentries = TMath::Min(maxEntries, nentries);
  tbrun->nevents = nentries;
  printf("... total entries  %llu looping over %llu \n ", rtree->GetEntries(), nentries);

  for (Long64_t entry = 0; entry < nentries; ++entry)
  {
    if (entry / 1000 * 1000 == entry)
      printf("... %llu \n", entry);
    getEvent(entry);
    anaEvent(entry);
  }

  // finder->hPeakCount->Print("all");
  printf(" tbrun events %llu \n",tbrun->nevents);
  fout->ls();
  fout->Write();
  fout->Close();
}
