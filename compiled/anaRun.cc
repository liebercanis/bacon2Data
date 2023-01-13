//////////////////////////////////////////////////////////
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
// root 
#include <TROOT.h>
#include <TKey.h>
#include <TBranch.h>
#include <TBranchElement.h>
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
// bobj classes 
#include "TBWave.hxx"
#include "TBRun.hxx"
#include "TBRawEvent.hxx"
#include "hitFinder.hxx"

class anaRun
{
public:
  TFile *fout;
  TFile *fin;
  vector<TTree *> treeList;
  vector<int> chanList;
  vector<TBRawEvent *> rawBr;
  TBRun *tbrun;
  vector<TH1D *> noiseHist;
  vector<TH1D *> skewHist;
  vector<TH1D *> sumWave;
  vector<TH1D *> valHist;
  TH1D *hEvRawWave;
  vector<TH1D *> hEvGaus;
  vector<double> ddigi;
  vector<double> hdigi;
  ofstream dumpFile;
  //vector<TBWave *> waveList;
  hitFinder *finder;
  TString tag;
  anaRun(const char *theTag = "run", Long64_t maxEntries = 0);
  virtual ~anaRun() { ; }
  bool openFile();
  unsigned getTrees();
  void getBranches();
  void getEvent(Long64_t entry);
  void anaEvent(Long64_t entry);
  Long64_t nentries;
  int nfill3;
};

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

  printf(" opened file %s\n", fileName.Data());
  return true;
}

unsigned  anaRun::getTrees()
{
  treeList.clear();
  chanList.clear();
  //
  TIter next(fin->GetListOfKeys());
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TTree"))
      continue;
    TTree *tree = (TTree *)key->ReadObj();
    //cout << tree->GetName() << " " << key->GetCycle() << endl;
    if (tree->GetEntries() < 1)
      continue;
    // pick up only last in cycle
    int ichan = TString(TString(tree->GetName())(5, 2)).Atoi();
    if (find(chanList.begin(), chanList.end(), ichan) != chanList.end())
      continue;
    cout << tree->GetName() << " chan " << ichan << endl;
    treeList.push_back(tree);
    chanList.push_back(ichan);
  }
  return treeList.size();
}

/* get rawBr */
void anaRun::getBranches()
{
  for (unsigned it = 0; it < treeList.size(); ++it)
  {
    TObjArray *brList = treeList[it]->GetListOfBranches();
    TIter next(brList);
    TBranchElement *aBranch = NULL;
    while ((aBranch = (TBranchElement *)next()))
    {
      TBRawEvent *chan = (TBRawEvent *)aBranch->GetObject();
      rawBr[it] = chan;
    }
  }
}

void anaRun::getEvent(Long64_t entry)
{
  for (unsigned it = 0; it < treeList.size(); ++it)
  {
    treeList[it]->GetEntry(entry);
    getBranches();
  }
}

/* analyze rawBr */
void anaRun::anaEvent(Long64_t entry)
{
  tbrun->clear();
  // loop over channels
  for (int ichan = 0; ichan < rawBr.size(); ++ichan)
  {
    int nbins = rawBr[ichan]->rdigi.size();
    if (nbins < 1)
      continue;

    // simple baseline
    //double base = std::accumulate(rawBr[ichan]->rdigi.begin(), rawBr[ichan]->rdigi.end(), 0) / rawBr[ichan]->rdigi.size();
    unsigned  baseLength = 100;

    double base = 0;
    for (unsigned j = 0; j < baseLength; ++j){
      base += (double) rawBr[ichan]->rdigi[j];
    }
    base /= double(baseLength);

    TString hname;
    hEvRawWave = NULL;

    // baseline correction from fitted Gaussian
    hEvGaus[ichan]->Reset();
    for (unsigned j = 0; j < rawBr[ichan]->rdigi.size(); ++j)
      hEvGaus[ichan]->Fill((double)rawBr[ichan]->rdigi[j] - base);

    hEvGaus[ichan]->Fit("gaus", "QO");
    TF1 *gfit = (TF1 *)hEvGaus[ichan]->GetListOfFunctions()->FindObject("gaus");
    double ave = gfit->GetParameter(1);
    double sigma = gfit->GetParameter(2);
    double skew = hEvGaus[ichan]->GetSkewness();
    // fitptr->Print();
    cout << hname << " " <<   " ave " << ave << " sigma " << sigma << " skew " << skew <<  endl;
    noiseHist[ichan]->Fill(sigma);
    skewHist[ichan]->Fill(skew);
    double sign = TMath::Sign(1., skew);
    TDet *idet = tbrun->getDet(ichan);
    if(idet!=NULL){
      idet->ave = ave;
      idet->sigma = sigma;
      idet->skew = skew;
      idet->event = entry;
      idet->base = base;
    }
    // fi
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
      for (unsigned j = 0; j < rawBr[ichan]->rdigi.size(); ++j)
      {
        //double val = sign * ((double)rawBr[ichan]->rdigi[j] - base - ave - valHist[ichan]->GetMean());
        double val = (double)rawBr[ichan]->rdigi[j] - base;
        hEvRawWave->SetBinContent(j + 1, val);
      }
    }

    // summed plots
    for (unsigned j = 0; j < rawBr[ichan]->rdigi.size(); ++j)
    {
      double val = sign * ((double)rawBr[ichan]->rdigi[j] - base);
      //double val = sign * ((double)rawBr[ichan]->rdigi[j] - base - ave) - valHist[ichan]->GetMean();
      sumWave[ichan]->SetBinContent(j + 1, sumWave[ichan]->GetBinContent(j + 1) + val);
      valHist[ichan]->Fill(val);
    }

    // hit finding
    /*
    ddigi.clear();
    for (unsigned j = 0; j < rawBr[ichan]->rdigi.size(); ++j)
    {
      double val = sign * ((double)rawBr[ichan]->rdigi[j] - base);
      ddigi.push_back(val);
    }

    finder->event(ichan, entry, ddigi);
    finder->plotEvent(ichan, entry);
    */
    
  } // channel loop
  // fill output Tree for each event
  tbrun->fill();
} // anaEvent


anaRun::anaRun(const char *theTag, Long64_t maxEntries)
{
  nfill3 = 0;
  tag = TString(theTag);
  cout << " starting anaRun entries = " << maxEntries << " tag =  " << tag << endl;
  if (!openFile())
  {
     printf("cannot open file!\n");
     return;
  }
  getTrees();
  cout << " treeList size " << treeList.size() << endl;
  if(treeList.size()<1)
     return;
  rawBr.resize(treeList.size());
  getBranches();

  fout = new TFile(Form("anaRun-%s-%llu.root", tag.Data(), maxEntries), "recreate");
  cout << " opened output file " << fout->GetName() << endl;
  tbrun = new TBRun(tag);
  for (unsigned it = 0; it < chanList.size(); ++it)
     tbrun->addDet(chanList[it]);

  //fout->Append(tbrun->btree);
  getEvent(0);
  for (unsigned i = 0; i< chanList.size(); ++i)
  {
     unsigned ichan = chanList[i];
     noiseHist.push_back(new TH1D(Form("noiseChan%i", ichan), Form("noiseChan%i", ichan), 100, 0, 10));
     skewHist.push_back(new TH1D(Form("skewChan%i", ichan), Form("skewChan%i", ichan), 2000, -10, 10));
     sumWave.push_back(new TH1D(Form("sumWave%i", ichan), Form("sumWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
     valHist.push_back(new TH1D(Form("valChan%i", ichan), Form("valChan%i", ichan), 1000, -200, 200));
     hEvGaus.push_back(new TH1D(Form("baseline%i", ichan), Form("baseline%i", ichan), 200, -100, 100));
  }
  //fout->ls();

  finder = new hitFinder(fout, tbrun, tag, rawBr[0]->rdigi.size(), chanList);
  Long64_t nentries = treeList[0]->GetEntries();
  if (maxEntries > 0)
    nentries = TMath::Min(maxEntries, nentries);
  printf("... total entries  %llu looping over %llu \n ", treeList[0]->GetEntries(), nentries);

  for (Long64_t entry = 0; entry < nentries; ++entry)
  {
    if (entry / 1000 * 1000 == entry)
      printf("... %llu \n", entry);
    getEvent(entry);
    /* test 
    for (unsigned it = 0; it < treeList.size(); ++it)
      printf("tree %u chan %i event %i time %lld \n", it, rawBr[it]->channel, rawBr[it]->trigger, rawBr[it]->time);
    */
    anaEvent(entry);
  }

  // finder->hPeakCount->Print("all");
  printf(" tbrun entries %llu \n",tbrun->btree->GetEntriesFast());
  tbrun->btree->GetListOfBranches()->ls();
  fout->ls();
  fout->Write();
  fout->Close();
}
