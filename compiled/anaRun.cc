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
  enum
  {
    UPCROSS,
    DOWNCROSS,
    DOUBLEUPCROSS,
    DOUBLEDOWNCROSS
  };
  TFile *fout;
  TFile *fin;
  vector<TTree *> treeList;
  vector<int> chanList;
  vector<TBRawEvent *> rawBr;
  TBRun *tbrun;
  TNtuple *ntChan;
  vector<TH1D *> noiseHist;
  vector<TH1D *> skewHist;
  vector<TH1D *> sumWave;
  vector<TH1D *> valHist;
  vector<TH1D *> threshHist;
  vector<TH1D *> crossHist;
  TH1D *hEvRawWave;
  vector<TH1D *> hEvGaus;
  TH1D *evCount;
  vector<double> digi;
  vector<double> ddigi;
  vector<double> hdigi;
  std::vector<unsigned> thresholds;
  std::vector<unsigned> crossings;
  std::vector<unsigned> crossingBin;
  std::vector<double> crossingTime;
  ofstream dumpFile;
  // vector<TBWave *> waveList;
  hitFinder *finder;
  TString tag;
  anaRun(const char *theTag = "run", Long64_t maxEntries = 0);
  virtual ~anaRun() { ; }
  bool openFile();
  unsigned getTrees();
  void getBranches();
  void getEvent(Long64_t entry);
  bool anaEvent(Long64_t entry);
  void differentiate(unsigned diffStep);
  void derivativeCount(TDet *idet, Double_t rms);
  void thresholdCount(TDet *idet, Double_t rms);
  TDirectory *badDir;
  TDirectory *evDir;
  TDirectory *sumDir;
  TDirectory *anaDir;
  Long64_t nentries;
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

unsigned anaRun::getTrees()
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
    // cout << tree->GetName() << " " << key->GetCycle() << endl;
    if (tree->GetEntries() < 1)
      continue;
    // pick up only last in cycle
    int ichan = TString(TString(tree->GetName())(5, 2)).Atoi();
    if (find(chanList.begin(), chanList.end(), ichan) != chanList.end())
      continue;
    cout << "tree name " << tree->GetName() << " chan " << ichan << endl;
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
    rawBr[it] = NULL;
    TObjArray *brList = treeList[it]->GetListOfBranches();
    TIter next(brList);
    TBranchElement *aBranch = NULL;
    while ((aBranch = (TBranchElement *)next()))
    {
      // cout << "tree " << it << " " << aBranch->GetName() << endl;
      TBRawEvent *chan = (TBRawEvent *)aBranch->GetObject();
      if (!chan)
      {
        cout << " \t\t !!!! NULL " << it << endl;
        continue;
      }
      rawBr[it] = chan;
      // cout << "rawBr name  ==>  " << rawBr[it]->GetName() << endl;
    }
  }
}

void anaRun::getEvent(Long64_t entry)
{
  for (unsigned it = 0; it < treeList.size(); ++it)
    treeList[it]->GetEntry(entry);
  getBranches();
}

/* analyze rawBr */
bool anaRun::anaEvent(Long64_t entry)
{
  tbrun->clear();
  // loop over channels
  for (int ib = 0; ib < rawBr.size(); ++ib)
  {
    unsigned ichan = chanList[ib];
    int nbins = rawBr[ib]->rdigi.size();
    // cout << ib << " nbins " << nbins << " max hist " << hEvGaus.size() << " rawBr.size() " << rawBr.size() << endl;
    if (nbins != 1024)
      continue;

    // simple baseline
    // double base = std::accumulate(rawBr[ichan]->rdigi.begin(), rawBr[ichan]->rdigi.end(), 0) / rawBr[ichan]->rdigi.size();
    unsigned baseLength = 30;

    double base = 0;
    for (unsigned j = 0; j < baseLength; ++j)
    {
      base += (double)rawBr[ib]->rdigi[j];
    }
    base /= double(baseLength);

    double peak = 0;
    unsigned peakWidth = 100;
    for (unsigned j = baseLength; j < baseLength + peakWidth; ++j)
    {
      peak += (double)rawBr[ib]->rdigi[j] - base;
    }
    // peak /= double(peakWidth);
    double sum2 = 0;
    for (unsigned j = baseLength + peakWidth; j < rawBr[ib]->rdigi.size(); ++j)
      sum2 += (double)rawBr[ib]->rdigi[j] - base;
    // sum2 /= double(rawBr[ib]->rdigi.size() - baseLength - peakWidth);

    double sum = 0;
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
      sum += (double)rawBr[ib]->rdigi[j] - base;
    // sum /= double(rawBr[ib]->rdigi.size());

    TString hname;
    hEvRawWave = NULL;

    // baseline correction from fitted Gaussian
    hEvGaus[ib]->Reset("ICES");
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
      hEvGaus[ib]->Fill((double)rawBr[ib]->rdigi[j] - base);
    hEvGaus[ib]->Fit("gaus", "QO");
    TF1 *gfit = (TF1 *)hEvGaus[ib]->GetListOfFunctions()->FindObject("gaus");
    double ave = hEvGaus[ib]->GetMean();
    double sigma = hEvGaus[ib]->GetRMS();
    double skew = 0;
    if (!isnan(hEvGaus[ib]->GetSkewness()))
      skew = hEvGaus[ib]->GetSkewness();
    if (gfit != nullptr)
    {
      ave = gfit->GetParameter(1);
      sigma = gfit->GetParameter(2);
    }

    // fitptr->Print();
    noiseHist[ib]->Fill(sigma);
    skewHist[ib]->Fill(skew);
    double sign = TMath::Sign(1., skew);
    TDet *idet = tbrun->getDet(ichan);
    if (idet == NULL)
    {
      printf("!!!!!NULL idet br %u ichan %i\n", ib, ichan);
      continue;
    }
    idet->ave = ave;
    idet->sigma = sigma;
    idet->skew = skew;
    idet->event = entry;
    idet->base = base;
    idet->peak = peak;
    idet->sum2 = sum2;
    idet->sum = sum;

    // cout << "tree " << ib << " ch " << ichan << " "
    //     << " ave " << ave << " sigma " << sigma << " skew " << skew << " base " << base << endl;

    /* fill baseline subtracted vector */
    digi.clear();
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = sign * ((double)rawBr[ib]->rdigi[j] - base);
      digi.push_back(val);
    }
    double thresholdValue = 10.0;
    Double_t cut = idet->sigma * thresholdValue;
    thresholdCount(idet, thresholdValue);
    unsigned diffStep = 3;
    differentiate(diffStep);     // fill ddigi
    derivativeCount(idet, 10.0); // use ddigi
   

    idet->thresholds = int(thresholds.size());
    idet->crossings = int(crossings.size());
    threshHist[ib]->Fill(idet->thresholds);
    crossHist[ib]->Fill(idet->crossings);

    ntChan->Fill(float(rawBr[ib]->trigger), float(ichan), float(ave), float(sigma), float(skew), float(base), float(peak), float(sum2), float(sum), float(crossings.size()), float(thresholds.size()), float(cut));
    // set cut 
    idet->pass = true;
    if (idet->crossings > 50)
      idet->pass = false;
    bool trig = ichan == 9 || ichan == 10 || ichan == 11;
    if (trig && idet->thresholds > 1)
      idet->pass == false;

    // plot some events
    if (!(idet->pass) && badDir->GetList()->GetEntries() < 100)
    {
      badDir->cd();
      hname.Form("EvRawWaveEv%ich%ithresh%icross%i", int(entry), ichan, idet->thresholds, idet->crossings);
      // cout << " failed " << hname << endl;
      hEvRawWave = new TH1D(hname, hname, nbins, 0, nbins);
      for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
      {
        // double val = sign * ((double)rawBr[ichan]->rdigi[j] - base - ave - valHist[ichan]->GetMean());
        double val = (double)rawBr[ib]->rdigi[j] - base;
        hEvRawWave->SetBinContent(j + 1, val);
      }
      fout->cd();
    }

    if (idet->pass && idet->thresholds > 0 && evDir->GetList()->GetEntries() < 100)
    {
      evDir->cd();
      hname.Form("EvRawWaveEv%ich%i", int(entry), ichan);
      hEvRawWave = new TH1D(hname, hname, nbins, 0, nbins);
      for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
      {
        // double val = sign * ((double)rawBr[ichan]->rdigi[j] - base - ave - valHist[ichan]->GetMean());
        double val = (double)rawBr[ib]->rdigi[j] - base;
        hEvRawWave->SetBinContent(j + 1, val);
      }
    }
    fout->cd();

  } // channel loop
  // fill output Tree for each event
  // cout << "finished " << entry << endl;
  tbrun->fill();

  bool eventPass = true;
  for (int ib = 0; ib < rawBr.size(); ++ib)
  {
    unsigned ichan = chanList[ib];
    TDet *idet = tbrun->getDet(ichan);
    if (!idet->pass)
      eventPass = false;
  }

  evCount->Fill(0);
  if (!eventPass)
    return eventPass;

  for (int ib = 0; ib < rawBr.size(); ++ib)
  {
    unsigned ichan = chanList[ib];
    int nbins = rawBr[ib]->rdigi.size();
    if (nbins != 1024)
      continue;
    TDet *idet = tbrun->getDet(ichan);
    if (idet == NULL)
    {
      printf("!!!!!NULL idet br %u ichan %i\n", ib, ichan);
      continue;
    }
    if (idet->thresholds < 1)
      continue;

    /* fill the summed histograms */

    evCount->Fill(ichan + 1);
    digi.clear();
    for (unsigned j = 0; j < rawBr[ib]->rdigi.size(); ++j)
    {
      double val = (double)rawBr[ib]->rdigi[j] - idet->base;
      digi.push_back(val);
      sumWave[ib]->SetBinContent(j + 1, sumWave[ib]->GetBinContent(j + 1) + val);
      valHist[ib]->Fill(val);
    }
    /* pulse finding  */
    finder->event(ichan, entry, digi, 3. ,3);
    finder->plotEvent(ichan,entry);
  } // second channel loop
  return eventPass;
} // anaEvent

// revised derivative Dec 8 2022 MG
void anaRun::differentiate(unsigned diffStep)
{
  ddigi.clear();
  Double_t sump = 0;
  Double_t summ = 0;
  unsigned nsamples = digi.size();
  ddigi.push_back(0); // first entry is zero
  for (unsigned i = 1; i < nsamples; ++i)
  {
    // sum limit
    int maxSum = diffStep;
    if (i < diffStep)
      maxSum = i;
    if (nsamples - 1 - i < diffStep)
      maxSum = nsamples - 1 - i;
    //
    sump = 0;
    for (unsigned j = 1; j < maxSum; ++j)
      sump += digi[i + j];
    summ = 0;
    for (unsigned j = 1; j < maxSum; ++j)
      sump -= digi[i - j];
    ddigi.push_back(sump - summ);
  }
}

void anaRun::thresholdCount(TDet *idet, Double_t rms)
{
  thresholds.clear();
  Double_t cut = idet->sigma * rms;
  unsigned step = 10;
  unsigned ibin = 0;
  while (ibin < digi.size() - step)
  {
    Double_t vi = digi[ibin];
    Double_t vj = digi[ibin + step];
    if (vi > 0 && vi < cut && vj > cut)
    {
      thresholds.push_back(UPCROSS);
    }
    ibin += step;
  }
}

void anaRun::derivativeCount(TDet *idet, Double_t rms)
{
  crossings.clear();
  crossingBin.clear();
  crossingTime.clear();
  unsigned vsize = ddigi.size();
  double microSec = 1.0E-3;
  double timeUnit = 1.0;
  Double_t cut = idet->sigma * rms;
  unsigned step = 1;
  // cout << " for det " << idet->channel  << " in derivative peaks >>>> rms " << rms << " cut " << cut << endl;
  Double_t ncut = -cut;
  // find all crossings
  for (unsigned ibin = step; ibin < vsize; ++ibin)
  {
    Double_t u = double(ibin) * timeUnit * microSec;
    Double_t vi = ddigi[ibin];
    Double_t vj = ddigi[ibin - step];
    unsigned ctype = 10;
    // if (idet == 5)
    // printf(" det %i  bin %i %f %f  \n", idet, ibin, vj, vi);
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
  }
  return;
}

anaRun::anaRun(const char *theTag, Long64_t maxEntries)
{
  tag = TString(theTag);
  cout << " starting anaRun entries = " << maxEntries << " tag =  " << tag << endl;
  if (!openFile())
  {
    printf("cannot open file!\n");
    return;
  }
  getTrees();
  cout << "treeList size " << treeList.size() << endl;
  if (treeList.size() < 1)
  {
    printf("EXIT no trees\n");
    return;
  }
  rawBr.clear();
  rawBr.resize(treeList.size());
  getEvent(0);

  fout = new TFile(Form("myData/anaRun-%s-%llu.root", tag.Data(), maxEntries), "recreate");
  evDir = fout->mkdir("evDir");
  badDir = fout->mkdir("badDir");
  fout->cd();
  cout << " opened output file " << fout->GetName() << endl;
  tbrun = new TBRun(tag);
  for (unsigned it = 0; it < chanList.size(); ++it)
    tbrun->addDet(chanList[it]);

  // fout->Append(tbrun->btree);
  ntChan = new TNtuple("ntChan", "channel ntuple", "trig:chan:ave:sigma:skew:base:peak:sum2:sum:crossings:thresholds:cut");
  evCount = new TH1D("eventCount", "event count", 14, 0, 14);

  anaDir = fout->mkdir("anaDir");
  anaDir->cd();

  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = chanList[i];
    noiseHist.push_back(new TH1D(Form("noiseChan%i", ichan), Form("noiseChan%i", ichan), 1000, 0, 1000));
    skewHist.push_back(new TH1D(Form("skewChan%i", ichan), Form("skewChan%i", ichan), 200, -3, 7));
    threshHist.push_back(new TH1D(Form("threshChan%i", ichan), Form("threshChan%i", ichan), 20, 0, 20));
    crossHist.push_back(new TH1D(Form("crossChan%i", ichan), Form("crossChan%i", ichan), 100, 0, 100));
    if (ichan > 8 && ichan < 12)
    {
      valHist.push_back(new TH1D(Form("valChan%i", ichan), Form("valChan%i", ichan), 1500, -500, 1000));
      hEvGaus.push_back(new TH1D(Form("evGaus%i", ichan), Form("evGaus%i", ichan), 200, -500, 500));
    }
    else
    {
      valHist.push_back(new TH1D(Form("valChan%i", ichan), Form("valChan%i", ichan), 1000, -200, 200));
      hEvGaus.push_back(new TH1D(Form("evGaus%i", ichan), Form("evGaus%i", ichan), 200, -100, 100));
    }
  }
  sumDir = fout->mkdir("sumDir");
  sumDir->cd();
  for (unsigned i = 0; i < rawBr.size(); ++i)
  {
    unsigned ichan = chanList[i];
    sumWave.push_back(new TH1D(Form("sumWave%i", ichan), Form("sumWave%i", ichan), rawBr[0]->rdigi.size(), 0, rawBr[0]->rdigi.size()));
  }
  fout->cd();
  // fout->ls();

  finder = new hitFinder(fout, tbrun, tag, rawBr[0]->rdigi.size(), chanList);
  Long64_t nentries = treeList[0]->GetEntries();
  if (maxEntries > 0)
    nentries = TMath::Min(maxEntries, nentries);
  printf("... total entries  %llu looping over %llu \n ", treeList[0]->GetEntries(), nentries);
  unsigned npass = 0;
  unsigned nfail = 0;
  for (Long64_t entry = 0; entry < nentries; ++entry)
  {
    if (entry / 1000 * 1000 == entry)
      printf("... %llu pass %u fail %u \n", entry, npass, nfail);
    getEvent(entry);
    /* test
    for (unsigned it = 0; it < treeList.size(); ++it)
      printf("tree %u chan %i event %i time %lld \n", it, rawBr[it]->channel, rawBr[it]->trigger, rawBr[it]->time);
    */
    if (anaEvent(entry))
      ++npass;
    else
      ++nfail;
  }

  // finder->hPeakCount->Print("all");
  printf(" npass %u nfail %u tbrun entries %llu fout %s \n", npass, nfail, tbrun->btree->GetEntriesFast(), fout->GetName());
  // tbrun->btree->GetListOfBranches()->ls();
  // fout->ls();
  fout->Write();
  fout->Close();
}
