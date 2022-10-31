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
#include "hitFinder.hxx"

class anaRun
{
public:
  enum
  {
    NDET = 4
  };

  int detId[NDET];
  ofstream dumpFile;
  TFile *fout;
  TString runTag;
  TBRun *tbrun;
  vector<TBWave *> waveList;
  vector<double> rdigi;
  hitFinder *finder;
  TString tag;
  anaRun(Long64_t maxEntries = 0, const char *theTag = "run");
  virtual ~anaRun() { ; }
  Long64_t nentries;
};

anaRun::anaRun(Long64_t maxEntries, const char *theTag)
{
  tag = TString(theTag);
  fout = new TFile(Form("anaRun-%s.root", tag.Data()), "recreate");

  TString dirName;
  dirName.Form("data/XenonDoped1_5PPMSiPM4/DAQ/%s/RAW", tag.Data());
  // dirName.Form("data/SiPM_1_inside/DAQ/%s/RAW",tag.Data());
  cout << tag << "  dirName " << dirName << endl;
  TSystemDirectory dir("RAW", dirName);
  TList *files = dir.GetListOfFiles();
  tbrun = new TBRun(tag);
  cout << " create TBRun  " << tbrun->GetName() << endl;

  TIter next(files);
  TSystemFile *file;
  // get files
  while ((file = (TSystemFile *)next()))
  {
    string name = string(file->GetName());
    cout << name << endl;
    string exten = name.substr(name.find_last_of(".") + 1);
    if (exten != string("root"))
      continue;
    string tag = name.substr(0, name.find_last_of("."));
    string tdet = name.substr(name.find_last_of(".")-1,1);
    // string chan = tag.substr( tag.find_first_of("CH"), tag.find_first_of("@") - tag.find_first_of("CH"));
    string chan = tag.substr(tag.find_first_of("det"), tag.find_first_of(".root") - tag.find_first_of("det"));
    // cout << " tag " << tag << "  " << tag.find_first_of("CH") <<  "   " << " chan  " << chan << endl;
    string fullName = string(dirName.Data()) + string("/") + name;
    cout << " tag " << tag << " open " << fullName << " tdet " << tdet << endl;
    TFile *f = new TFile(fullName.c_str(), "readonly");
    if (!f)
      continue;
    TTree *dtree = NULL;
    f->GetObject("Data_R", dtree);
    if (!dtree)
      continue;
    TBranch *b_Samples = dtree->FindBranch("Samples");
    if (!b_Samples)
      continue;
    if (b_Samples->GetEntries() < 1)
    {
      cout << " skipping  " << tag << endl;
      continue;
    }
    cout << f->GetName() << "  " << b_Samples->GetEntries() << endl;
    TBWave *wave = new TBWave(tag.c_str());
    TTree *tree = NULL;
    f->GetObject("Data_R", tree);
    wave->Init(tree);
    waveList.push_back(wave);
    detId[waveList.size() - 1] = stoi(tdet);
    cout << " file " << wave->GetName() << " open " << waveList.size() << " detId " << detId[waveList.size() - 1] << endl;
    tbrun->addDet(Form("DET%i", detId[waveList.size() - 1]), wave->GetName());
  }

 
  cout << " got " << waveList.size() << " files " << endl;

  if (waveList.size() < 1)
  {
    fout->Close();
    return;
  }

  for (unsigned ifile = 0; ifile < waveList.size(); ++ifile)
  {
    Long64_t nentries = waveList[ifile]->fChain->GetEntries();
    cout << ifile << " for wave  " << waveList[ifile]->GetName() << " det entries  " << nentries << endl;
  }

// loop over files  
  finder = NULL;
  for (unsigned ifile = 0; ifile < waveList.size(); ++ifile)
  {
    cout << " run finder " << ifile << " for wave  " << waveList[ifile]->GetName() << " det " << tbrun->detList[ifile]->GetName() << endl;
    Long64_t nbytes = 0, nb = 0;
    Long64_t ientry = waveList[ifile]->LoadTree(0);
    nb = waveList[ifile]->fChain->GetEntry(0);
    int nsamples = waveList[ifile]->Samples->GetSize();
    if(!finder) finder = new hitFinder(fout, tbrun, tag,nsamples);

    fout->cd();

    Long64_t nentries = waveList[ifile]->fChain->GetEntries();
    if (maxEntries > 0)
      nentries = maxEntries;
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
      if (jentry / 1000 * 1000 == jentry)
        printf("\t ... det %u event %lld \n", ifile , jentry);
      ientry = waveList[ifile]->LoadTree(jentry);
      if (ientry < 0)
        break;
      nb = waveList[ifile]->GetEntry(jentry);
      nbytes += nb;
      rdigi.clear();
      for (int i = 0; i < nsamples; ++i)
        rdigi.push_back(waveList[ifile]->Samples->GetAt(i));
      //
      finder->event(ifile, jentry, rdigi);

      if (1) { //finder->detHits.size() >0 ) {
        finder->plotWave(ifile, jentry);
      }
      finder->plotEvent(ifile, jentry);
    }
  }
  for (int i = 0; i < tbrun->detList.size(); ++i)
    printf(" %i detId %i %s \n", i, detId[i],tbrun->detList[i]->GetName());
  // fout->ls();
  //
  
  printf("... ending ");
  
  finder->hPeakCount->Print("all");
  //fout->ls();
  fout->Write();
  fout->Close();
}
