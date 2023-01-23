/*  read raw data from SIS3316 125MHz MG, Dec 15 2022 */
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TString.h>
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TBranchElement.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TKey.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TString.h"
#include "TObjString.h"
#include "TSystemDirectory.h"
#include "TFile.h"
#include "TBRawRun.hxx"

using namespace std;

int maxEvents;
TFile *fin;
TFile *fout;

vector<TTree *> treeList;
vector<int> chanNum;
vector<TBRawEvent *> rawBr;
vector<TBRawEvent *> rawEv;

void getTrees() {
  //
  TIter next(fin->GetListOfKeys());
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TTree"))
      continue;
    TTree *tree  = (TTree *)key->ReadObj();
    //cout << tree->GetName() << " " << key->GetCycle() << endl;
    if (tree->GetEntries() < 1)
      continue;
    // pick up only last in cycle 
    int ichan = TString(TString(tree->GetName())(5, 2)).Atoi();
    if( find(chanNum.begin(),chanNum.end(),ichan ) !=chanNum.end())
      continue;
    treeList.push_back(tree);
    chanNum.push_back(ichan);
  }
}

void getBranches()
{
  for (unsigned it = 0; it < treeList.size(); ++it)
  {
    TObjArray *brList  = treeList[it]->GetListOfBranches();
    TIter next(brList);
    TBranchElement *aBranch = NULL;
    while ((aBranch = (TBranchElement *)next()))
    {
      TBRawEvent *chan = (TBRawEvent *)aBranch->GetObject();
      rawBr[it] = chan;
      ;
    }
  }
}

void getEvent(Long64_t entry) {
  for (unsigned it = 0; it < treeList.size(); ++it)
  {
    treeList[it]->GetEntry(entry);
    getBranches();
  }
}

int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << " argc " << argc << endl;
  if (argc < 2)
  {
    printf(" usage: build <tag>  <number  of events>  \n ");
    exit(0);
  }

  TString tag(argv[1]);
  maxEvents = 0;
  if (argc > 2)
    maxEvents = atoi(argv[2]);
  //
  TString fileName = TString("data/rootData/") + tag + TString(".root");
  fin = new TFile(fileName, "readonly");
  TString foutName = TString("data/rootData/build-") + tag + TString(".root");
  fout = new TFile(foutName, "recreate");


  getTrees();
  for (unsigned it = 0; it < treeList.size(); ++it)
    printf("%i tree %s entries %lld chan %i \n",it,treeList[it]->GetName(), treeList[it]->GetEntries(),chanNum[it]);

  

  rawBr.resize(treeList.size());
  for (Long64_t entry = 0; entry < 10; ++entry)
  {
    getEvent(entry);
    // look at data
    for (unsigned it = 0; it < treeList.size(); ++it)
    {
      printf("tree %u chan %i event %i time %lld \n", it, rawBr[it]->channel, rawBr[it]->trigger, rawBr[it]->time);
    }
  }

  fin->Close();
  fout->Write();
  fout->Close();
  return 0;
}
