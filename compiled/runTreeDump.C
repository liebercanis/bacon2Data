#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TBranchElement.h"
#include <vector>
#include "TDet.hxx"

std::vector<TDet *> detList;
std::vector<long> hocc;
TFile *fin = NULL;
TTree *RunTree = NULL;

bool openFile(TString fileName)
{
  // open input file and make some histograms
  printf(" looking for file %s\n", fileName.Data());

  bool exists = false;
  FILE *aFile;
  aFile = fopen(fileName.Data(), "r");
  if (aFile)
  {
    fclose(aFile);
    exists = true;
  }
  if (!exists)
  {
    printf(" couldnt open file %s\n", fileName.Data());
    return false;
  }

  exists = true;
  fin = new TFile(fileName, "readonly");
  printf(" opened file %s\n", fileName.Data());
  RunTree = NULL;
  fin->GetObject("RunTree", RunTree);
  if (!RunTree)
  {
    printf(" no TTree* RunTree in file %s\n", fileName.Data());
    exists = false;
  }
  return exists;
}

/* get rawBr */
unsigned long getBranches()
{
  detList.clear();
  TObjArray *brList = RunTree->GetListOfBranches();
  TIter next(brList);
  TBranchElement *aBranch = NULL;
  while ((aBranch = (TBranchElement *)next()))
  {
    // aBranch->SetAddress(0);
    if (TString(aBranch->GetName()).Contains("eventData"))
      continue;
    int ichan = TString(TString(aBranch->GetName())(4, 2)).Atoi();
    // cout << "branch " << aBranch->GetName() << "..." << aBranch->GetClass()->GetName() << " idet " << ichan << endl;
    if (ichan > 12)
      break;
    TDet *det = (TDet *)aBranch->GetObject();
    det->SetName(aBranch->GetName());
    // printf(" %ld det %s \n", detList.size(), det->GetName());
    detList.push_back(det);
  }
  return detList.size();
}

void dump(int ichan = 8, Long64_t nev = 0)
{
  if (nev == 0)
    nev = RunTree->GetEntries();
  else
    nev = min(nev, RunTree->GetEntries());

  // loop over events
  for (Long64_t entry = 0; entry < nev; ++entry)
  {

    RunTree->GetEntry(entry);
    getBranches();

    // if first time, init vector
    if (hocc.size() == 0)
      hocc.resize(detList.size(), 0);

    // sum hits for all channel
    for (unsigned id = 0; id < detList.size(); ++id)
    {
      hocc[id] += detList[id]->hits.size();
      // printf("event %lld chan %u %s hits %lu tot %lu \n", entry, id, detList[id]->GetName(), detList[id]->hits.size(), hocc[id]);
    }

    TDet *det = detList[ichan + 1];
    if (det->hits.size() > 0)
      printf("dump event %lld chan %u hits %ld \n", entry, ichan, det->hits.size());

    for (unsigned ih = 0; ih < det->hits.size(); ++ih)
    {
      TDetHit thit = det->hits[ih];
      printf("\t ih %d start %f qpeak %f \n", ih, thit.startTime, thit.qpeak);
    }
  }
}

// caenData soft link to appropriate directory
void runTreeDump(TString fileName = "caenData/anaCRun-run-09_01_2024-file-0.root")
{
  if (!openFile(fileName))
  {
    return;
  }
  printf("successfully opened %s with RunTree entries %lld \n", fin->GetName(), RunTree->GetEntries());
  int ichan = 8;
  Long64_t nev = RunTree->GetEntries();
  dump(ichan, nev);
  printf("hit summary %lld \n", nev);
  for (unsigned id = 0; id < detList.size(); ++id)
  {
    printf("chan %i occ %ld \n", id, hocc[id]);
  }
}
