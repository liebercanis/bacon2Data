#include <TDet.hxx>
#include <TDetHit.hxx>
#include <TBranchElement.h>
TFile *fin;
std::vector<unsigned> vchan;
std::vector<TDet *> detList;
TTree *RunTree;
std::vector<TH1D *> hMult;
std::vector<TH1D *> hQPeak;
std::vector<TH1D *> hQSum;
vector<TH1D *> sumHitWave;
std::map<int, int> chanMap;
TH1D *hTimeDiff;

bool exists(TString name)
{
  ifstream f(name.Data());
  return f.good();
  f.close();
}

/* get rawBr */
unsigned getBranches()
{
  detList.clear();
  TObjArray *brList = RunTree->GetListOfBranches();
  TIter next(brList);
  TBranchElement *aBranch = NULL;
  while ((aBranch = (TBranchElement *)next()))
  {
    if (TString(aBranch->GetName()).Contains("event"))
      continue;
    int ichan = TString(TString(aBranch->GetName())(4, 2)).Atoi();
    cout << "branch " << aBranch->GetName() << "..." << aBranch->GetClass()->GetName() << " idet " << ichan << endl;
    vchan.push_back(ichan);
    aBranch->SetAddress(0);
    chanMap.insert(std::pair<int, int>(ichan, vchan.size() - 1));
  }
  return vchan.size();
}

void anaEvent(Long64_t entry)
{
  RunTree->GetEntry(entry);
  // RunTree->GetListOfBranches()->ls();
  //  find branches and cast as TDet each is a channel
  TIter next(RunTree->GetListOfBranches());
  TBranchElement *aBranch = NULL;
  detList.clear();
  while ((aBranch = (TBranchElement *)next()))
  {
    TDet *det = (TDet *)aBranch->GetObject();
    if (TString(det->GetName()).Contains("event"))
      continue;
    detList.push_back(det);
  }

  if (entry == 0)
    for (unsigned id = 0; id < detList.size(); ++id)
      printf("\t\t %u %s \n", id, detList[id]->GetName());

  double nominalGain = 227.4;

  for (unsigned id = 0; id < detList.size(); ++id)
  {
    // loop over hits
    double cut = 100;
    if (id == 12)
      cut = 200.;
    int nhits = 0;
    // printf("chan %u %s hits %lu \n", id, detList[id]->GetName(), detList[id]->hits.size());
    for (unsigned ih = 0; ih < detList[id]->hits.size(); ++ih)
    {
      hQPeak[id]->Fill(detList[id]->hits[ih].qpeak);
      if (detList[id]->hits[ih].qpeak > cut)
        ++nhits;
    }
    // fill mult for this event
    hMult[id]->Fill(double(nhits));
  }
}

void hits(TString fileName = TString("anaCRun-run-09_10_2024-file_9-0.root"), Long64_t maxEntries = 0)
{
  int nsamples = 1024;
  gStyle->SetOptStat(1001101);
  TString fullName = TString("caenData/") + fileName;
  cout << " looking for " << fullName << endl;
  TString tag("unknown"); // needs to be parsed fileName !
  if (!exists(fullName))
  {
    cout << " didnt find " << fullName << endl;
    return;
  }
  fin = new TFile(fullName, "readonly");

  cout << " opened  " << fileName << endl;
  RunTree = NULL;
  fin->GetObject("RunTree", RunTree);
  if (!RunTree)
  {
    cout << " no RunTree " << endl;
    return;
  }

  Long64_t ntriggers = RunTree->GetEntries();
  cout << " ntriggers " << ntriggers << endl;

  TString foutName = TString("hits-") + tag + TString(".root");
  TFile *fout = new TFile(foutName, "recreate");

  hTimeDiff = new TH1D("TimeDiff", "hit time diff", 1000, 0, 1000);

  unsigned ndets = getBranches();
  cout << " number of dets " << ndets << endl;

  /*for (unsigned index = 0; index < vchan.size(); ++index)
    {
    int id = vchan[index];
    chanMap.insert(std::pair<int, int>(id, index));
    }
    */

  // printf(" channel mapping \n");
  for (unsigned index = 0; index < vchan.size(); ++index)
  {
    int id = chanMap.at(vchan[index]);
    printf("index %i chan %i mapped to index  %i \n", index, vchan[index], id);
  }

  double limit;
  fout->cd();
  hQPeak.clear();
  for (unsigned ichan = 0; ichan < vchan.size(); ++ichan)
  {
    // cout << " hist " << ichan << endl;
    int limit = 1000;
    hQPeak.push_back(new TH1D(Form("QPeakChan%i", ichan), Form("QPeakChan%i", ichan), 1000, 0, limit));
  }
  hMult.clear();
  for (unsigned ichan = 0; ichan < vchan.size(); ++ichan)
  {
    hMult.push_back(new TH1D(Form("HitMultChan%i", ichan), Form("HitMultChan%i", ichan), 10, 0, 10));
  }

  for (Long64_t entry = 0; entry < ntriggers; ++entry)
    anaEvent(entry);

  // stop here
  fout->Write();
  hMult[12]->Print("all");
  return;

  for (int id = 0; id < sumHitWave.size(); ++id)
  {
    printf("id %i entries %f\n", id, sumHitWave[id]->GetEntries());
    sumHitWave[id]->SetMarkerStyle(21);
    sumHitWave[id]->SetMarkerSize(0.4);

    if (sumHitWave[id]->GetEntries() > 0)
      sumHitWave[id]->Scale(1. / sumHitWave[id]->GetEntries());
  }

  TCanvas *can = new TCanvas("qsum", "qsum");
  can->Divide(3, 3);
  unsigned ican = 0;
  vector<double> qpe;
  vector<double> eqpe;
  qpe.resize(vchan.size());
  eqpe.resize(vchan.size());
  for (unsigned j = 0; j < vchan.size(); ++j)
  {
    unsigned id = chanMap.at(vchan[j]);
    qpe[id] = -10E9;
    eqpe[id] = -10E9;
    if (vchan[j] < 9)
    {
      ++ican;
      can->cd(ican);
      hQSum[id]->Fit("gaus", "", "", 1500, 4000);
      TF1 *g = (TF1 *)hQSum[id]->GetListOfFunctions()->FindObject("gaus");

      if (g)
      {
        hQSum[id]->Draw();
        qpe[id] = g->GetParameter(1);
        eqpe[id] = g->GetParError(1);
      }
    }
  }

  TCanvas *canw = new TCanvas("wavesum", "wavesum");
  canw->Divide(3, 3);
  for (unsigned j = 3; j < 12; ++j)
  {
    canw->cd(j - 3 + 1);
    gPad->SetLogy();
    sumHitWave[j]->Draw();
  }

  // fout->ls();
  for (unsigned j = 0; j < vchan.size(); ++j)
    if (qpe[j] > -9E9)
      printf(" chan %i mean %f %f \n", vchan[j], qpe[j], eqpe[j]);

  fout->Write();
}
