#include <TDet.hxx>
TFile *fin;
std::vector<unsigned> vchan;
std::vector<TDet*> detList;
TTree *RunTree;
std::vector<TH1D*> hQSum;
std::vector<TH1D*> hQPeak;
vector<TH1D *> sumHitWave;
std::map<int, int> chanMap;

bool exists(TString name) {
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
    int ichan = TString(TString(aBranch->GetName())(4, 2)).Atoi();
    cout << "branch " << aBranch->GetName() << "..." << aBranch->GetClass()->GetName() << " idet " << ichan << endl;
    vchan.push_back(ichan); 
    aBranch->SetAddress(0);
    chanMap.insert(std::pair<int, int>(ichan, vchan.size()-1));
  }
  return vchan.size();
}


void anaEvent(Long64_t entry) {
  RunTree->GetEntry(entry);
  // find branches and cast as TDet each is a channel
  TIter next(RunTree->GetListOfBranches());
  TBranchElement *aBranch = NULL;
  detList.clear();
  while ((aBranch = (TBranchElement *) next())) {
    TDet *det = (TDet *) aBranch->GetObject();
    //if(entry==0) printf(" det %s \n", det->GetName());
    detList.push_back(det);
  }

  for(unsigned id=0; id<detList.size(); ++id) {
    for(unsigned ih=0; ih<detList[id]->hits.size(); ++ih) {
      if (detList[id]->hits[ih].qsum < 1000) // low threshold cut
        continue;
      hQSum[id]->Fill(detList[id]->hits[ih].qsum);
      hQPeak[id]->Fill(detList[id]->hits[ih].qpeak);
      TDetHit thit = detList[id]->hits[ih];
      sumHitWave[id]->SetBinContent(thit.peakBin + 1, sumHitWave[id]->GetBinContent(thit.peakBin + 1) + thit.qsum);
    }
  }
}

void hits(TString fileName = TString("anaRun-run-03_01_2023-file14-0.root"), Long64_t maxEntries = 0)
{
  int nsamples = 1024;
  gStyle->SetOptStat(1001101);
  TString fullName = TString("myData/") + fileName;
  TString tag("unknown"); // needs to be parsed fileName !
  if(!exists(fullName)) {
    cout << " didnt find " << fullName << endl;
    return;
  }
  fin = new TFile(fullName, "readonly");

  cout << " opened  " << fileName << endl;
  RunTree=NULL;
  fin->GetObject("RunTree",RunTree);
  if(!RunTree) {
    cout << " no RunTree " << endl;
    return;
  }

  Long64_t ntriggers = RunTree->GetEntries();
  cout << " ntriggers " << ntriggers << endl;

  TString foutName = TString("hits-") + tag + TString(".root");
  TFile *fout = new TFile(foutName, "recreate");

  unsigned ndets = getBranches();
  cout << " number of dets " << ndets << endl;

  /*for (unsigned index = 0; index < vchan.size(); ++index)
    {
    int id = vchan[index];
    chanMap.insert(std::pair<int, int>(id, index));
    }
    */

  printf(" channel mapping \n");
  for (unsigned index = 0; index < vchan.size(); ++index)
  {
    int id = chanMap.at(vchan[index]);
    printf("index %i chan %i mapped to index  %i \n", index, vchan[index], id);
  }

  double limit;
  for (unsigned ichan = 0; ichan < vchan.size(); ++ichan)
  {
    bool trigger = vchan[ichan] == 9 || vchan[ichan] == 10 || vchan[ichan] == 11;
    if (trigger)
      limit = 200000;
    // else if(vchan[ichan]==12) limit = 10000;
    else
      limit = 20000;

    hQSum.push_back(new TH1D(Form("QSumChan%i", vchan[ichan]), Form("QSumChan%i", vchan[ichan]), 1000, 0, limit));
    hQPeak.push_back(new TH1D(Form("QPeakChan%i", vchan[ichan]), Form("QPeakChan%i", vchan[ichan]), 1000, 0, limit));
    sumHitWave.push_back(new TH1D(Form("sumHitWave%i", vchan[ichan]), Form("sumHitWave%i", vchan[ichan]), nsamples, 0, nsamples));
  }

  for (Long64_t entry = 0; entry < ntriggers; ++entry)
    anaEvent(entry);

  TCanvas *can = new TCanvas("qsum", "qsum");
  can->Divide(3, 3);
  unsigned  ican = 0;
  vector<double> qpe;
  vector<double> eqpe;
  qpe.resize(vchan.size());
  eqpe.resize(vchan.size());
  for (unsigned j = 0; j < vchan.size(); ++j)
  {
    unsigned id = chanMap.at(vchan[j]);
    qpe[id]=-10E9;
    eqpe[id]=-10E9;
    if(vchan[j]<9) {
      ++ican;
      can->cd(ican);
      hQSum[id]->Fit("gaus", "", "", 1500, 4000);
      TF1 *g = (TF1 *)hQSum[id]->GetListOfFunctions()->FindObject("gaus");
      
      hQSum[id]->Draw();
      qpe[id] = g->GetParameter(1);
      eqpe[id] = g->GetParError(1);
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
    if(qpe[j]>-9E9) printf(" chan %i mean %f %f \n",vchan[j],qpe[j],eqpe[j]);

  fout->Write();
}
