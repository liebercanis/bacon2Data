TChain *RunTree;
TFile *fout;
TString tag;
TNtuple* ntSum;
enum
{
  NCHAN = 12
};
// histograms as vectors
std::vector<TH1D *> hTotSum;
std::vector<TH1D *> hPreSum;
std::vector<TH1D *> hTrigSum;
std::vector<TH1D *> hLateSum;

void loop(Long64_t maxEntry)
{
  // loop over entries
  for (Long64_t entry = 0; entry < maxEntry; ++ entry)
  {
    // get entry 
    if(entry == entry/10000*10000)
      printf("...at event %llu\n", entry);
    RunTree->GetEntry(entry);
    // get branch pointers and save in detList
    TIter next(RunTree->GetListOfBranches());
    TBranchElement *aBranch = NULL;
    // loop over branches
    double trigPre = 0;
    double trigTrig = 0;
    double trigLate = 0;
    double sipmPre = 0;
    double sipmTrig = 0;
    double sipmLate = 0;


    while ((aBranch = (TBranchElement *)next()))
    {
      int id = TString(TString(aBranch->GetName())(4, 2)).Atoi();
      if (id >= NCHAN) // skip PMT
        continue;
      bool trig = false; // define trigger sipms
      if(id==9||id==10||id==11)
        trig = true;
      if (TString(aBranch->GetName()) == TString("eventData")) // skip this branch
        continue;
      TDet *det = (TDet *)aBranch->GetObject();
      //if(id==0) cout << "branch " << aBranch->GetName() << " entry " << entry << " TotSum " << det->totSum << endl;
      if(det->totPeakSum>0) hTotSum[id]->Fill(det->totPeakSum);
      if(det->prePeakSum > 0) hPreSum[id]->Fill(det->prePeakSum);
      if(det->trigPeakSum > 0) hTrigSum[id]->Fill(det->trigPeakSum);
      if(det->latePeakSum > 0) hLateSum[id]->Fill(det->latePeakSum);
      if(trig) {
        trigPre  +=det->prePeakSum;
        trigTrig +=det->trigPeakSum;
        trigLate +=det->latePeakSum;
      } else {
        sipmPre += det->prePeakSum;
        sipmTrig += det->trigPeakSum;
        sipmLate += det->latePeakSum;
      }
    }
    if(trigTrig>0&&sipmTrig>0) 
      ntSum->Fill(trigPre, trigTrig, trigLate, sipmPre, sipmTrig, sipmLate);
  }
}

void post(TString tag = TString("11_26_2023"), Long64_t maxEntry=0)
{
  gStyle->SetOptStat(1001101);
  cout << " the tag is  " << tag << endl;
  /* get RunTree */
  RunTree = new TChain("RunTree");
  TString name;
  name.Form("caenData/anaCRun*%s*.root", tag.Data());
  printf("open chain with %s \n", name.Data());
  RunTree->Add(name);
  if (!RunTree)
    return;
  printf("files in chain:\n");
  RunTree->GetListOfFiles()->Print();
  Long64_t ntriggers = RunTree->GetEntries();
  printf(" total triggers in this chain %lld \n", ntriggers);
  if(maxEntry==0)
    maxEntry = ntriggers;
  RunTree->GetListOfBranches()->ls();
  TString sentries;
  sentries.Form("-%llu",maxEntry);
  fout = new TFile(TString("post-") + tag + sentries + TString(".root"), "recreate");

  // make histos
  for (unsigned i = 0; i < NCHAN; ++i)
  {
    double limit = 2000;
    int nbins = 2000.;
    if (i > 8){
      limit = 5000.;
      nbins = 5000;
    }
    hTotSum.push_back(new TH1D(Form("TotPeakSumChan%i", i), Form("tot peak sum chan %i", i), nbins, 0, limit));
    hPreSum.push_back(new TH1D(Form("PrePeakSumChan%i", i), Form("pre peak sum chan %i", i), nbins, 0, limit));
    hTrigSum.push_back(new TH1D(Form("TrigPeakSumChan%i", i), Form("trig peak sum chan %i", i), nbins, 0, limit));
    hLateSum.push_back(new TH1D(Form("LatePeakSumChan%i", i), Form("late peak sum chan %i", i), nbins, 0, limit));
  }

  ntSum = new TNtuple("ntSum", " ADC sums ", "trigPre:trigTrig:trigLate:sipmPre:sipmTrig:sipmLate");

  loop(maxEntry);

  fout->Write();
  //fout->ls();
  }
