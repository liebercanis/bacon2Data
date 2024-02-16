using namespace TMath;
TChain *RunTree;
TFile *fout;
TString tag;
TNtuple *ntSum;
TNtuple *ntHit;
enum
{
  NCHAN = 14
};
// hSipmLatePeaksFitograms as vectors
std::vector<TH1D *> hTotSum;
std::vector<TH1D *> hPreSum;
std::vector<TH1D *> hTrigSum;
std::vector<TH1D *> hLateSum;

double y[NCHAN];
double nominalGain;

TH1D *hTrigEventSumArea;
TH1D *hTrigEventSumPeak;
TH2D *hTrigHitPeakTime;
TH2D *hTrigHitPeakTimeCoarse;
TH2D *hAllHitPeakTimeCoarse;
TH2D *hAllHitPeakTime;

TH1D *hTrigLatePeaks;
TH1D *hSipmLatePeaks;
TH1D *hTrigLatePeaksFit;
TH1D *hSipmLatePeaksFit;
unsigned sipmCut;
unsigned trigCut;

void loop(Long64_t maxEntry)
{
  sipmCut = 0;
  trigCut = 0;
  // loop over entries
  for (Long64_t entry = 0; entry < maxEntry; ++entry)
  {
    // get entry
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
    double trigEventSumArea = 0;
    double trigEventSumPeak = 0;
    double trigTotPeak = 0;

    bool trigIsCut = false;
    bool sipmIsCut = false;

    while ((aBranch = (TBranchElement *)next()))
    {
      int id = TString(TString(aBranch->GetName())(4, 2)).Atoi();
      if (id >= NCHAN) // skip PMT
        continue;
      bool trig = false; // define trigger sipms
      if (id == 9 || id == 10 || id == 11)
        trig = true;
      if (TString(aBranch->GetName()) == TString("eventData")) // skip this branch
        continue;
      TDet *det = (TDet *)aBranch->GetObject();
      // if(id==0) cout << "branch " << aBranch->GetName() << " entry " << entry << " TotSum " << det->totSum << endl;
      hTotSum[id]->Fill(det->totPeakSum / nominalGain);
      hPreSum[id]->Fill(det->prePeakSum / nominalGain);
      hTrigSum[id]->Fill(det->trigPeakSum / nominalGain);
      hLateSum[id]->Fill(det->latePeakSum / nominalGain);
      if (trig)
      {
        trigPre += det->prePeakSum / nominalGain;
        trigTrig += det->trigPeakSum / nominalGain;
        trigLate += det->latePeakSum / nominalGain;
        trigEventSumPeak += det->trigPeakSum / nominalGain;
        trigEventSumArea += det->trigSum;
        trigTotPeak += (det->trigPeakSum + det->latePeakSum) / nominalGain;
      }
      else
      {
        sipmPre += det->prePeakSum / nominalGain;
        sipmTrig += det->trigPeakSum / nominalGain;
        sipmLate += det->latePeakSum / nominalGain;
      }
      // loop over hits
      for (unsigned ihit = 0; ihit < det->hits.size(); ++ihit)
      {
        TDetHit hiti = det->hits[ihit];
        // printf(" det %i time %.0f qpeak %f \n",id,det->hits[ihit].startTime, det->hits[ihit].qpeak);
        if (trig)
          hTrigHitPeakTime->Fill(det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);
        if (trig)
          hTrigHitPeakTimeCoarse->Fill(det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);
        hAllHitPeakTime->Fill(det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);
        hAllHitPeakTimeCoarse->Fill(det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);
        ntHit->Fill(double(entry), double(id), det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);

        if (trig && det->hits[ihit].startTime > 800)
          hTrigLatePeaks->Fill(det->hits[ihit].qpeak / nominalGain);
        if (!trig && det->hits[ihit].startTime > 800)
          hSipmLatePeaks->Fill(det->hits[ihit].qpeak / nominalGain);

        if (trig && (det->hits[ihit].startTime > 800 || det->hits[ihit].startTime < 600) && det->hits[ihit].qpeak / nominalGain > 10)
          trigIsCut = true;
        if (!trig && (det->hits[ihit].startTime > 800 || det->hits[ihit].startTime < 600) && det->hits[ihit].qpeak / nominalGain > 10)
          sipmIsCut = true;

        if ((det->hits[ihit].startTime > 800 || det->hits[ihit].startTime < 600) && det->hits[ihit].qpeak / nominalGain > 10)
          printf("xxxxx event %llu chan %i time %f qpeak %f\n", entry, id, det->hits[ihit].startTime, det->hits[ihit].qpeak / nominalGain);
      }
    }
    if (trigEventSumPeak < 0)
      trigEventSumPeak = 0;
    ntSum->Fill(trigEventSumArea, trigEventSumPeak, trigTotPeak, trigPre, trigTrig, trigLate, sipmPre, sipmTrig, sipmLate);
    if (entry == entry / 10000 * 10000)
      printf("... event %llu trig sum %.2E %.2E \n", entry, trigEventSumArea, trigEventSumPeak);
    hTrigEventSumArea->Fill(trigEventSumArea);
    hTrigEventSumPeak->Fill(trigEventSumPeak);
    if (trigIsCut)
      ++trigCut;
    if (sipmIsCut)
      ++sipmCut;
  }
}

void post(TString tag = TString("12_28_2023"), Long64_t maxEntry = 0)
{
  /*gains-2024-02-01-17-06.root*/
  y[0] = 229.5;
  y[1] = 221;
  y[2] = 236;
  y[3] = 200.6;
  y[4] = 229;
  y[5] = 229; // not fit same as 4
  y[6] = 231.6;
  y[7] = 237.2;
  y[8] = 232.6;
  y[9] = 646.9;
  y[10] = 619.5;
  y[11] = 605;
  nominalGain = 270.5;
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
  if (maxEntry == 0)
    maxEntry = ntriggers;
  RunTree->GetListOfBranches()->ls();
  TString sentries;
  sentries.Form("-%llu", maxEntry);
  fout = new TFile(TString("post-") + tag + sentries + TString(".root"), "recreate");

  // make hSipmLatePeaksFitos
  for (unsigned i = 0; i < NCHAN; ++i)
  {
    double limit = 40;
    int nbins = 400.;
    // normalized to SPE
    hTotSum.push_back(new TH1D(Form("TotPeakSumChan%i", i), Form("tot peak sum chan %i", i), nbins, 0, limit));
    hPreSum.push_back(new TH1D(Form("PrePeakSumChan%i", i), Form("pre peak sum chan %i", i), nbins, 0, limit));
    hTrigSum.push_back(new TH1D(Form("TrigPeakSumChan%i", i), Form("trig peak sum chan %i", i), nbins, 0, limit));
    hLateSum.push_back(new TH1D(Form("LatePeakSumChan%i", i), Form("late peak sum chan %i", i), nbins, 0, limit));
  }
  hTrigEventSumArea = new TH1D("TrigEventSumArea", "trig sipm trigger window sum area ", 600, 0, 6.E5);
  hTrigEventSumPeak = new TH1D("TrigEventSumPeak", "trig sipm trigger window sum peak ", 400, 0, 40);
  hTrigHitPeakTime = new TH2D("TrigHitPeakTime", "trig summed trig peak versus time ", 7500, 0, 7500, 400, 0, 40);
  hAllHitPeakTime = new TH2D("AllHitPeakTime", "trig summed trig peak versus time ", 7500, 0, 7500, 400, 0, 40);
  hTrigHitPeakTimeCoarse = new TH2D("TrigHitPeakTimeCoarse", "trig summed trig peak versus time ", 75, 0, 7500, 400, 0, 40);
  hAllHitPeakTimeCoarse = new TH2D("AllHitPeakTimeCoarse", "trig summed trig peak versus time ", 75, 0, 7500, 400, 0, 40);

  ntSum = new TNtuple("ntSum", " ADC sums ", "trigSumArea:trigSumPeak:trigTotPeak:trigPre:trigTrig:trigLate:sipmPre:sipmTrig:sipmLate");
  ntHit = new TNtuple("ntHit", " hits ", "event:chan:time:qpeak");

  hTrigLatePeaks = new TH1D("TrigEventSumPeak", "trig sipm trigger late peaks ", 1000, 0, 100);
  hTrigLatePeaksFit = new TH1D("TrigEventSumPeakFit", "trig sipm trigger late peaks ", 1000, 0, 100);
  hSipmLatePeaks = new TH1D("SipmEventSumPeak", "non-trig sipm late peaks ", 1000, 0, 100);
  hSipmLatePeaksFit = new TH1D("SipmEventSumPeakFit", "non-trig sipm late peaks ", 1000, 0, 100);

  //
  loop(maxEntry);

  // poisson fit
  TF1 *fpoi1 = new TF1("fpoi1", "[1]*pow([0],x)*Exp(-[0])/Gamma(x+1.)", 0, 10);
  double norm = hSipmLatePeaks->Integral();
  // set non-zero initial values for parameters
  fpoi1->SetParameter(0, 1);
  fpoi1->SetParameter(1, norm);
  hSipmLatePeaks->Fit("fpoi1", "R");
  for (int ib = 1; ib < hSipmLatePeaksFit->GetNbinsX(); ++ib)
  {
    double xbin = hSipmLatePeaksFit->GetBinCenter(ib) - 0.5;
    double fbin = fpoi1->Eval(xbin);
    hSipmLatePeaksFit->SetBinContent(ib, fbin);
    hSipmLatePeaksFit->SetBinError(ib, 0);
    hSipmLatePeaksFit->GetYaxis()->SetTitle("yield");
    hSipmLatePeaksFit->GetXaxis()->SetTitle("SPE");
  }
  hSipmLatePeaksFit->SetLineColor(kRed);

  TCanvas *canSipmLate = new TCanvas("sipmLate", "sipmLate");
  hSipmLatePeaks->Draw();
  hSipmLatePeaksFit->Draw("same");

  // R" = fit between "xmin" and "xmax" of the "f1"

  // poisson fit
  TF1 *fpoi2 = new TF1("fpoi2", "[1]*pow([0],x)*Exp(-[0])/Gamma(x+1.)", 0, 10);
  norm = hTrigLatePeaks->Integral();
  // set non-zero initial values for parameters
  fpoi1->SetParameter(0, 1);
  fpoi1->SetParameter(1, norm);
  hTrigLatePeaks->Fit("fpoi1", "R");
  for (int ib = 1; ib < hTrigLatePeaksFit->GetNbinsX(); ++ib)
  {
    double xbin = hTrigLatePeaksFit->GetBinCenter(ib) - 0.5;
    double fbin = fpoi1->Eval(xbin);
    hTrigLatePeaksFit->SetBinContent(ib, fbin);
    hTrigLatePeaksFit->SetBinError(ib, 0);
    hTrigLatePeaksFit->GetYaxis()->SetTitle("yield");
    hTrigLatePeaksFit->GetXaxis()->SetTitle("SPE");
  }
  hTrigLatePeaksFit->SetLineColor(kRed);

  TCanvas *canTrigLate = new TCanvas("trigLate", "trigLate");
  hTrigLatePeaks->Draw();
  hTrigLatePeaksFit->Draw("same");

  printf("total %llu trig cut %u (%f) sipm cut %u (%f) \n",
         maxEntry, trigCut, double(trigCut) / double(maxEntry), sipmCut, double(sipmCut) / double(maxEntry));

  fout->Write();
  // fout->ls();
}
