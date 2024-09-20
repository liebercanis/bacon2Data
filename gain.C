/*
   read and fit gains from summary June 24 2024
   updated sept 16 2024
*/
#include <ctime>
#include <iostream>
#include "TDirectory.h"
#include "TGraphErrors.h"
double nominalGain = 227.4; // average
std::vector<std::vector<TF1 *>> vfit;
std::vector<double> fPositionX;
std::vector<double> fPositionY;
std::vector<double> fFitADCY;
std::vector<double> vchan;
std::vector<int> fFitBin;
std::vector<double> fFitADC;
std::vector<double> fFitADCError;
std::vector<double> fSpeNumber;
std::vector<double> fSpeNumberError;
std::vector<double> sipmSavedGain;
std::vector<double> sipmSavedGainError;
std::vector<double> sipmGain;
std::vector<double> sipmGainError;
std::vector<double> absoluteGain;
std::vector<double> absoluteGainError;
std::vector<double> sipmNumber;
std::vector<double> sipmNumberError;
TGraphErrors *gSavedGain;
std::string sdate;

enum
{
  CHANNELS = 14,
  NONSUMCHANNELS = CHANNELS - 1
};

enum
{
  NPEAKS = 4
};

std::vector<TH1D *> hlist;
int colors[11] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 30, 40};

/*  start of code */
string currentDate()
{
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  char output[30];
  strftime(output, 30, "%Y-%m-%d-%H-%M", timeinfo);
  return string(output);
}

bool readGains(TString fileName)
{
  TFile *gfin = new TFile(fileName, "readonly");
  if (gfin->IsZombie())
  {
    cout << "Error opening file" << fileName << endl;
    return false;
  }
  cout << " opened sipm gain file " << fileName << endl;
  gSavedGain = NULL;
  gfin->GetObject("gGain", gSavedGain);
  if (gSavedGain == NULL)
  {
    cout << "no gGain in file " << endl;
    return false;
  }
  cout << "found graph named " << gSavedGain->GetName() << " in file " << fileName << endl;
  sipmSavedGain.clear();
  sipmSavedGainError.clear();
  sipmSavedGain.resize(NONSUMCHANNELS);
  sipmSavedGainError.resize(NONSUMCHANNELS);
  for (unsigned long j = 0; j < sipmGain.size(); ++j)
  {
    sipmSavedGain[j] = nominalGain;
    sipmSavedGainError[j] = sqrt(nominalGain);
  }
  for (int i = 0; i < gSavedGain->GetN(); ++i)
  {
    int index = int(gSavedGain->GetPointX(i));
    sipmSavedGain[index] = gSavedGain->GetPointY(i);
    sipmSavedGainError[index] = gSavedGain->GetErrorY(i);
  }

  printf("\t\t\t stored gains %lu \n", sipmSavedGain.size());
  for (unsigned long j = 0; j < sipmSavedGain.size(); ++j)
  {
    printf(" %lu  saved gain %.4f error %.4f   \n", j, sipmSavedGain[j], sipmSavedGainError[j]);
  }
  return true;
}

// Define a linear fit
double fline(double *x, double *par)
{
  return par[0] + x[0] * par[1];
}

/* main */
void gain()
{
  sdate = currentDate();
  printf(" making cans on %s \n", sdate.c_str());
  TString tag0;
  TString tag1;
  TString tagDate;
  // tag0=TString("01_25_2024");
  // tag1=TString("01_25_2024");
  tag0 = TString("09_01_2024");
  tag1 = TString("09_04_2024");
  tagDate = TString("2024-09-06-13-17");
  int nfiles = 41;
  TString tagNfiles;
  tagNfiles.Form("-nfiles-%i", nfiles);

  TString fileName = TString("summary-") + tag0 + TString("-") + tag1 + tagNfiles + TString("-created-") + tagDate + TString(".root");
  cout << " gains from file " << fileName << endl;

  // TH1D *hFit = new TH1D("GausFit", "GausFit", 4000, -20.E3, 200.E3);
  // TFile *fin = new TFile("summary-12_28_2023-01_25_2024-nfiles-45-created-2024-06-27-14-43.root", "readonly");
  TFile *fin = new TFile(fileName, "readonly");
  if (fin->IsZombie())
    return;
  std::string sdate = currentDate();
  TFile *fout = new TFile(Form("gains-%s.root", sdate.c_str()), "recreate");
  TF1 *line = new TF1("myLine", fline, 0, 2.E5, 2);

  // read old gain file
  TString savedGainTag = TString("gains-2024-02-15-17-26-save");
  TString gainFileName = TString(getenv("BOBJ")) + savedGainTag + TString(".root");
  printf("read gains from file %s \n", gainFileName.Data());
  if (!readGains(gainFileName))
  {
    printf("no gain file %s so exit \n", gainFileName.Data());
    exit(0);
  }
  gSavedGain->GetHistogram()->GetYaxis()->SetTitle("absolute saved gain ");
  gSavedGain->GetHistogram()->GetXaxis()->SetTitle("channel");
  gSavedGain->SetName(savedGainTag);
  gSavedGain->SetTitle(Form("saved gains %s ", savedGainTag.Data()));
  fout->Append(gSavedGain);

  // get histos from file
  TDirectory *gainSumDir = nullptr;
  fin->GetObject("gainSumDir", gainSumDir);

  //
  TIter next(gainSumDir->GetListOfKeys());
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());

    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    TString hname(h->GetName());
    int ichan = TString(hname(hname.Last('n') + 1, hname.Length())).Atoi();
    // struggle to get only last hist cycle
    TString lastName;
    if (hlist.size() > 0)
      lastName = hlist[hlist.size() - 1]->GetName();
    // cout << " xxx " << h->GetName() << " " << lastName << endl;
    if (TString(h->GetName()) == lastName)
      continue;
    if (ichan > 11)
      continue;
    if (ichan == 5)
      continue;
    cout << "hist name " << h->GetName() << " chan " << ichan << " cycle " << key->GetCycle() << endl;
    hlist.push_back(h);
    vchan.push_back(ichan);
  }
  printf(" have %lu gain hists \n", hlist.size());

  // fit peaks
  vfit.resize(hlist.size());

  double width = nominalGain / 2.0;
  unsigned i = 4;
  sipmGain.clear();
  sipmGainError.clear();
  for (unsigned i = 0; i < hlist.size(); ++i)
  {
    fFitADC.clear();
    fFitADCY.clear();
    fFitADCError.clear();
    fSpeNumber.clear();
    fSpeNumberError.clear();

    for (unsigned j = 0; j < NPEAKS; ++j)
    {
      TH1D *h = hlist[i];
      double nominalPeak = nominalGain * double(j + 1);
      double fitStart = nominalGain * double(j + 1) - width;
      double fitEnd = nominalGain * double(j + 1) + width;
      int bin = h->FindBin(nominalPeak);
      double xbin = h->GetBinLowEdge(bin);
      double peakVal = h->GetBinContent(bin);
      if (peakVal < 3)
        continue;
      // fit here
      h->GetListOfFunctions()->Clear();
      // printf("fit to hist %s peak # %u bin %i xbin %f val %f from %f to %f \n",
      //     h->GetName(), j+1, bin,xbin,peakVal,fitStart,fitEnd);
      h->Fit("gaus", "QQSS", "", fitStart, fitEnd);
      TF1 *gFit = (TF1 *)h->GetListOfFunctions()->FindObject("gaus");
      if (gFit != nullptr)
      {
        vfit[i].push_back(gFit);
        // scale to nominal
        double mean = (gFit->GetParameter(1) - nominalGain) / nominalGain;
        double meanError = gFit->GetParError(1) / nominalGain;
        int jbin = h->FindBin(mean);
        double val = h->GetBinContent(jbin);
        printf("fit to hist %s point %u  nominal bin %i nominal x %.2f from %.0f to %.0f  mean %.2f error %.2f peak bin %i val %.2f  \n",
               h->GetName(), j + 1, bin, xbin, fitStart, fitEnd, mean, meanError, jbin, val);
        fFitADC.push_back(mean);
        fFitADCY.push_back(val);
        fFitADCError.push_back(meanError);
        fSpeNumber.push_back(j + 1);
        fSpeNumberError.push_back(0);
        // printf("\t fit to hist %i point %lu  (%f,%f) bin %i x %f xadc %f  \n", i, j,fPositionX[j],fPositionY[j],bin,
        // fFitADC[j],fFitADCY[j]);
      }
    } // loop over peaks

    printf(" fit for channel %i n= %lu \n", i, fFitADC.size());
    for (unsigned long j = 0; j < fFitADC.size(); ++j)
      printf("  point %lu  ADC %.2f +/- %.2f  y %.2f\n",
             j, fFitADC[j], fFitADCError[j], fFitADCY[j]);

    if (fFitADC.size() < 2)
      continue;
    TGraphErrors *g = new TGraphErrors(fFitADC.size(), &fSpeNumber[0], &fFitADC[0], &fSpeNumberError[0], &fFitADCError[0]);
    // g->Print();
    g->SetMarkerStyle(23);
    g->SetMarkerColor(kRed);
    g->SetMarkerSize(1.3);

    // add marker
    TPolyMarker *pmold = (TPolyMarker *)hlist[i]->GetListOfFunctions()->FindObject("TPolyMarker");
    if (pmold)
    {
      hlist[i]->GetListOfFunctions()->Remove(pmold);
      delete pmold;
    }
    TPolyMarker *pm = new TPolyMarker(fFitADC.size(), &fFitADC[0], &fFitADCY[0]);
    hlist[i]->GetListOfFunctions()->Add(pm);
    double *xp = pm->GetX();
    double *yp = pm->GetY();
    for (int ip = 0; ip < pm->GetN(); ++ip)
      printf("det %i poly %i %f %f \n", i, ip, xp[ip], yp[ip]);
    pm->SetMarkerStyle(23);
    pm->SetMarkerColor(kRed);
    pm->SetMarkerSize(1.5);

    // make graph

    TString gname;
    gname.Form("gainChan%i", i);
    TString gtitle;
    gtitle.Form(" gain channel  %i ;  number SPE ; ADC/SPE", i);
    g->SetName(gname);
    g->SetTitle(gtitle.Data());
    line->SetParameters(0.5, 0);
    g->Fit("myLine");
    // g->GetHistogram()->GetListOfFunctions()->ls();
    TF1 *gFit = g->GetFunction("myLine");
    if (gFit == nullptr)
      continue;
    TCanvas *gcan = new TCanvas(Form("GainMarkerChan%i", i), Form("chan%i", i));
    gPad->SetLogy(0);
    gStyle->SetOptFit();
    // g->GetHistogram()->GetXaxis()->SetRangeUser(0, 4);
    // g->GetHistogram()->GetYaxis()->SetRangeUser(0, 1.5E5);
    gFit->SetLineStyle(5);
    gFit->SetLineWidth(5);
    g->Draw("APE1");
    gFit->Draw("same");
    gPad->SetGrid();
    gcan->Print(".png");
    fout->Add(g);
    if (gFit->GetParameter(1) < 0)
      continue;
    sipmGain.push_back(gFit->GetParameter(1));
    sipmGainError.push_back(gFit->GetParError(1));
    sipmNumber.push_back(vchan[i]);
    sipmNumberError.push_back(0);

    // plot
    TCanvas *can = new TCanvas(Form("PeaksChan%i", i), Form("chan%i", i));
    gPad->SetLogy(1);
    hlist[i]->Draw();
    fout->Add(hlist[i]);
    can->Print(".png");
    // fout->Add(can);
  } // llop over hists

  // final plots
  for (unsigned i = 0; i < hlist.size(); ++i)
    printf(" good fits to %s = %lu \n",
           hlist[i]->GetName(), vfit[i].size());

  // plot
  for (unsigned i = 0; i < hlist.size(); ++i)
  {
    hlist[i]->GetListOfFunctions()->Clear();
    TString cname;
    cname.Form("Gain%s", hlist[i]->GetName());
    TCanvas *c = new TCanvas(cname, cname);
    gPad->SetLogy(1);
    hlist[i]->Draw();
    for (unsigned ifit = 0; ifit < vfit[i].size(); ++ifit)
      vfit[i][ifit]->Draw("same");
    cout << " printing " << c->GetName() << endl;
    c->Print(".pdf");
    fout->Add(c);
  }

  // plot
  TString cname;
  cname.Form("GainAll");
  TCanvas *c2 = new TCanvas(cname, cname);
  gStyle->SetOptStat(0);
  auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
  for (int i = hlist.size() - 1; i >= 0; --i)
  {
    cout << " plot " << hlist[i]->GetName() << endl;
    legend->AddEntry(hlist[i], Form("%s", hlist[i]->GetName()), "f");
    gPad->SetLogy(1);
    hlist[i]->SetLineColor(colors[i]);
    if (i == hlist.size() - 1)
      hlist[i]->Draw();
    else
      hlist[i]->Draw("sames");
  }
  legend->Draw();
  fout->Add(c2);

  // plot
  TString graphTitle;
  graphTitle = Form("gain variation %s to %s ", tag0.Data(), tag1.Data());
  TGraphErrors *gGain = new TGraphErrors(sipmGain.size(), &sipmNumber[0], &sipmGain[0], &sipmNumberError[0], &sipmGainError[0]);
  gGain->SetName("gGain");
  gGain->SetTitle(graphTitle);
  gGain->SetMarkerStyle(23);
  gGain->SetMarkerColor(kRed);
  gGain->SetMarkerSize(1.3);
  gGain->GetHistogram()->GetYaxis()->SetTitle(Form("gain as fraction of nominal %.2f ", nominalGain));
  gGain->GetHistogram()->GetXaxis()->SetTitle("channel");

  TString canTitle;
  TString canName;
  canName = Form("gainVariation-%s-%s", tag0.Data(), tag1.Data());
  canTitle = Form("gain variation %s to %s ", tag0.Data(), tag1.Data());
  TCanvas *canGain = new TCanvas(canName, canTitle);
  gPad->SetGrid();
  gGain->Draw("AP");
  cout << " printing " << canGain->GetName();
  canGain->Print(".png");

  // fill absolute
  absoluteGain.resize(sipmSavedGain.size());
  absoluteGainError.resize(sipmSavedGain.size());
  for (unsigned long j = 0; j < sipmSavedGain.size(); ++j)
  {
    absoluteGain[j] = sipmGain[j] * sipmSavedGain[j];
    absoluteGainError[j] = sipmGainError[j];
  }

  TString absoluteTitle;
  absoluteTitle = Form("absolute gain  %s to %s ", tag0.Data(), tag1.Data());
  TGraphErrors *absGain = new TGraphErrors(absoluteGain.size(), &sipmNumber[0], &absoluteGain[0], &sipmNumberError[0], &absoluteGainError[0]);
  TString absoluteName;
  absoluteName = Form("gains-%s-%s", tag0.Data(), tag1.Data());
  absGain->SetName("absoluteName");
  absGain->GetHistogram()->GetYaxis()->SetTitle("absolute gain");
  absGain->GetHistogram()->GetXaxis()->SetTitle("channel");
  absGain->SetTitle(absoluteTitle);
  absGain->SetMarkerStyle(23);
  absGain->SetMarkerColor(kRed);
  absGain->SetMarkerSize(1.3);
  absGain->GetHistogram()->GetYaxis()->SetTitle("gain as fraction of nominal");
  absGain->GetHistogram()->GetXaxis()->SetTitle("channel");

  fout->Add(absGain);

  canName = Form("absoluteGain-%s-%s", tag0.Data(), tag1.Data());
  canTitle = Form("absolute gain  %s to %s ", tag0.Data(), tag1.Data());
  TCanvas *canAbsGain = new TCanvas(canName, canTitle);
  gPad->SetGrid();
  absGain->Draw("AP");
  cout << " printing " << canGain->GetName();
  canAbsGain->Print(".png");

  canName = Form("savedGain-%s", savedGainTag.Data());
  canTitle = Form("saved gain file %s ", savedGainTag.Data());
  TCanvas *canSavedGain = new TCanvas(canName, canTitle);
  gPad->SetGrid();
  gSavedGain->Draw("AP");
  cout << " printing " << canSavedGain->GetName();
  canSavedGain->Print(".png");

  fout->Write();
}