//  fit singelet to landau
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "modelFit.hh"
std::string sdate;
// time is in microseconds
using namespace TMath;
TFile *fin;
TFile *fout;
bool fixTriplet = true;
bool restrictive = false;

std::ofstream dataFile;
std::ofstream inteFile;
std::ofstream fitFile;

std::vector<double> vchan;
std::vector<double> vchanErr;
std::vector<double> vinte;
std::vector<double> vinteCorr;
std::vector<double> veff;
std::vector<double> vnorm1;
std::vector<double> vnorm1Err;
std::vector<double> vabs;
double fbackStart = 8000.; // ns
double fbackEnd = 15000.;  // ns

double effAndFill = SiPMQE128Ham * fillFactor;
double fitLow, fitHigh;

TDirectory *runSumDir;
TString fileName;
TCanvas *canModel;
double ppm = 0;
double binwidth = theBinWidth;
double tiny = 1.0E-11;
int colors[NTYPES] = {1, 3, 2, 4, 7, kMagenta};
int chanMap[12] = {0, 0, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9};

// 3603795;
modelFit *models[NTYPES];
TH1D *hmodel[NTYPES];
TF1 *gback = NULL;

TF1 *fitBack(TH1D *hist, double &ave, double &cnst, double &tau)
{
  TF1 *gfit = NULL;
  int backStart = hist->FindBin(fbackStart);
  int backEnd = hist->FindBin(fbackEnd);
  TF1 *funcBack = new TF1("funcBack", "[0]*TMath::Exp(-x/[1])", fbackStart, fbackEnd);
  funcBack->SetParameter(0, 1E-5);
  funcBack->SetParameter(1, 1E4);
  funcBack->SetParLimits(0, 1.E-9, 1E2);
  funcBack->SetParName(0, "const");
  funcBack->SetParName(1, "tau [ns]");

  hist->Fit("funcBack", "LF0+", " ", fbackStart, fbackEnd);
  gfit = (TF1 *)hist->GetListOfFunctions()->FindObject("funcBack");
  if (gfit)
  {
    cnst = gfit->GetParameter(0);
    ave = gfit->Eval(fbackEnd);
    tau = gfit->GetParameter(1);
  }
  else
    printf("funcBack Fit to hist fails \n");
  double aveb = hist->Integral(backStart, backEnd) / double(backEnd - backStart);
  printf("FITBACK: %s npar %i ave %E aveb %E \n\n", gfit->GetName(), gfit->GetNpar(), ave, aveb);
  return gfit;
}

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

void errorProp(double a, double ea, double b, double eb, double &c, double &ce)
{
  c = a * b;
  ce = c * sqrt(pow(ea / a, 2.) + pow(eb / b, 2.));
}

void setModels(int setModel, int ichan)
{
  // set all model parameters
  for (int j = 0; j < NTYPES; ++j)
  {
    printf("setModels setModel %i for type %i %s chan %i \n", setModel, j, modelNames[j].Data(), models[j]->theChan);
    for (int ii = 0; ii < NPARS; ++ii)
    {
      models[j]->fp->SetParameter(ii, models[setModel]->fp->GetParameter(ii));
    }
    models[j]->fp->SetParameter(TYPE, j);
    /*for (int ii = 0; ii < NPARS; ++ii)
    {
      printf(" ... %i %s par %f \n", ii, models[setModel]->fp->GetParName(ii), models[j]->fp->GetParameter(ii));
    }
    */
  }
}

// plot models using graphs
void plotModels(int setModel, int ichan, TH1D *hWave)
{
  // set all model parameters
  setModels(setModel, ichan);
  printf("gback pars %i \n", gback->GetNpar());

  double off = 1.E-4; // for log plotting
  // clone histograms
  for (int i = 0; i < NTYPES; ++i)
  {
    TString graphName;
    graphName.Form("modelSet%sChan%i%s", modelNames[setModel].Data(), ichan, modelNames[i].Data());
    // hWave->GetListOfFunctions()->Clear();
    hmodel[i] = (TH1D *)hWave->Clone(graphName);
    hmodel[i]->Reset("ICE");
    hmodel[i]->SetTitle(graphName);
    // printf("created %s gback pars %i \n", hmodel[i]->GetName(), gback->GetNpar());
    //  hmodel[i] = new TH1D(graphName, graphName, hWave->GetNbinsX(), hWave->GetBinLowEdge(0), hWave->GetBinLowEdge(hWave->GetNbinsX()));

    for (int ibin = startBin; ibin < endBin; ++ibin)
    {
      double val = models[i]->fp->Eval(hWave->GetBinCenter(ibin));

      if (val > 0)
        val = TMath::Max(val, 1.E-12); // avoid very small values
      // if (isnan(val))
      // if (ibin > 600 && ibin < 750)
      // printf("............. model %i bin %i xbin %f val %E !!!\n", i, ibin, hWave->GetBinCenter(ibin), val);

      if (!isnan(val))
      {
        hmodel[i]->SetBinContent(ibin, val);
        hmodel[i]->SetBinError(ibin, 0);
        // if (i == 0 && val > off)
        // printf("............. %i bin %i xbin %f  val %E \n", i, ibin, hWave->GetBinCenter(ibin), val);
      }
    }

    TString graphTitle;
    graphTitle.Form("modelset%stype%s", modelNames[setModel].Data(), modelNames[i].Data());
    // hmodel[i]->GetXaxis()->SetRangeUser(startTime, 4000.);
    // hmodel[i]->GetYaxis()->SetRangeUser(10.E-5, 1.E-1);
    hmodel[i]->GetYaxis()->SetMoreLogLabels();

    hmodel[i]->GetXaxis()->SetTitle("time [ns] ");
    hmodel[i]->GetYaxis()->SetTitle(" yield [SPE]");
    hmodel[i]->SetTitle(graphTitle);
    hmodel[i]->SetLineColor(colors[i]);
    hmodel[i]->SetMarkerColor(colors[i]);
    hmodel[i]->SetMarkerStyle(20 + i);
    hmodel[i]->SetMarkerSize(0.2);
    fout->Add(hmodel[i]);
    printf("PPPPPPP model %i chan %i histo %s \n", i, models[i]->theChan, hmodel[i]->GetName());
  }
  TString canName;
  canName.Form("plotModels%sChan%iStart-PPM-%.3f", modelNames[setModel].Data(), ichan, ppm);
  canModel = new TCanvas(canName, canName);
  canModel->SetLogy();
  hmodel[NTYPES - 1]->Draw("HIST"); // model all
  int LAST = NTYPES - 1;
  if (setModel == MODELSINGLET)
    LAST = 3;

  for (int j = 0; j < LAST; ++j)
    hmodel[j]->Draw("HISTsame");
  canModel->BuildLegend();
}

// plot models using hWave clones
void printModels(TH1D *hWave)
{
  int maxBin = hWave->GetNbinsX();
  for (int ibin = 700; ibin < 710; ++ibin)
    models[MODELALL]->fp->Eval(hWave->GetBinCenter(ibin));
}

bool openFile()
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

  fin = new TFile(fileName, "readonly");
  printf(" opened file %s\n", fileName.Data());
  runSumDir = NULL;
  fin->GetObject("runSumDir", runSumDir);
  if (!runSumDir)
  {
    printf(" no runSumDir in file %s\n", fileName.Data());
  }

  // output ascii file append
  TString dataFileName;
  sdate = currentDate();
  dataFileName.Form("singletNorm.dat");
  dataFile.open(dataFileName.Data(), std::ios_base::app);
  if (!dataFile)
  {
    std::cout << dataFileName << " failed to open. Closing " << std::endl;
  }

  // output ascii file append
  TString fitFileName;
  sdate = currentDate();
  fitFileName.Form("fitResults.dat");
  fitFile.open(fitFileName.Data(), std::ios_base::app);
  if (!fitFile)
  {
    std::cout << fitFileName << " failed to open. Closing " << std::endl;
  }

  // remake
  TString inteFileName;
  inteFileName.Form("singletInte.dat");
  inteFile.open(inteFileName.Data(), std::ios::out | std::ios::in | std::ios::trunc);
  if (!inteFile)
  {
    std::cout << inteFileName << " failed to open. Closing " << std::endl;
    return false;
  }
  // runSumDir->ls();/

  // get integrals
  TIter next(runSumDir->GetListOfKeys());
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    TString name = TString(h->GetName());
    if (!name.Contains("RunPeakWave"))
      continue;
    std::string sname(name.Data());
    std::string chan = sname.substr(sname.find_last_of("n") + 1);
    int ichan = std::atoi(chan.c_str());
    if (ichan > 11)
      break;
    vchan.push_back(ichan);
    vchanErr.push_back(0.);
    vinte.push_back(h->Integral(startBin, singletEndBin));
    double dist = distanceLevel[level(ichan)];
    double abs = Absorption(ppm, dist);
    double effgeo = effGeoFunc(ichan); // dont forget
    vabs.push_back(abs);
    veff.push_back(effgeo);
    printf("!!!!PUSH BACK chan %i level %i dist %.3f eff %.3E abm %.3E eff*ab %.3E total SPE %.3E \n ", ichan, level(ichan), dist, effgeo, abs, effgeo * abs, vinte[vinte.size() - 1]);
  }

  return true;
}

void tbFit1(unsigned ichan = 8, TString flag = "A", double thePPM = 2.)
{
  ppm = thePPM;
  TString fileName0 = TString("summary-11_28_2023-01_07_2024-nfiles-349-created-2024-07-25-15-45.root");
  TString fileName2 = TString("summary-03_15-2024-03_28_2024-nfiles-125-created-2024-08-06-12-34.root");
  fileName = fileName2;

  // make individual light curves
  for (int ifit = 0; ifit < NTYPES; ++ifit)
  {
    models[ifit] = new modelFit(ifit, ichan, ppm);
  }

  if (!openFile())
    return;

  printf("tfFit1: opened file %s chan %i \n", fileName.Data(), ichan);

  for (int ichan = 0; ichan < 12; ++ichan)
  {
    int ilevel = -1;
    if (ichan == 6 || ichan == 7 || ichan == 8 || ichan == 9 || ichan == 10 || ichan == 11)
      ilevel = 0;
    else if (ichan == 3 || ichan == 4 || ichan == 5)
      ilevel = 1;
    else if (ichan == 0 || ichan == 1 || ichan == 2)
      ilevel = 2;
    else if (ichan == 12)
      ilevel = 3;
    printf(" chan %i level %i dist %.3f eff %.3E \n ", ichan, ilevel, distanceLevel[ilevel], effGeoFunc(ichan));
  }

  vinteCorr.resize(vinte.size());
  printf(">>>> channels %lu  %lu  %lu \n", vchan.size(), vinte.size(), veff.size());
  TString inteElement;
  for (int jc; jc < vchan.size(); ++jc)
  {
    vinteCorr[jc] = vinte[jc] / veff[jc] / vabs[jc] / effAndFill;
    // vinteCorr[jc] = vinte[jc]/veff[jc]/SiPMQE128Ham;
    printf(" chan %i int %.3E eff corrected int %.3E eff %.3E vabs %.3f \n", int(vchan[jc]), vinte[jc], vinteCorr[jc], veff[jc], vabs[jc]);
    inteElement.Form(" %E %E  %i \n", vinte[jc], vinteCorr[jc], int(vchan[jc]));
    inteFile << inteElement.Data();
  }
  inteFile.close();

  printf(">>>>  graph #points %lu \n", vchan.size());
  TGraph *gnorm = new TGraph(vchan.size(), &vchan[0], &vinte[0]);
  gnorm->SetTitle("integral");
  gnorm->GetXaxis()->SetTitle("channel ");
  gnorm->GetYaxis()->SetTitle(" singlet integral");
  gnorm->SetMarkerStyle(21);
  // gnorm->GetYaxis()->SetRangeUser(0,80);

  TCanvas *cnorm1 = new TCanvas("integral", "integral");
  gnorm->Draw("apm");

  TGraph *gnorme = new TGraph(vchan.size() - 2, &vchan[0], &vinteCorr[0]);
  gnorme->SetTitle("eff corrected integral");
  gnorme->GetXaxis()->SetTitle("channel ");
  gnorme->GetYaxis()->SetTitle(" correctd singlet integral");
  gnorme->SetMarkerStyle(21);
  gnorme->SetMarkerColor(kBlue);
  gnorme->SetMarkerSize(1.);
  // gnorme->GetYaxis()->SetRangeUser(1E4,4E+05);
  // gnormf->GetYaxis()->SetRangeUser(50,2000);

  fout = new TFile("tbFit.root", "recreate");
  TString hname;
  hname.Form("RunPeakWaveChan%i", ichan);
  // get lightcurve to fit
  TH1D *hWaveIn = NULL;
  runSumDir->GetObject(hname, hWaveIn);
  if (!hWaveIn)
    printf("no hist %s \n", hname.Data());
  else
    printf("got %s \n", hWaveIn->GetName());

  // make some histograms for later
  TH1D *hWave = (TH1D *)hWaveIn->Clone();
  hWave->GetListOfFunctions()->Clear();
  /* */
  TH1D *hSinglet = (TH1D *)hWaveIn->Clone(Form("SingletChan%i", ichan));
  hSinglet->SetTitle(Form("SingletChan%i", ichan));
  TH1D *hSingletFit = (TH1D *)hWaveIn->Clone(Form("SingletFitChan%i", ichan));
  hSingletFit->SetTitle(Form("SingletFitChan%i", ichan));
  TH1D *hWaveTriplet = (TH1D *)hWaveIn->Clone(Form("Triplet%s", hWaveIn->GetName()));
  hWaveTriplet->SetTitle(Form("WaveTripletChan%i", ichan));
  // TH1D *hFullFit = (TH1D *)hWaveIn->Clone(Form("FullFitChan%i", ichan));
  // hFullFit->SetTitle(Form("FullFitChan%i", ichan));
  //  double xWaveUpper = hWave->GetBinLowEdge(hWave->GetNbinsX()) + hWave->GetBinWidth(0);
  //  TH1D *hSingletFit = new TH1D("hSingletFit", "Wave Singlet Fit", hWave->GetNbinsX(), 0, xWaveUpper);
  //  TH1D *hFullFit = new TH1D("hFullFit", "Wave Full Fit", hWave->GetNbinsX(), 0, xWaveUpper);

  // fout->Add(hWave);

  // plot with defaults
  models[MODELALL]->show();
  // plotModels(ichan, hWave);
  //  PPPPPPP just single

  // ******* fit singlet first *****
  TF1 *fp1 = models[MODELSINGLET]->fp;
  double bAve;
  double bCnst;
  double bTau;
  gback = fitBack(hWave, bAve, bCnst, bTau);
  // gback->Print();
  printf("background chan %i ave %E tau %E  over %f ns  rate %E kHz \n", ichan, bAve, bTau, fbackEnd - fbackStart, bAve * 1.E6);
  printf("BKG %s ichan=%i; back=%E; rate=%E; tau=%E; \n", gback->GetName(), ichan, bAve, bAve * 1.E6, bTau);

  TCanvas *bcan = new TCanvas("backFit", "backFit");
  hWave->GetXaxis()->SetTitle("time [ns] ");
  hWave->GetYaxis()->SetTitle(" yield/2ns [SPE]");
  gStyle->SetOptFit(1111);
  bcan->SetLogy();
  hWave->Draw();
  gback->Draw("same");
  bcan->SetGrid();
  bcan->Print(".pdf");
  // return;

  // hWave->GetListOfFunctions()->Clear();

  Double_t binwidth = hWave->GetBinWidth(1); // ns
  int maxBin = hWave->GetMaximumBin();
  double mpv = hWave->GetBinCenter(maxBin);
  double effGeo = effGeoFunc(ichan);
  double detNorm1 = hWave->Integral(startBin, singletEndBin);
  double detNorm3 = hWave->Integral(singletEndBin, endBin);
  double norm1Start = detNorm1 / effGeo / effAndFill;
  double norm3Start = detNorm3 / effGeo / effAndFill;

  double effChan = effGeoFunc(ichan);
  // set parameters
  double normStart = norm1Start + norm3Start;
  double sfracStart = norm1Start / normStart;
  // fp1->SetParameter(NORM, normStart);
  // double tnorm = 1./effChan/SiPMQE128Ham/sfracStart;
  double tnorm = normStart;

  printf(">>>> setting normStart %E <<<<< \n", normStart);

  fp1->SetParameter(NORM, normStart);
  fp1->SetParLimits(NORM, 1.0E-3 * normStart, 1.0E1 * normStart);
  fp1->SetParameter(MPV, mpv);
  fp1->SetParLimits(MPV, 0.9 * mpv, 1.1 * mpv);
  fp1->SetParameter(SINGLETSIGMA, singletSigma0);
  fp1->SetParLimits(SINGLETSIGMA, 0.1 * singletSigma0, 100. * singletSigma0);
  // fp1->FixParameter(ARECOMB, 100 * Arecomb0);
  // fp1->SetParLimits(TRECOMB, Trecomb0 / 10., 10. * Trecomb0);
  if (ichan == 3)
  {
    fp1->SetParLimits(SINGLETTAU, 0.1 * singletTau0, 2. * singletTau0);
    // fp1->FixParameter(SINGLETTAU, singletTau0);
  }

  // fix parametgers
  // fp1->FixParameter(EXPFRAC, 0.); /// fix or set?
  if (ichan == 1) // just landau
  {
    fp1->FixParameter(RECOMBPOWER, RecombPower0); // not fitting
    fp1->FixParameter(ARECOMB, 0);
    fp1->FixParameter(TRECOMB, Trecomb0);
    fp1->FixParameter(SINGLETTAU, singletTau0); // not fitting
    fp1->FixParameter(EXPFRAC, 0.0);            // all Landau
    fp1->FixParameter(SFRAC, sfrac0);
  }
  else if (ichan == 4)
  {
    fp1->FixParameter(SINGLETSIGMA, singletSigma0); // not fitting
    fp1->SetParameter(SINGLETTAU, singletTau0);     // not fitting
    fp1->SetParLimits(SINGLETTAU, 0.1 * singletTau0, 100. * singletTau0);
    // fp1->FixParameter(SINGLETTAU, singletTau0); // not fitting
    //  fp1->FixParameter(MPV, mpv);
    fp1->FixParameter(EXPFRAC, 1.0); // just Exp
    fp1->FixParameter(SFRAC, sfrac0);
    fp1->SetParameter(ARECOMB, Arecomb0);
    // fp1->SetParLimits(ARECOMB, 1.E-1 * Arecomb0, 10 * Arecomb0);
    fp1->SetParLimits(RECOMBPOWER, 1.E-4, 1.E2);
  }
  else
  {
    fp1->FixParameter(SINGLETSIGMA, singletSigma0); // not fitting
    fp1->FixParameter(RECOMBPOWER, RecombPower0);
    fp1->FixParameter(ARECOMB, Arecomb0);
    // fp1->SetParLimits(ARECOMB, 0, 1000. * Arecomb0);
    fp1->FixParameter(EXPFRAC, 1.0);
    // fp1->SetParameter(EXPFRAC, expFrac0);
    fp1->FixParameter(TRECOMB, Trecomb0);
    // fp1->SetParameter(EXPFRAC, expFrac0);
    // fp1->SetParLimits(EXPFRAC, 1.0E-9, 1.);
    // fp1->FixParameter(EXPFRAC, 1.0); // just Exp
    // fp1->SetParameter(SINGLETTAU, singletTau0); /// fix or set?
    fp1->SetParLimits(SINGLETTAU, 0.1 * singletTau0, 2.0 * singletTau0);
  }
  fp1->FixParameter(PPM, ppm); /// fix or set?
  // fp1->SetParameter(SFRAC, sfracStart);
  fp1->FixParameter(SFRAC, sfrac0);
  fp1->FixParameter(BKGCONST, bCnst);
  fp1->FixParameter(BKGTAU, bTau);
  fp1->FixParameter(AB, vabs[ichan]);
  fp1->FixParameter(TAU3, tTriplet0);
  fp1->FixParameter(KQ, kq0);
  fp1->FixParameter(KDIFFUSION, kdiffusion0); // this is starting value Eur. Phys. J. C (2013) 73:2618
  fp1->FixParameter(TAUM, tauM0);

  // fix parameters not in singlet fit
  // fp1->SetParLimits(PPM,0.0,100* ppm);
  // fp1->SetParLimits(SINGLETSIGMA, .001*singletSigma0, 100*singletSigma0);

  models[MODELSINGLET]->show();
  printf("MMMMMMMMMMMMM chan %i maxBin %i mpv %f \n", ichan, maxBin, mpv);
  printf(" \n****** singletStart %.3f bin %i singletEnd %.3f bin %i last bin %i\n", startTime, startBin, singletEndTime, singletEndBin, endBin);
  printf("\n******   MPV-nominalTrigger  %.1f effGeo %.3E detNorm3 %.3E norm1  %.3E norm3 %.3E  sfracStart %.3f *****\n", mpv, effGeo, detNorm3, norm1Start, norm3Start, sfracStart);

  // histogram for peak singlet fit
  hSinglet->SetTitle(Form("SingletWave%i", ichan));
  hSinglet->SetMarkerColor(kBlue);
  hSinglet->SetLineColor(kBlue);

  /* do the fit here L=likelihood M=improved E=better errors 0 = do not draw S return full results + add fit to hist */
  double singletFitEnd = singletEndTime; // 3000.;
  if (ichan == 4)
    singletFitEnd = 6000.;

  double singletEndPlot = 1500; // 3000.;
  if (ichan == 1 || ichan == 4)
    singletEndPlot = 3000.;

  fp1->SetRange(startTime, singletFitEnd);
  TFitResultPtr fptr;

  if (ichan == 1 || ichan == 4)
    fptr = hSinglet->Fit(fp1, "RS0+", "", startTime, singletFitEnd);
  else
    fptr = hSinglet->Fit(fp1, "WLRS0+", "", startTime, singletFitEnd);
  int fitResult = fptr;
  if (fitResult != 0)
  {
    printf("!!!!!! fit to singlet fails %i  status = migradStatus + 10*minosStatus + 100*hesseStatus + 1000*improveStatus \n", fitResult);
    models[MODELSINGLET]->show();
    // still fill singlet
    hSingletFit->SetTitle(Form("SingletFitChan%i", ichan));
    hSingletFit->SetMarkerColor(colors[0]);
    hSingletFit->SetLineColor(colors[0]);
    fout->Add(hSingletFit);
    for (int ibin = 0; ibin < hSingletFit->GetNbinsX(); ++ibin)
    {
      double val = fp1->Eval(hSingletFit->GetBinCenter(ibin));
      if (val < tiny || isnan(val))
        val = 0;
      // printf("PPPPPPP fill models for chan %i type %i \n",ichan, models[MODELSINGLET]->theChan);
      hSingletFit->SetBinContent(ibin, val);
    }
    gStyle->SetOptStat();
    gStyle->SetOptFit(1111);

    hSinglet->GetXaxis()->SetRangeUser(startTime, singletEndPlot);
    hSingletFit->GetXaxis()->SetRangeUser(startTime, singletEndPlot);
    printf("fit singlet in range  %f - %f \n", fitLow, fitHigh);
    plotModels(MODELSINGLET, ichan, hWave);
    TCanvas *canSingletFit = new TCanvas(Form("SingletFailedChan%i-%.3f-PPM", ichan, ppm), Form("SingletFailedFit-%.3f-PPM", ppm));
    canSingletFit->SetLogy(1);
    hSinglet->Draw("HIST");
    hSingletFit->Draw("SAMEHISTC");
    hSingletFit->Draw("SAMEHISTC");
    printf("MMMMMMMMMMMMM chan %i maxBin %i mpv %f \n", ichan, maxBin, mpv);
    fp1->GetRange(fitLow, fitHigh);

    return;
  }
  TMatrixDSym cov = fptr->GetCorrelationMatrix();
  printf(" correlation chan %i cov(2,3) %f \n", ichan, cov(2, 3));
  fptr->Print("V");

  // fill fitted histogran
  hSingletFit->SetMarkerColor(colors[0]);
  hSingletFit->SetLineColor(colors[0]);
  fout->Add(hSingletFit);
  for (int ibin = 0; ibin < hSingletFit->GetNbinsX(); ++ibin)
  {
    double val = fp1->Eval(hSingletFit->GetBinCenter(ibin));
    // if (ibin < 1000)
    //   printf("FFFFF filling hSingletFit bin %i val %E \n", ibin, val);
    if (val < tiny)
      val = 0;
    hSingletFit->SetBinContent(ibin, val);
  }

  // get fit results

  mpv = fp1->GetParameter(MPV);
  ppm = fp1->GetParameter(PPM);
  double norm = fp1->GetParameter(NORM);
  double sfrac = fp1->GetParameter(SFRAC);
  double singletSigma = fp1->GetParameter(SINGLETSIGMA);
  double expFrac = fp1->GetParameter(EXPFRAC);
  double singletTau = fp1->GetParameter(SINGLETTAU);
  double mpvErr = fp1->GetParError(MPV);
  double singletSigmaErr = fp1->GetParError(SINGLETSIGMA);
  double singletTauErr = fp1->GetParError(SINGLETTAU);
  double normErr = fp1->GetParError(NORM);
  double sfracErr = fp1->GetParError(SFRAC);
  double arecomb1 = fp1->GetParameter(ARECOMB);
  double trecomb1 = fp1->GetParameter(TRECOMB);
  double recombPower1 = fp1->GetParameter(RECOMBPOWER);

  // results
  cout << " ----------------------------------------------------------  " << endl;
  models[MODELSINGLET]->show();
  cout << " ----------------------------------------------------------  " << endl;

  double detTotalUn = hWave->Integral(startBin, endBin);
  double detTotalNorm = detTotalUn / veff[ichan] / effAndFill;
  double fitNorm1Err;
  double fitNorm1 = hSingletFit->IntegralAndError(startBin, singletEndBin, fitNorm1Err);
  double fitNorm1Corr = fitNorm1 / veff[ichan] / effAndFill;
  double fitNorm1CorrErr = fitNorm1Err / veff[ichan] / effAndFill;

  double norm1, norm1Err;
  errorProp(norm, normErr, sfrac, sfracErr, norm1, norm1Err);

  // printf("fit to singlet channel # , norm, sfrac  MPV , lsigma , gsigma \n %i & %.2E $\\pm$ %.2E & %.2E $\\pm$ %.2E & %.2E $\\pm$ %.2E & %.2E $\\pm$ %.2E & %.2E $\\pm$ %.2E & %.2E $\\pm$ %.2E \n", ichan, norm, normErr, sfrac0, sfracErr, mpv, mpvErr, singletSigma,singletSigmaErr, singletTau, singletTauErr,expFrac,expFracErr);

  // do a fit to triplet
  // hWave->GetListOfFunctions()->Clear();
  double tripletFitStart = 2000;
  double tripletFitEnd = 4000;
  TF1 *fTriplet = new TF1("ftriplet", "[0]*Exp(-x/[1])", tripletFitStart, tripletFitEnd);
  fTriplet->SetParameter(0, 0.1);
  fTriplet->SetParameter(1, 1000.);
  fTriplet->SetParName(0, "const");
  fTriplet->SetParName(1, "tau");
  fTriplet->Print();

  hWaveTriplet->Fit(fTriplet, "R0+", "", tripletFitStart, tripletFitEnd);
  TF1 *gTriplet = (TF1 *)hWaveTriplet->GetListOfFunctions()->FindObject("ftriplet");
  double cnst = 0;
  double tau3 = 0;
  if (gTriplet)
  {
    cnst = gTriplet->GetParameter(0);
    tau3 = gTriplet->GetParameter(1);
  }
  // *** triplet fit canvas ***
  TCanvas *canTrip = new TCanvas(Form("tripletExpChan%i", ichan), Form("tripletExpChan%i", ichan));
  canTrip->cd();
  // hWave->GetYaxis()->SetMoreLogLabels();
  hWaveTriplet->GetXaxis()->SetRangeUser(startTime, 2. * tripletFitEnd);
  hWaveTriplet->Draw("");
  gTriplet->Draw("same");
  gStyle->SetOptFit(1111);
  canTrip->SetLogy();
  canTrip->SetGrid();
  canTrip->Print(".pdf");

  gStyle->SetOptStat();
  gStyle->SetOptFit(1111);
  // f (ichan == 1 || ichan == 4)
  // singletEndPlot = 3000.;
  hSinglet->GetXaxis()->SetRangeUser(startTime, singletEndPlot);
  hSingletFit->GetXaxis()->SetRangeUser(startTime, singletEndPlot);
  /* ************* */

  fp1->GetRange(fitLow, fitHigh);
  printf("fit singlet in range  %f - %f \n", fitLow, fitHigh);
  plotModels(MODELSINGLET, ichan, hWave);
  models[MODELSINGLET]->show();
  // ****** singlet fit canvas ******
  TCanvas *canSingletFit = new TCanvas(Form("SingletFitChan%i-%.3f-PPM", ichan, ppm), Form("SingletFit-%.3f-PPM", ppm));
  canSingletFit->cd();
  canSingletFit->SetLogy();
  // if (ichan == 1 || ichan == 4)
  canSingletFit->SetLogy();
  // new TPaveStats(0.243553,0.147651,0.6088825,0.4608501,"brNDC");
  hSinglet->Draw("PE0");
  hSingletFit->Draw("SAMEHISTC");
  // hmodel[MODELRECOMB]->Draw("SAMEHISTC");
  // hmodel[MODELALL]->Draw("SAMEHISTC");
  canSingletFit->Update();
  if (ichan != 1 && ichan != 4)
  {
    TPaveStats *s = (TPaveStats *)canSingletFit->GetPrimitive("stats");
    s->SetX1NDC(0.243553);
    s->SetY1NDC(0.147651);
    s->SetX2NDC(0.6002865);
    s->SetY2NDC(0.4608501);
  }
  canSingletFit->Print(".pdf");
  fout->Add(canSingletFit);
  // canSingletFit->Close();

  printf("background chan %i ave %E tau %E  over %f ns  rate %E kHz \n", ichan, bAve, bTau, fbackEnd - fbackStart, bAve * 1.E6);
  printf("BKG %s npar %i ichan=%i; back=%E; rate=%E; tau=%E; \n", gback->GetName(), gback->GetNpar(), ichan, bAve, bAve * 1.E6, bTau);

  printf("TTTTTT tripletFit: from %f to %f ichan %i; cnst=%E ; tau=%E ;\n", tripletFitStart, tripletFitEnd, ichan, cnst, tau3);

  printf("TTTTT ichan %i  normStart %.3E norm %.3E fitNorm1 %.3E +/-  %.3E chan %i norm1Start  %.3E norm1Start+norm3Start %.3E hist integral %.3E  detTotalNorm %.3E tau3 %E sfrac %.3f effChan %.3E \n ",
         ichan, normStart, norm, fitNorm1, fitNorm1Err, int(chanMap[ichan]), norm1Start, norm1Start + norm3Start, detTotalUn, detTotalNorm, tau3, sfrac, effChan);

  TString dataElement;
  dataElement.Form(" %E %E %E %E %E %E %E %E %i \n", norm1, norm1Err, fitNorm1, fitNorm1Err, fitNorm1Corr, fitNorm1CorrErr, sfracStart, tau3, ichan);
  dataFile << dataElement.Data();
  dataFile.close();

  printf("early fit is valid\n");

  // TTTTTT
  if (flag.Contains("T"))
    return;

  /* fit the recomb here */

  /****************
   * *************fit remaining model ************
   * ***************/
  printf(" FIT MODEL ALL \n");
  hWave->GetListOfFunctions()->Clear();
  TF1 *fpAll = models[MODELALL]->fp;
  if (restrictive)
    fpAll->SetRange(5000, endTime); // axis values
  else
    fpAll->SetRange(startTime, endTime); // axis values

  /* set parameters */
  if (fixTriplet)
  {
    fpAll->FixParameter(TAU3, tau3);
    hWaveTriplet->SetName(Form("WaveTriplet%.0fChan%i", tau3, ichan));
    hWaveTriplet->SetTitle(Form("WaveTriplet%.0fChan%i", tau3, ichan));
  }
  else
  {
    fpAll->SetParameter(TAU3, tTriplet0);
    fpAll->SetParLimits(TAU3, 0.5 * tTriplet0, 1.5 * tTriplet0);
  }

  // fpAll->FixParameter(SFRAC, sfracStart);
  if (!restrictive)
  {
    fpAll->SetParameter(NORM, norm);
    fpAll->SetParameter(SFRAC, sfrac);
    // fpAll->SetParLimits(SFRAC, 0., 1.);
    fpAll->FixParameter(TAUM, tauM0);
    fpAll->SetParameter(ARECOMB, arecomb1); // arecomb1
    fpAll->SetParameter(TRECOMB, trecomb1);
    fpAll->SetParameter(RECOMBPOWER, recombPower1);
    fpAll->SetParLimits(RECOMBPOWER, 0., 10.);
    fpAll->FixParameter(SINGLETSIGMA, singletSigma);
    fpAll->FixParameter(SINGLETTAU, singletTau);
    fpAll->FixParameter(EXPFRAC, 1.);
    fpAll->FixParameter(BKGCONST, bCnst);
    fpAll->FixParameter(BKGTAU, bTau);
    fpAll->FixParameter(AB, vabs[ichan]);
    fpAll->FixParameter(PPM, thePPM);
    fpAll->FixParameter(KQ, kq0);
    fpAll->FixParameter(KDIFFUSION, kdiffusion0);
  }
  else
  {
    fpAll->SetParameter(NORM, norm);
    fpAll->FixParameter(MPV, mpv);
    fpAll->FixParameter(SFRAC, sfracStart);
    fpAll->FixParameter(TAUM, tauM0);
    fpAll->FixParameter(ARECOMB, arecomb1);
    fpAll->FixParameter(TRECOMB, trecomb1);
    fpAll->FixParameter(RECOMBPOWER, recombPower1);
    fpAll->FixParameter(SINGLETSIGMA, singletSigma);
    fpAll->FixParameter(SINGLETTAU, singletTau);
    fpAll->FixParameter(EXPFRAC, 0);
    // fpAll->FixParameter(EXPFRAC, expFrac);
    fpAll->FixParameter(BKGCONST, bCnst);
    fpAll->FixParameter(BKGTAU, bTau);
    fpAll->FixParameter(AB, vabs[ichan]);
    fpAll->FixParameter(PPM, thePPM);
    fpAll->FixParameter(KQ, kq0);
    fpAll->FixParameter(KDIFFUSION, kdiffusion0);
  }

  /* limited parameters */
  // fpAll->SetParLimits(SFRAC, 0., 1.);
  // fpAll->SetParLimits(PPM, 0, 100. * thePPM);
  // fpAll->SetParLimits(AB, 0, 10. * vabs[ichan]);

  /*
  fpAll->SetParLimits(TAUM, .1 * tauM0, 100. * tauM0);
  fpAll->SetParLimits(KQ, .1 * kq0, 100. * kq0);
  fpAll->SetParLimits(KDIFFUSION, 0, 0.1);

  fpAll->SetParLimits(TAU3, 200., 2000.);
  fpAll->FixParameter(TAU3, tTriplet0);
  */

  // fpAll->FixParameter(KQ, kq0);
  // fpAll->FixParameter(KDIFFUSION, kdiffusion0); // this is starting value Eur. Phys. J. C (2013) 73:2618
  // fpAll->SetParameter(TAUM, tauM0);

  printf("\n\n  >>> modelFit initial value parameters fit chan %i ppm %.3f \n", (int)fpAll->GetParameter(CHAN), fpAll->GetParameter(PPM));
  printf(">>> effGeo %E norm %E  effGeo*norm %E sfrac %E \n", effChan, nPhotons, nPhotons * effChan, sfracStart);

  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  param %i %s %.3E +/- %.3E \n", ii, fpAll->GetParName(ii), fpAll->GetParameter(ii), fpAll->GetParError(ii));
  }
  cout << " ---------------  " << endl;

  cout << " ----------------------------------------------------------  " << endl;
  models[MODELALL]->show();
  cout << " ----------------------------------------------------------  " << endl;

  // PPPPPPP return;
  if (flag.Contains("P"))
    return;

  gStyle->cd();
  gStyle->SetOptFit(1111111);

  /*
  TCanvas *canShow = new TCanvas(Form("WaveStart-%.3f-PPM", ppm), Form("WaveStart%.3f-PPM", ppm));

  // fpAll->SetRange(startTime, endTime); // axis values
  canShow->SetLogy();
  hmodel[MODELALL]->Draw();
  hFullFit->Draw("sameHIST");
  */

  /*
  hmodel[MODELALL]->Draw("sameHIST");
  hmodel[MODELSINGLET]->Draw("sameHIST");
  hmodel[MODELRECOMB]->Draw("sameHIST");
  hmodel[MODELTRIPLET]->Draw("sameHIST");
  */
  // fpAll->SetLineColor(kMagenta);
  // fpAll->Draw("sames");
  //  PPPPPPP   return;

  // SSSSSSS show start return
  if (flag.Contains("S"))
    return;

  /* *****  do fit here */
  fpAll->GetRange(fitLow, fitHigh);
  printf("fit all in range  %f - %f \n", fitLow, fitHigh);
  fptr = hWave->Fit(fpAll, "WLRS0+", "", fitLow, fitHigh);
  // fptr = hWave->Fit(fpAll, "LE0SR+", "", fitStart, fitEnd);
  int fullResult = fptr;
  if (fullResult != 0)
  {
    printf("!!!!!! full fit fails %i  status = migradStatus + 10*minosStatus + 100*hesseStatus + 1000*improveStatus \n", fullResult);
    return;
  }

  cov = fptr->GetCorrelationMatrix();
  printf(" correlation chan %i cov(2,3) %f \n", ichan, cov(2, 3));
  fptr->Print("V");

  /* fill fitted histogram
  hFullFit->GetXaxis()->SetRangeUser(startTime, endTime);
  hFullFit->SetMarkerColor(kRed);
  hFullFit->SetLineColor(kRed);
  fout->Add(hSingletFit);
  for (int ibin = 0; ibin < hFullFit->GetNbinsX(); ++ibin)
  {
    double val = fpAll->Eval(hFullFit->GetBinCenter(ibin));
    if (val < tiny)
      val = 0;
    hFullFit->SetBinContent(ibin, val);
  }
  */
  // TPaveStats *st = (TPaveStats *)hWave->GetListOfFunctions()->FindObject("stats");
  // hWave->GetListOfFunctions()->Clear();

  fout->Append(hWave);
  fout->Append(fpAll);
  /* set parameters after fit */
  plotModels(MODELALL, ichan, hWave);

  gStyle->SetOptFit(1111111);
  gStyle->SetOptStat(0);
  TString canFitName;
  if (fixTriplet)
    canFitName.Form("WaveTripletFitChan%i-dopant-%.3f-PPM", ichan, ppm);
  else
    canFitName.Form("WaveFitChan%i-dopant-%.3f-PPM", ichan, ppm);
  TCanvas *canFit = new TCanvas(canFitName, canFitName);
  canFit->SetLogy();

  if (fixTriplet)
  {
    hWave->SetName(Form("WaveTriplet%.0fChan%i", tau3, ichan));
    hWave->SetTitle(Form("WaveTriplet%.0fChan%i", tau3, ichan));
  }
  else
    hWave->SetTitle(canFitName);

  hWave->GetYaxis()->SetRangeUser(1E-8, 2.E-1);
  hWave->GetXaxis()->SetRangeUser(startTime, 7500.);
  hWave->SetMarkerSize(.2);

  // gStyle->SetOptStat(); // for example
  hWave->Draw("p");
  hmodel[MODELALL]->Draw("sameHIST");
  hmodel[MODELTRIPLET]->Draw("sameHIST");
  hmodel[MODELSINGLET]->Draw("sameHIST");
  hmodel[MODELRECOMB]->Draw("sameHIST");
  //  fpAll->SetLineColor(kRed);
  //  fpAll->SetLineStyle(5);
  //  fpAll->SetLineWidth(4);
  canFit->Modified();
  canFit->SetGrid();
  gStyle->SetOptFit(1); // for example
  canFit->Update();
  canFit->BuildLegend();
  canFit->Print(".pdf");

  norm = fpAll->GetParameter(NORM);
  normErr = fpAll->GetParError(NORM);
  sfrac = fpAll->GetParameter(SFRAC);
  sfracErr = fpAll->GetParError(SFRAC);
  double arecomb = fpAll->GetParameter(ARECOMB);
  double trecomb = fpAll->GetParameter(TRECOMB);
  double recombPower = fpAll->GetParameter(RECOMBPOWER);
  double tripletTau = fpAll->GetParameter(TAU3);
  double mixedTau = fpAll->GetParameter(TAUM);
  double tripletTauErr = fpAll->GetParError(TAU3);
  double mixedTauErr = fpAll->GetParError(TAUM);
  double ab = fpAll->GetParError(AB);
  int itype = int(fpAll->GetParError(TYPE));

  dataElement.Form(" %.3E %.3E %.3E %.3E %.3E %.3E %.3E %.3E %i \n", norm, normErr, sfrac, sfracErr, tripletTau, tripletTauErr, mixedTau, mixedTauErr, ichan);
  cout << dataElement;

  fpAll->GetRange(fitLow, fitHigh);
  printf("fit All in range  %f - %f \n", fitLow, fitHigh);

  double transmission = 1. - ab;
  // double normalization = binwidth * sfrac * effGeoFunc(ichan) * fillFactor * transmission * SiPMQ128;
  double normalization = norm * effGeo * effAndFill * SiPMQ128;
  double alpha1 = normalization * sfrac;
  // double alpha1 = binwidth * norm * sfrac * effGeoFunc(ichan) * fillFactor * transmission * SiPMQ128;
  // double alpha3 = binwidth * norm * (1. - sfrac) * effGeoFunc(ichan) * fillFactor * transmission * SiPMQ128;
  double recombCut = recombCutParameter * trecomb;

  double rNorm = binwidth * effGeo * effAndFill * transmission * SiPMQE150;

  double alphaR = rNorm;
  double alpha3 = alpha1 * (1. - sfrac);
  double alphaSum = alpha1 + alpha3;
  printf(" final values for fit chan %i type %i norm %E (%E) recombCut %E sfrac %.4f alpha1 %E alpha3 %E alpha sum %E alphaR %E \n", ichan, itype, norm, normalization, recombCut, sfrac, alpha1, alpha3, alphaSum, alphaR);
  //
  fitFile << dataElement.Data();

  double singletSum = hmodel[MODELSINGLET]->Integral(startBin, endBin);
  double tripletSum = hmodel[MODELTRIPLET]->Integral(startBin, endBin);
  double recombSum = hmodel[MODELRECOMB]->Integral(startBin, endBin);
  double allSum = hmodel[MODELALL]->Integral(startBin, endBin);
  double singletSumFrac = singletSum / (singletSum + tripletSum);

  double singletSum0 = hmodel[MODELSINGLET]->Integral(startBin, singletEndBin);
  double tripletSum0 = hmodel[MODELTRIPLET]->Integral(startBin, singletEndBin);
  double recombSum0 = hmodel[MODELRECOMB]->Integral(startBin, singletEndBin);
  double allSum0 = hmodel[MODELALL]->Integral(startBin, singletEndBin);

  cout << hmodel[MODELSINGLET]->GetName() << " , " << hmodel[MODELTRIPLET]->GetName() << " , " << hmodel[MODELRECOMB]->GetName() << endl;
  printf(" startBin %i endBin %i \n singlet %.4f \n triplet %.4f \n recomb %.4f \n all %.4f \n sfrac %.4f \n",
         startBin, endBin, singletSum, tripletSum, recombSum, allSum, singletSumFrac);

  printf(" singlet region sums  %i %i : \n singlet %.4f \n triplet %.4f \n recomb %.4f \n all %.4f  \n", startBin, singletEndBin, singletSum0, tripletSum0, recombSum0, allSum0);

  printf("fitResult  = %i \n", fitResult);
  fitFile.close();

  fout->Write();

  // for(int ifit=0; ifit<NTYPES; ++ifit ) delete models[ifit];
} //***************************** end of routine ***********
