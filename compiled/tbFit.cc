#include <iostream>
#include <fstream>
#include "modelFit.hh"
// time is in microseconds
using namespace TMath;
TFile *fin;
TDirectory *runSumDir;

int nominalTrigger = 729;
double startTime = 700.; // hWave->GetBinLowEdge(maxBin) + hWave->GetBinWidth(maxBin) / 2.;
double singletStart = 700;
double singletEnd = 750;
// 3603795;

double
fitBack(TH1D *hist)
{
  double lowCut = 12000;
  double highCut = 15000;
  TF1 *gfit = NULL;
  double ave = 0;
  int lowBins = hist->FindBin(lowCut);
  int highBins = hist->FindBin(highCut);
  auto fitBack = new TF1("fitBack", "pol0", lowCut, highCut);
  fitBack->SetParameter(0,1E-1);
  fitBack->SetParLimits(0,1.E-9,1E2);
  
  hist->Fit("fitBack", "LF", " ", lowCut, highCut);
  gfit = (TF1 *)hist->GetListOfFunctions()->FindObject("fitBack");
  if (gfit)
  {
    ave = gfit->GetParameter(0);
  }
  else
    printf("P1 Fit to hist fails \n");
  double aveb = hist->Integral(lowCut,highCut) / double(highBins-lowBins);
  printf(" \t\t ave %E aveb %E \n\n", ave, aveb);
  return ave;
}

static double singletPeak(double *xx, double *par)
{
  double t = xx[0];
  double norm = par[0];
  double mean = par[1];
  double sigma = par[2];
  double arg = (t - mean) / sigma;
  double binwidth = 2.; // 2 ns
  double g = binwidth * norm / sqrt(TMath::TwoPi()) / sigma * TMath::Exp(-0.5 * arg * arg);
  return g;
}

bool openFile(TString fileName = "summary-01_25_2024-nfiles-10-dir-caenDataGoldSave-2024-05-20-14-17.root")
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
    return false;
  }
  runSumDir->ls();
  return true;
}

void tbFit(int ichan = 8)
{

  if (!openFile())
    return;
  TFile *fout = new TFile("tbFit", "recreate");
  TString hname;
  hname.Form("RunHitWaveChan%i", ichan);
  TH1D *hWave = NULL;
  runSumDir->GetObject(hname, hWave);
  if (!hWave)
    return;

  hWave->GetListOfFunctions()->Clear();

  cout << " got " << hWave->GetName() << endl;
  double ppm = 0.05;

  // peak singlet fit
  TH1D *hSinglet = (TH1D *)hWave->Clone("SingletWave");
  gStyle->SetOptStat();
  gStyle->SetOptFit(1111);
  auto fsinglet = new TF1("singlet", singletPeak, singletStart,singletEnd, 3);
  fsinglet->SetParameters(1,double(nominalTrigger), 10.);
  fsinglet->SetParNames("norm", "mean", "sigma");

  TCanvas *canSinglet = new TCanvas(Form("SingletFit-%.3f-PPM", ppm), Form("SingletFit-%.3f-PPM", ppm));
  hSinglet->Fit("singlet");
  hSinglet->GetXaxis()->SetRangeUser(700,800);
  hSinglet->Draw();
  double sumSinglet = hWave->Integral(hWave->FindBin(700), hWave->FindBin(800));

  TCanvas *canSingletFunc = new TCanvas(Form("SingletFitFunc-%.3f-PPM", ppm), Form("SingletFitFunc-%.3f-PPM", ppm));
  fsinglet->Draw();

  printf(" singlet  integral 700 to 800 = %.3E\n", sumSinglet);

  // fit background
  //TH1D *hBack = (TH1D *)hWave->Clone("BackWave");
  //double back = fitBack(hBack);
  //cout << " fitted back  " << back << endl;

  //hWave->GetXaxis()->SetRangeUser(0, 15000);
  double back = 0;
  double markerSize = 0.5;

  /***** fitting *****/
  hWave->GetListOfFunctions()->Clear();
  Double_t binwidth = hWave->GetBinWidth(1);
  int maxBin = hWave->GetMaximumBin();
  double hnorm = hWave->Integral(1, 7500);

  printf("startTime %f hnorm %f \n",startTime,hnorm);

  int theFit = 4;
  modelFit *model = new modelFit(theFit, ichan, ppm);
  TF1 *fp = model->fp;
  double sFrac = 0.2;

  int ilevel = level(ichan);
  double dist = distanceLevel[ilevel];
  double ab = Absorption(ppm, dist);
  double kplus = 1;
  double kPrime = 1;

  /*

  fp->SetParName(0, "norm");
  fp->SetParName(1, "PPM");
  fp->SetParName(2, "tau3");
  fp->SetParName(3, "kp");
  fp->SetParName(4, "sfrac");
  fp->SetParName(5, "rfrac"); // recombination this is starting value Eur. Phys. J. C (2013) 73:2618
  fp->SetParName(6, "ab");
  fp->SetParName(7, "kxprime");
  fp->SetParName(8, "tmix");
  fp->SetParName(9, "bgk");
  fp->SetParName(10, "chan");
  fp->SetParName(11, "binw");
  fp->SetParName(12, "type");
  */

  fp->SetParameter(0,nPhotons);
  // up to 8 parameters
  fp->FixParameter(1, ppm);
  fp->SetParLimits(2, 800.,2000. );
  //fp->FixParameter(2, tTriplet);
  fp->FixParameter(3, kplus);
  //fp->FixParameter(4, sFrac); // not fitting singlet
  fp->FixParameter(5, 1.E-5);
  fp->SetParLimits(5, 1.E-6 , 1.E-4);
  fp->FixParameter(6,ab);
  // fp->SetParLimits(6, 1.E-9, 1.);
  //fp->FixParameter(7, 2 * kPrime);
  fp->FixParameter(8, tMix);
  fp->FixParameter(9, back);
  //fp->FixParameter(12,0); // model fit

  printf("  >>> modelFit initial value parameters fit chan %i ppm %.3f \n", (int)fp->GetParameter(10), fp->GetParameter(1));
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  param %i %s %.3E +/- %.3E \n", ii, fp->GetParName(ii), fp->GetParameter(ii), fp->GetParError(ii));
  }
  cout << " ---------------  " << endl;

  gStyle->cd();
  gStyle->SetOptFit(1111111);
  TCanvas *canShow = new TCanvas(Form("WaveShow-%.3f-PPM", ppm), Form("WaveShow-%.3f-PPM", ppm));
  canShow->SetLogy();
  hWave->Draw();
  fp->Draw("sames");
  fp->SetRange(xMin, xMax);

  /* do the fit here */
  TFitResultPtr fptr = hWave->Fit(fp, "RLE0S+", "", 0, xMax);
  TMatrixDSym cov = fptr->GetCorrelationMatrix();
  printf(" correlation chan %i cov(2,3) %f \n", ichan, cov(2, 3));
  fptr->Print("V");

  gStyle->cd();
  gStyle->SetOptFit(1111111);
  TCanvas *canFit = new TCanvas(Form("WaveFit-%.3f-PPM", ppm), Form("WaveFit-%.3f-PPM", ppm));
  canFit->SetLogy();
  // hWave->GetYaxis()->SetRangeUser(1E-1,1E3);
  hWave->SetMarkerSize(0.2);
  fout->Append(hWave);
  fout->Append(fp);
  hWave->Draw("p");
  fp->SetLineColor(kRed);
  fp->SetLineStyle(5);
  fp->SetLineWidth(4);
  fp->Draw("same");
  TPaveStats *st = (TPaveStats *)hWave->GetListOfFunctions()->FindObject("stats");
  gStyle->SetOptFit(); // for example
  canFit->Modified();
  canFit->Update();
  hWave->GetListOfFunctions()->ls();
  double xlow, xhigh;
  fp->GetRange(xlow, xhigh);
  cout << " trigger Time " <<  nominalTrigger << " hist integral  " << hnorm << " fit from " << xlow << " to " << xhigh << endl;
  printf(" singlet  integral 1360 to 1460 = %.3E\n", sumSinglet);
  model->show();
  //model->showEff();
}
