#include <iostream>
#include <fstream>
#include "modelAFit.hh"
// time is in microseconds
using namespace TMath;
TFile *fin;
TDirectory *runSumDir;
vector<TH1D *> vhist;

double fitBack(TH1D *hist)
{
  double lowCut = 12000;
  double highCut = 15000;
  TF1 *gfit = NULL;
  double ave = 0;
  int lowBins = hist->FindBin(lowCut);
  int highBins = hist->FindBin(highCut);
  auto fitBack = new TF1("fitBack", "pol0", lowCut, highCut);
  fitBack->SetParameter(0, 1E-1);
  fitBack->SetParLimits(0, 1.E-9, 1E2);

  hist->Fit("fitBack", "LF", " ", lowCut, highCut);
  gfit = (TF1 *)hist->GetListOfFunctions()->FindObject("fitBack");
  if (gfit)
  {
    ave = gfit->GetParameter(0);
  }
  else
    printf("P1 Fit to hist fails \n");
  double aveb = hist->Integral(lowCut, highCut) / double(highBins - lowBins);
  printf(" \t\t ave %E aveb %E \n\n", ave, aveb);
  return ave;
}

bool openFile(TString fileName = "summary-type-1-dir-caenData-2023-07-13-13-37.root")
{
  // open input file
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
  for (int ic = 0; ic < 13; ++ic)
  {
    if (ic == 5 || ic == 6 || ic == 3 || ic > 8)
      continue;
    TString hname;
    hname.Form("RunSumHitWaveChan%i", ic);
    TH1D *hWave = NULL;
    runSumDir->GetObject(hname, hWave);
    if (hWave != NULL)
      vhist.push_back(hWave);
  }
  return true;
}

void tbFitAll()
{

  TFile *fout = new TFile("tbFitAll", "recreate");

  if (!openFile())
    return;
  // fill buff
  cout << " got " << vhist.size() << endl;
  for (unsigned ih = 0; ih < vhist.size(); ++ih)
  {
    for (int ib = 1; ib < vhist[ih]->GetNbinsX(); ++ib)
      buff[ih][ib] = vhist[ih]->GetBinContent(ib);
  }

  /* parameter definitions
  fp->SetParName(0, "norm");
  fp->SetParName(1,"ppm");
  fp->SetParName(2, "tau3");
  fp->SetParName(3, "kp");
  fp->SetParName(4, "sfrac");
  fp->SetParName(5, "rfrac");
  fp->SetParName(6, "kxprime");
  fp->SetParName(7, "tmix");
  fp->SetParName(8, "bkg");
  */

  // Set starting values and step sizes for parameters
  static Double_t vstart[NPARS];
  static Double_t step[NPARS];
  TString parNames[NPARS];
  parNames[0] = TString("norm");
  parNames[1] = TString("ppm");
  parNames[2] = TString("tau3");
  parNames[3] = TString("kp");
  parNames[4] = TString("sfrac");
  parNames[5] = TString("rfrac");
  parNames[6] = TString("kxprime");
  parNames[7] = TString("tmix");
  parNames[8] = TString("bkg");
  vstart[0] = 4.3344E+04;
  vstart[1] = 0.05;
  vstart[2] = 1.1777E+03;
  vstart[3] = kplus;
  vstart[4] = 0.886;
  vstart[5] = 0.001;
  vstart[6]=6.2;
  vstart[7] =4.7E3;
  vstart[8] = 2.4E-6;

  TMinuit *gMinuit = new TMinuit(NPARS); // initialize TMinuit with a maximum of 5 params
  gMinuit->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 0.5; // for likelihood
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  for (unsigned j = 0; j < NPARS; ++j)
  {
    step[j] = 1.E-2 * vstart[j];
    gMinuit->mnparm(j, parNames[j].Data(), vstart[j], step[j], 0, 0, ierflg);
  }
  double currentValue;
  double currentError;

  double fval;
  double gin[NPARS];
  int npar = NPARS;
  fcn(npar, gin, fval, vstart, 0);
  printf(" >>>>   fval %E \n", fval);

  // Now ready for minimization step
  arglist[0] = 1000; // maxcalls
  arglist[1] = 1.E-2;   // tolerance

  /* MIGrad[maxcalls][tolerance]*/
  //gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // Print results
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  /*
  Prints the values of the parameters at the time of the call.
  According to the value of IKODE, the printout is: IKODE=INKODE=

0 only info about function value
1 parameter values, errors, limits
2 values, errors, step sizes, internal values
3 values, errors, step sizes, first derivs.
4 values, parabolic errors, MINOS errors
when INKODE=5, MNPRIN chooses IKODE=1,2, or 3, according to fISW[1]
  */
  gMinuit->mnprin(4, amin);

  // put parameters into model par array
  for (int k = 0; k < NPARS; ++k)
  {
    gMinuit->GetParameter(k, currentValue, currentError);
    lpar[k] = currentValue;
  }
  int ichan[6] = {0, 1, 2, 4, 7, 8};
  TF1 *ffit[6];
  // create functions to plot
  double xlow = 0;
  double xhigh = 15000;
  gStyle->SetOptFit(1111111);
  for (int ifit = 0; ifit < 6; ++ifit)
  {
    TString cname;
    cname.Form("FitChan-%i-Dopant-%.f", ichan[ifit], lpar[1]);
    TCanvas *can = new TCanvas(cname, cname);
    can->SetLogy();
    ffit[ifit] = new TF1(Form("ModelFitChan-%i", ichan[ifit]), modelFunc, xlow, xhigh, 1);
    ffit[ifit]->SetParameter(0, ichan[ifit]);
    vhist[ifit]->GetListOfFunctions()->Clear();
    vhist[ifit]->Draw();
    ffit[ifit]->Draw("same");
    fout->Add(ffit[ifit]);
    fout->Add(vhist[ifit]);
    can->Print(".png");
    show();
  }
}
