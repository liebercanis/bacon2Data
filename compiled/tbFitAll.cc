#include <iostream>
#include <fstream>
#include "modelAFit.hh"
// time is in microseconds
using namespace TMath;
TFile *fin;
TDirectory *runSumDir;
std::vector<TH1D *> vhist;

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
  vhist.resize(13);
  for (unsigned k = 0; k < vhist.size(); ++k)
    vhist[k] = NULL;
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
    if (!goodChannel(ic))
      continue;
    TString hname;
    hname.Form("RunSumHitWaveChan%i", ic);
    TH1D *hWave = NULL;
    runSumDir->GetObject(hname, hWave);
    if (hWave != NULL)
    {
      vhist[ic]=hWave;
      background[ic] = fitBack(hWave);
    }
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
    if(vhist[ih]==NULL)
      continue;
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
  */

  // Set starting values and step sizes for parameters
  static Double_t vstart[NPARS];
  static Double_t step[NPARS];
  // fit starting values  
  vstart[0] = 4.3344E+04;
  vstart[1] = 0.05;
  vstart[2] = 1.1777E+03;
  vstart[3] = kplus;
  vstart[4] = 0.886;
  vstart[5] = 0.001;
  vstart[6] = 6.2;
  vstart[7] = 4.7E3;

  TMinuit *gMinuit = new TMinuit(NPARS); // initialize TMinuit with a maximum of 5 params
  gMinuit->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 0.5; // for likelihood
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  setParNames();
  for (unsigned j = 0; j < NPARS; ++j)
  {
    step[j] = 1.E-6 * vstart[j];
    gMinuit->mnparm(j, lparNames[j].Data(), vstart[j], step[j], 0, 0, ierflg);
    lpar[j] = vstart[j];
  }


  double fval=0;
  double gin[NPARS];
  int npar = NPARS;
  fcn(npar, gin, fval, vstart, 0);
  printf(" >>>>   fval %E \n", fval);

  double currentValue;
  double currentError;
  printf(" \n\n >>> modelFit start parameters fit ppm %f \n", vstart[1]);
  for (int ii = 0; ii < NPARS; ++ii)
  {
    gMinuit->GetParameter(ii, currentValue, currentError);
    printf("\t  param %i %s %.4E  \n", ii, lparNames[ii].Data(),currentValue);
  }

  // fix ppm
  gMinuit->FixParameter(1);

  // minimize with MIGRAD
  // Now ready for minimization step
  arglist[0] = 100000;  // maxcalls
  arglist[1] = 1.E-2; // tolerance

  /* MIGrad[maxcalls][tolerance]*/
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

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
  // create functions to plot
  int myColor[13] = {kRed, kBlue-9, kGreen, kTeal-4,kOrange, kBlue, kAzure+3, kCyan-3, kGreen-4, kGreen, kSpring, kYellow, kOrange};
  int myStyle[13] = {21, 22, 23, 24, 25, 26, 21, 22, 23, 31, 32, 33, 34};

  for (int ifit = 0; ifit < 13; ++ifit)
  {
    if(vhist[ifit]==NULL)
      continue;
    vhist[ifit]->SetLineColor(kWhite);
    vhist[ifit]->SetMarkerColor(myColor[ifit]);
    vhist[ifit]->SetMarkerStyle(myStyle[ifit]);
    vhist[ifit]->GetYaxis()->SetRangeUser(1.E-7,1.);
  }

  createFunctions();

  TString cname;
  gStyle->SetOptFit(1111111);
  for (int k = 0; k < 13; ++k )
  {
    if (!goodChannel(k))
      continue;
    cname.Form("FitChan-%i-Dopant-%.3f",k, lpar[1]);
    TCanvas *can = new TCanvas(cname, cname);
    can->SetLogy();
    ffit[k]->SetLineColor(kBlack);
    vhist[k]->GetListOfFunctions()->Clear();
    vhist[k]->Draw();
    ffit[k]->Draw("same");
    fout->Add(ffit[k]);
    fout->Add(vhist[k]);
    can->Print(".png");
  }

  show();


  cname.Form("DataDopant-%.f", lpar[1]);
  TCanvas *cand = new TCanvas(cname, cname);
  cand->SetLogy();
  for (int  k = 0; k < 13; ++k )
  {
    if (!goodChannel(k))
      continue;
    if (k == 0)
      vhist[k]->Draw();
    else
      vhist[k]->Draw("sames");
  }
  cand->BuildLegend();
  cand->Print(".png");

  cname.Form("FitDopant-%.f", lpar[1]);
  TCanvas *canf = new TCanvas(cname, cname);
  canf->SetLogy();
  for (int k = 0; k < 13; ++k)
  {
    if (!goodChannel(k))
      continue;
    //ffit[ifit]->SetLineColor(myColor[ifit]);
    ffit[k]->GetYaxis()->SetRangeUser(1.E-7, 1.);
    if (k == 0)
      ffit[k]->Draw();
    else
      ffit[k]->Draw("sames");
  }
  canf->BuildLegend();
  canf->Print(".png");
}
