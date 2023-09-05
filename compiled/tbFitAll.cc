#include <iostream>
#include <fstream>
#include "modelAFit.hh"
// time is in microseconds
using namespace TMath;
TFile *fin;
TDirectory *runSumDir;
std::vector<TH1D *> vhist;
TString histSet;
double dopant[10];

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

int openFile(int fileNum=0)
{
  TString summaryFile[10];
  summaryFile[0] = TString("summary-type-1-dir-caenDataZeroPPM-2023-08-27-13-39.root"); // 05_22_2023 05_30_2023
  summaryFile[1] = TString("summary-type-1-dir-caenDataPointOnePPM-2023-08-27-13-30.root"); //06_01_2023 06_02_2023 
  summaryFile[2] = TString("summary-type-1-dir-caenDataPointTwoPPM-2023-08-27-13-31.root");//07_06_2023  07_07_2023
  dopant[0] = 0.05;
  dopant[1] = 0.1;
  dopant[2] = 0.2;
  TString fileName = summaryFile[fileNum];

  // open input filex
  vhist.resize(13);
  for (unsigned k = 0; k < vhist.size(); ++k)
    vhist[k] = NULL;
  printf(" looking for summary file %s\n", fileName.Data());

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
  // file exists, so open with TFile
  fin = new TFile(fileName, "readonly");
  printf(" opened file %s\n", fileName.Data());
  runSumDir = NULL;
  // get subdirectory pointer
  fin->GetObject("runSumDir", runSumDir);
  if (!runSumDir)
  {
    printf(" no runSumDir in file %s\n", fileName.Data());
    return false;
  }
  runSumDir->ls();
  cout << " making histogram list with " << histSet << "  " << runSumDir->GetName()  << endl;
  int iGot=0;
  for (int ic = 0; ic < 13; ++ic)
  {
    if (!goodChannel(ic))
      continue;
    TString hname;
    hname.Form("RunSumWaveChan%i", ic);
    TH1D *hWave = NULL;
    printf("\t look for  chan %i %s \n",ic,hname.Data());
    // get histogram by name
    runSumDir->GetObject(hname, hWave);
    if (hWave != NULL)
    {
      vhist[ic] = hWave;
      printf("\t got chan %i %s \n",ic,hWave->GetName());
      background[ic] = fitBack(hWave);
      ++iGot;
    }
  }
  return iGot;
}

void tbFitAll(int fileNum=0)
{

  TFile *fout = new TFile("tbFitAll.root", "recreate");
  // get data histograms from file
  int iGot = openFile(fileNum);
  printf(" i fot %i \n",iGot);
  if (iGot<1)
    return;
  // fill buff

  // shift PMT by 46 ns = 23 bins
  int binShift = 23;
  for (int ibin = 0; ibin < vhist[12]->GetNbinsX() - binShift; ++ibin)
  {
    vhist[12]->SetBinContent(ibin,vhist[12]->GetBinContent(ibin+binShift)) ;
    vhist[12]->SetBinError(ibin,vhist[12]->GetBinError(ibin+binShift)) ;
  }

  for (unsigned ih = 0; ih < vhist.size(); ++ih)
    {
      if (vhist[ih] == NULL)
        continue;
      cout << ".... for channel %i " << ih << "  " << vhist[ih]->GetName() << " maximum bin " << vhist[ih]->GetMaximumBin() << endl;
      fout->Add(vhist[ih]);
      // fill data buffer
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
  vstart[0] = 2.06322e+03;
  vstart[1] = dopant[fileNum];
  vstart[2] = 1.1777E+03;
  vstart[3] = kplus; // defined in  modelAFit.hh
  vstart[4] = 0.2;  //0.886;
  vstart[5] = 1.E-3; // recomb  1.E-3;
  vstart[6] = 2.0;
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

  double fval = 0;
  double gin[NPARS];
  int npar = NPARS;
  fcn(npar, gin, fval, vstart, 0);
  printf(" >>>>   fval %E \n", fval);

  double currentValue;
  double currentError;
  

  // fix ppm
  //gMinuit->FixParameter(1);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = 5; // par sfrac
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = 2; // par PPM
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = 1;  // par kp
  arglist[1] = 100.; // low
  arglist[2] = 1.E10; // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = 4;  // par kp
  arglist[1] = 0.; // low
  arglist[2] = 10.; // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

  arglist[0] = 6;  // par rfrac
  arglist[1] = 0.; // low
  arglist[2] = 1.E-2; // high
  gMinuit->mnexcm("SET LIM", arglist, 3 , ierflg);

  arglist[0] = 7;  // kxprime
  arglist[1] = 0.; // low
  arglist[2] = 10.E-2; // high
  gMinuit->mnexcm("SET LIM", arglist, 3 , ierflg);



  printf(" \n\n >>> modelFit start parameters fit ppm %f \n", vstart[1]);
  for (int ii = 0; ii < NPARS; ++ii)
  {
    gMinuit->GetParameter(ii, currentValue, currentError);
    printf("\t  param %i %s %.4E  \n", ii, lparNames[ii].Data(), currentValue);
  }

  // gMinuit->FixParameter(6);
  // create functions to plot
  int myColor[13] = {kRed, kBlue - 9, kGreen, kOrange, kBlue, kAzure + 3, kCyan - 3, kGreen - 4, kGreen, kSpring, kYellow, kBlue + 9, kTeal - 4};
  int myStyle[13] = {21, 22, 23, 24, 25, 26, 21, 22, 23, 31, 32, 33, 34};
  TString cname;

  for (int ifit = 0; ifit < 13; ++ifit)
  {
    if (vhist[ifit] == NULL)
      continue;
    vhist[ifit]->GetListOfFunctions()->Clear();
    vhist[ifit]->SetLineColor(myColor[ifit]);
    vhist[ifit]->SetMarkerColor(myColor[ifit]);
    vhist[ifit]->SetMarkerStyle(myStyle[ifit]);
    //if (histSet == TString("HitWave"))
      //vhist[ifit]->GetYaxis()->SetRangeUser(1.E-7, 1.);
    //else
    //vhist[ifit]->GetYaxis()->SetRangeUser(1.E-3, 1.E3);
  }

  cname.Form("DataDopant-%.f", lpar[1]);
  TCanvas *cand = new TCanvas(cname, cname);
  cand->SetLogy();
  for (int k = 0; k < 13; ++k)
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
  createFunctions();
  //show();

  cname.Form("FitChan-Dopant-%.3f", lpar[1]);
  TCanvas *canAllChan = new TCanvas(cname, cname);
  canAllChan->SetLogy();
  for (int k = 0; k<13; ++k)
  {
    if (!goodChannel(k))
      continue;
    ffit[k]->SetLineColor(myColor[k]);
    //ffit[k]->GetYaxis()->SetRangeUser(1.E-10, 1.);
    if(k==0) ffit[k]->Draw();
    else
      ffit[k]->Draw("sames");
  }
  canAllChan->BuildLegend();
  canAllChan->Print(".png");

  // PMT components
  cname.Form("ModelPmt-Dopant-%.3f-recombination-%.2E", lpar[1], vstart[5]);
  TCanvas *canPmt = new TCanvas(cname, cname);
  canPmt->SetLogy();
  for (int k = 0; k < 6; ++k)
  {
    ffitPmt[k]->SetLineColor(myColor[k]);
    //ffitPmt[k]->GetYaxis()->SetRangeUser(1.E-7, 3E-3);
    fout->Add(ffitPmt[k]);
    if (k == 0)
      ffitPmt[k]->Draw();
    else
      ffitPmt[k]->Draw("sames");
  }
  canPmt->BuildLegend();
  canPmt->Print(".png");

  // channel7  components
  cname.Form("ModelChan7Dopant-%.3f", lpar[1]);
  TCanvas *canChan = new TCanvas(cname, cname);
  canChan->SetLogy();
  for (int k = 0; k < 6; ++k)
  {
    //ffitChan[k]->GetYaxis()->SetRangeUser(1.E-10, 3E-1);
    fout->Add(ffitChan[k]);
    ffitChan[k]->Print();
    ffitChan[k]->SetLineColor(myColor[k]);
    ffitChan[k]->GetYaxis()->SetRangeUser(1.E-4, 1.);
    printf(" comp %i\n", int(ffitChan[k]->GetParameter(1)));
    if (k == 0)
      ffitChan[k] ->Draw();
    else
      ffitChan[k]->Draw("sames");
  }
  canChan->BuildLegend();
  canChan->Print(".png");


  /* if just plotting model return here*/
  //return;

  // minimize with MIGRAD
  // Now ready for minimization step
  arglist[0] = 100000; // maxcalls
  arglist[1] = 1.E-2;  // tolerance

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
    printf("\t copy %s new value %f \n", lparNames[k].Data(), lpar[k]);
  }


  cname.Form("ModelAllChannelsDopant-%.f", lpar[1]);
  TCanvas *canf = new TCanvas(cname, cname);
  canf->SetLogy();
  for (int k = 0; k < 13; ++k)
  {
    if (!goodChannel(k))
      continue;
    fmodel[k]->SetLineColor(myColor[k]);
    fmodel[k]->GetYaxis()->SetRangeUser(1.E-4,1.);
    if (k == 0)
      fmodel[k]->Draw();
    else
      fmodel[k]->Draw("sames");
  }
  canf->BuildLegend();
  canf->Print(".png");

  gStyle->SetOptFit(1111111);
  for (int k = 0; k < 13; ++k)
  {
    if (!goodChannel(k))
      continue;
    cname.Form("FitToChan%i-Dopant-%.3f-recombination-%.2E", k,lpar[1], vstart[5]);
    TCanvas *can = new TCanvas(cname, cname);
    can->SetLogy();
    ffit[k]->SetLineColor(kBlack);
    ffit[k]->SetLineWidth(4);
    vhist[k]->SetMarkerSize(0.4);
    //vhist[k]->GetYaxis()->SetRangeUser(1E-4,1.);
    vhist[k]->Draw("");
    ffit[k]->Draw("same");
    fout->Add(ffit[k]);
    can->Print(".png");
  }
  fout->Write();
}
