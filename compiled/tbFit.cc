#include <iostream>
#include <fstream>
#include "modelFit.hh"
// time is in microseconds
using namespace TMath;
TFile *fin;
TDirectory *runSumDir;

bool openFile(TString fileName = "summary-type-1-dir-caenData-2023-06-20-14-09.root")
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


void tbFit(int ichan =8)
{ 

  if(!openFile())
    return;
  TString hname;
  hname.Form("RunSumHitWaveChan%i",ichan);
  TH1D *hWave=NULL;
  runSumDir->GetObject(hname,hWave);
  if(!hWave)
    return;

  cout <<" got " << hWave->GetName() << endl;

  double markerSize = 0.5;
  /***** fitting *****/
  Double_t binwidth = hWave->GetBinWidth(1);
  int maxBin = hWave->GetMaximumBin();
  double startTime = hWave->GetBinLowEdge(maxBin) + hWave->GetBinWidth(maxBin) / 2.;
  cout << startTime << endl;

  int ifit = 4;
  double ppm = 0.07;
  modelFit *model = new modelFit(ifit, 0);
  TF1 *fp = model->fp;
  double norm = 1.138E6;
  double sFrac = 0.2;
  fp->FixParameter(1, norm);

  fp->FixParameter(0, binwidth);
  fp->SetParameter(1, norm);
  fp->FixParameter(2, ppm);
  fp->SetParameter(3, tTriplet);
  fp->FixParameter(5, sFrac);
  fp->FixParameter(8, ifit);
  fp->FixParameter(9, tMix);
  fp->Print();

  /* do the fit here */
  double xstart = 100.;
  double xstop = 14000;
  TFitResultPtr fptr = hWave->Fit(fp, "LE0S+", "", xstart, xstop);
  TMatrixDSym cov = fptr->GetCorrelationMatrix();
  printf(" correlation chan %i cov(2,3) %f \n", ichan, cov(2, 3));
  fptr->Print("V");

  

  gStyle->cd();
  gStyle->SetOptFit(1111111);
  TCanvas *canFit = new TCanvas(Form("WaveFit-%.f-PPM",0.07), Form("WaveFit-%.f-PPM",0.07));
  canFit->SetLogy();
  // hWave->GetYaxis()->SetRangeUser(1E-1,1E3);
  hWave->SetMarkerSize(0.5);
  hWave->Draw("p");
  fp->SetLineColor(kTeal - 6);
  fp->SetLineStyle(5);
  fp->SetLineWidth(4);
  fp->Draw("sames");
  TPaveStats *st = (TPaveStats *)hWave->GetListOfFunctions()->FindObject("stats");
  gStyle->SetOptFit(); // for example
  canFit->Modified();
  canFit->Update();

  
}
