#include <iostream>
#include <fstream>
#include "modelFit.hh"
// time is in microseconds
using namespace TMath;
TFile *fin;
TDirectory *runSumDir;

bool openFile(TString fileName = "summary-type-1-dir-caenData-2023-06-28-15-09.root")
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
  double hnorm = hWave->Integral(1,maxBin);
  cout << startTime << " hnorm " << hnorm << endl;

  int ifit = 4;
  double ppm = 0.05;
  int theFit=4;
  modelFit *model = new modelFit(theFit,ichan,ppm);
  TF1 *fp = model->fp;
  double sFrac = 0.2;

  int ilevel = level(ichan);
  double dist = getDistance(ilevel);
  double ab = Absorption(ppm, dist);
  double kplus = 1;
  double kPrime = 1;



  fp->FixParameter(0, binwidth);
  fp->SetParameter(1, hnorm);
  fp->FixParameter(2, ppm);
  fp->FixParameter(3, 1300);
  fp->FixParameter(4, kplus);
  fp->SetParameter(5, sFrac);
  //fp->SetParLimits(5, .01, .5);
  fp->FixParameter(6, ab);
  //fp->SetParLimits(6, 1.E-9, 1.);
  fp->FixParameter(7, 2 * kPrime);
  fp->FixParameter(11,10);

  
  cout << " initial values " << endl;
   printf(" \n\n >>> modelFit fitted parameters fit chan %i ppm %f \n", (int)fp->GetParameter(11), fp->GetParameter(2));
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  param %i %s %.3E +/- %.3E \n", ii, fp->GetParName(ii), fp->GetParameter(ii), fp->GetParError(ii));
  }

  cout << " ---------------  " << endl;

  /* do the fit here */
  double xstart = 1400.;
  double xstop = 75000;
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
