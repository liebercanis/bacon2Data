/*********
 * fit one or all curves with theFitChannel = -1 this routine superceeds tbFitAll.cc which is deprecated
 modified to use modelAllFit.hh
 M.Gold Feb 17 2026 *
 ******************/
#include <iostream>
#include <fstream>
#include "TGraph.h"
#include "TMinuit.h"
#include "modelAllFit.hh"
// time is in microseconds
using namespace TMath;
TFile *fin;
TFile *fout;
TNtuple *ntParScan;
std::vector<TH1D *> hnorm;
std::vector<TH1D *> hcurve;
std::vector<TH1D *> hffit;
std::vector<TH1D *> hmodel;
std::vector<TH1D *> hfitModel;

int nominalTrigger = 729;
double singletStart = 700;
double singletEnd = 750;
// 3603795;
static Double_t vstart[NPARS];
static Double_t step[NPARS];

int colors[NCHAN] = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink, kGray};

TCanvas *makeCanFit(int i1, int i2, TString canName)
{
  printf(" makeCanFit %s from %i to %i size %lu \n", canName.Data(), i1, 12, hnorm.size());
  bool firstPlot = true;
  TCanvas *can = new TCanvas(canName, canName);
  can->Divide(2, 2);
  int ipanel = 0;
  for (int i = i1; i <= i2; ++i)
  {
    ++ipanel;
    if (isBadChannel(i))
      continue;
    hnorm[i]->GetXaxis()->SetRangeUser(1000, 4000);
    hfitModel[i]->GetXaxis()->SetRangeUser(1000, 4000);
    hfitModel[i]->SetLineWidth(2);
    can->cd(ipanel);
    gPad->SetLogy();
    hnorm[i]->Draw("");
    hfitModel[i]->Draw("HISTSAME");
  }
  // can->BuildLegend();
  can->SetLogy();
  can->Print(".pdf");
  return can;
}

/* fill background array lateBkg in modelAllFit.hh */
void fillLateBkg()
{

  for (unsigned ih = 0; ih < hnorm.size(); ++ih)
  {
    // integrate bin ranges
    lateBkg[ih] = hnorm[ih]->Integral(6000, 7500) / double(1500);
    printf("%u %s late int inte %E per bin \n", ih, hnorm[ih]->GetName(), lateBkg[ih]);
  }
}

TGraph *myScan(int thePar, double xlow, double xhigh)
{
  int maxPoints = 100;
  std::vector<double> xval;
  std::vector<double> yval;

  // get min parameters
  double fitVal[NPARS];
  double fitErr[NPARS];
  for (int ipar = 0; ipar < NPARS; ++ipar)
  {
    gMinuit->GetParameter(ipar, fitVal[ipar], fitErr[ipar]);
    // printf("line141 ipar %i par %f err %f \n", ipar, fitVal[ipar], fitErr[ipar]);
  }

  double *fGin;
  double nLL;
  for (int i = 0; i < maxPoints; ++i)
  {
    double x = xlow + double(i) * (xhigh - xlow) / double(maxPoints);
    fitVal[thePar] = x;
    gMinuit->Eval(thePar, fGin, nLL, &fitVal[0], 4);
    xval.push_back(x);
    yval.push_back(nLL);
    // printf("line53 i %i par nph %f r %f theta %f  phi %f \n", i, fitVal[0], fitVal[1], fitVal[2], fitVal[3]);
    ntScan->Fill(nLL, fitVal[1], fitVal[2], fitVal[3]);
    // printf("mySCAN par %i x= %f nLL %E \n", i, x, nLL);
  }
  // make and return graph
  return new TGraph(maxPoints, &xval[0], &yval[0]);
}

/* stored in modelAllFit..
 static double fitWave[NCHAN][MAXSAMPLE];
static double fitComp[NCHAN][NUMCOMP][MAXSAMPLE];
*/

void fillFitWave(int ichan, TH1D *hist)
{
  std::cout << " fillFitWave " << ichan << "  " << hist->GetName() << std::endl;
  hist->Reset("ICES");
  for (int ib = 1; ib < hist->GetNbinsX(); ++ib)
  {
    double val = max(fitWave[ichan][ib], 1.E-9);
    hist->SetBinContent(ib, val);
    hist->SetBinError(ib, 0);
    hist->GetYaxis()->SetTitle("yield");
    hist->GetXaxis()->SetTitle("time [ns]");
    hffit[ichan] = hist;
  }
}
void fillCompWave(int ichan, int icomp, TH1D *hist)
{
  // std::cout << " fillFitWave " << ichan << "  " << hist->GetName() << std::endl;
  hist->Reset("ICES");
  for (int ib = 1; ib < hist->GetNbinsX(); ++ib)
  {
    double val = max(fitComp[ichan][icomp][ib], 1.E-9);
    // if (ichan == 8 && ib == 1500)
    //   printf("!!!! chan %i sample %i val %E \n", ichan, ib, val);
    hist->SetBinContent(ib, val);
    hist->SetBinError(ib, 0);
    hist->GetYaxis()->SetTitle("yield");
    hist->GetXaxis()->SetTitle("time [ns]");
  }
  std::cout << std::endl;
}

void getCurves()
{

  TIter next(fin->GetListOfKeys());
  TKey *key;
  int ifile = 0;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());

    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();

    if (TString(h->GetName()).Contains("CurveChan"))
      hcurve.push_back(h);

    if (TString(h->GetName()).Contains("NormChan"))
    {
      hnorm.push_back(h);
      fout->Append(h);
    }
  }
}

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

bool openFile(TString fileName)
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
  return true;
}

/* fit to channel -1 = ALL */
void tbFit(int theFitChannel = -1)
{

  hffit.resize(NCHAN);
  hmodel.resize(NCHAN);

  TString inputFile = TString("post-anaCRun-btbSimNEW-2026-02-13-100000-7857.root");
  inputFile = TString("post-11_19_2025-11_19_2025-1371746.root");

  if (!openFile(inputFile))
    return;

  double dopant = 1.E-2; // PPM

  /*
    do all setups here
  */
  setupModelAllFit();

  // is this simulation?
  bool isSim = false;
  if (inputFile.Contains("btb"))
  {
    isSim = true;
  }
  if (isSim)
    fout = new TFile(Form("tbFitSimPPM%.2f.root", dopant), "recreate");
  else
    fout = new TFile(Form("tbFitPPM%.2f.root", dopant), "recreate");

  if (geoVersionOld)
    printf("OLD level distances 0 = %.3f 1= %.3f 2= %.3f 3 %.3f 4 %.3f \n", distanceLevel[0], distanceLevel[1], distanceLevel[2], distanceLevel[3], distanceLevel[4]);
  else
    printf("NEW level distances 0 = %.3f 1= %.3f 2= %.3f 3 %.3f 4 %.3f \n", distanceLevel[0], distanceLevel[1], distanceLevel[2], distanceLevel[3], distanceLevel[4]);

  /* read in all needed histogrms */
  getCurves();

  for (unsigned ih = 0; ih < hnorm.size(); ++ih)
  {
    printf("%u %s inte %E \n", ih, hnorm[ih]->GetName(), hnorm[ih]->Integral());
    hnorm[ih]->GetListOfFunctions()->Clear();
  }

  /* add constant late average background */
  fillLateBkg();

  /* fill buffer */
  printf("fill buff \n");
  for (unsigned ichan = 0; ichan < NCHANPMT; ++ichan)
  {
    printf(".... fill buffer for channel %i  hist %s maximum bin %i \n", ichan, hnorm[ichan]->GetName(), hnorm[ichan]->GetMaximumBin());
    // fill data buffer
    for (int isample = 1; isample < MAXSAMPLE; ++isample)
    {
      buff[ichan][isample - 1] = hnorm[ichan]->GetBinContent(isample); // C starts array from zero
    }
  }

  /* setup fit */
  TMinuit *gMinuit = new TMinuit(NPARS); // initialize TMinuit with a maximum of 5 params
  gMinuit->SetFCN(fcn);

  double currentValue;
  double currentError;

  /* total photons per event LY defined in modelAllFit.hh */
  double startNorm = 60. * LY;
  Double_t arglist[10];
  int ierflg = 0;
  arglist[0] = 0.5; // for likelihood up from minimum for 1 sigma errors
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  /**  fit starting values **/
  vstart[NORM] = startNorm;
  vstart[TRIGSTART] = 2. * hnorm[9]->GetMaximumBin();
  vstart[SFRAC] = 0.14; //;0.23;   // Segretto PHYSICAL REVIEW D 103, 043001 (2021)
  vstart[PPM] = dopant;
  vstart[TAU3] = tTriplet0;
  vstart[TAUM] = 4700.0;
  vstart[BKGCONST] = 4.0E-6;
  vstart[BKGTAU] = 5000.;
  vstart[THECHANNEL] = theFitChannel;
  /*
  printf("starting parameter values \n");
  for (int ip = 0; ip < NPARS; ++ip)
    printf(" par %i %s start val %.3f \n", ip, lparNames[ip].Data(), vstart[ip]);
    */

  // copy into Minuit
  /* have to put some errors here otherwise it will be constant*/
  for (unsigned j = 0; j < NPARS; ++j)
  {
    step[j] = 1.E-6 * vstart[j];
    gMinuit->mnparm(j, lparNames[j].Data(), vstart[j], step[j], 0.1 * vstart[j], 10. * vstart[j], ierflg);
    lpar[j] = vstart[j];
  }

  /******************************/
  // fix parameters minuit is fortran hence +1!
  /******************************/

  arglist[0] = TRIGSTART + 1; // trigger
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = TAUM + 1; // par tau3
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = BKGCONST + 1; // par
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  arglist[0] = BKGTAU + 1; // par
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  // arglist[0] = SFRAC + 1; // kp
  // gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  /******************************
      set limits ... here par starts with 1 so add 1
   *****************************/
  arglist[0] = NORM + 1;         // par
  arglist[1] = 0.01 * startNorm; // low
  arglist[2] = 10. * startNorm;  // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  // gMinuit->mnexcm("FIX", arglist, 3, ierflg);

  arglist[0] = SFRAC + 1; // par
  arglist[1] = 0.01;      // low
  arglist[2] = 1.0;
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  // gMinuit->mnexcm("FIX", arglist, 3, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = TAU3 + 1;         // par
  arglist[1] = 0.01 * tTriplet0; // low
  arglist[2] = 2.0 * tTriplet0;  // high
  gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

  // set limits ... here par starts with 1 so add 1
  arglist[0] = PPM + 1; // par
  arglist[1] = 0.0;     // low
  arglist[2] = 100.;    // high
  // gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
  gMinuit->mnexcm("FIX", arglist, 1, ierflg);

  // printModel(2. * hnorm[9]->GetMaximumBin(), &lpar[0]);
  // return;

  /*
printf(" \n\n >>> tbFit modelFit start parameters fit ppm %f \n", vstart[PPM]);
for (int ii = 0; ii < NPARS; ++ii)
{
  gMinuit->GetParameter(ii, currentValue, currentError);
  printf("\t  param %i %s %.4E  \n", ii, lparNames[ii].Data(), currentValue);
}
*/
  double amin;
  printf("call mnprin starting values \n");
  gMinuit->mnprin(1, amin);

  /* look at sarting model */

  // void printModel(int ibin, Double_t *par, double *fsChan, double *ftChan)
  // printModel(hnorm[ichan]->GetMaximumBin(), &vstart[0]);

  /******************
   *  now fit
   *****************/

  double fval = 0;
  double gin[NPARS];
  int npar = NPARS;
  // Call the function once and get the return value.
  int llist = NPARS; // Number of parameters
  fcn(llist, gin, fval, lpar, ierflg);
  printf(" starting value >>>>   fval %E \n", fval);
  double fvalStart = fval;
  if (isnan(fval))
  {
    printf("gMinuit returns NAN\n");
    return;
  }

  // minimize with MIGRADfill
  // Now ready for minimization step
  arglist[0] = 1000000; // maxcalls
  arglist[1] = 1.E-5;   // tolerance

  /* MIGrad[maxcalls][tolerance]*/
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  printf("\n...  after fit, call to  MIGRAD returns ierflg %i \n", ierflg);

  // Print results
  Double_t edm, errdef;
  Int_t nvpar, nparx, icstat;
  printf("...  call mnstat \n");
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
  printf("\n...  call mnprin \n");
  gMinuit->mnprin(1, amin);

  /* fit model waves */
  hfitModel.resize(NCHANPMT);

  for (unsigned ic = 0; ic < NCHANPMT; ++ic)
  {
    TH1D *hFit = (TH1D *)hnorm[ic]->Clone(Form("fitWaveFitChan%i", ic));
    hFit->Reset("ICES");
    hFit->SetTitle((Form("fitWaveFitChan%i", ic)));
    hFit->SetLineColor(colors[ic]);
    hfitModel[ic] = hFit;
    fillFitWave(ic, hFit);
  }

  TDirectory *compDir = fout->mkdir("components");
  compDir->cd();

  // plot by channel first
  if (theFitChannel == -1)
  {
    for (unsigned ic = 0; ic < NCHANPMT; ++ic)
    {
      for (int icomp = 0; icomp < NUMCOMP; ++icomp)
      {
        TH1D *hFit = (TH1D *)hnorm[ic]->Clone(Form("fit%sChan%i", compNames[icomp].Data(), ic));
        hFit->Reset("ICES");
        hFit->SetTitle((Form("fit%sChan%i", compNames[icomp].Data(), ic)));
        hFit->SetLineColor(colors[ic]);
        fillCompWave(ic, icomp, hFit); // only need one of these
      }
    }
  }
  else
  {
    for (int icomp = 0; icomp < NUMCOMP; ++icomp)
    {
      TH1D *hFit = (TH1D *)hnorm[theFitChannel]->Clone(Form("fit%sChan%i", compNames[icomp].Data(), theFitChannel));
      hFit->Reset("ICES");
      hFit->SetTitle((Form("fit%sChan%i", compNames[icomp].Data(), theFitChannel)));
      hFit->SetLineColor(colors[theFitChannel]);
      fillCompWave(theFitChannel, icomp, hFit);
    }
  }

  ntParScan = new TNtuple("ntParScan", "parameter scan", "nll:fitVal1:fitVal2:fitVal3");

  int thePar = TAU3;
  printf("scan parameter %i %s from %f to %f \n", thePar, lparNames[thePar].Data(), 0.001 * lpar[thePar], 2. * lpar[thePar]);

  TGraph *graph = myScan(thePar, 0.001 * lpar[thePar], 2. * lpar[thePar]);
  graph->SetName(Form("ScanPar%i", thePar));
  graph->SetTitle(Form("ScanPar%i", thePar));
  graph->GetYaxis()->SetTitle("FCN likelihood value");
  graph->GetXaxis()->SetTitle(Form("parameter %s", lparNames[thePar].Data()));
  fout->Add(graph);

  if (theFitChannel == -1)
  {
    makeCanFit(0, 2, TString("canFitLevel0"));
    makeCanFit(3, 5, TString("canFitLevel1"));
    makeCanFit(6, 8, TString("canFitLevel2"));
    makeCanFit(9, 11, TString("canFitTrig"));
  }

  printf("\n...  finished tbFit \n");
}
