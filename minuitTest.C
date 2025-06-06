#include "TMinuit.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TRandom3.h"

TRandom3 *ran;
TMinuit *gMinuit;

int ngen;
// global scope data arrays
static std::vector<double> xval;
static std::vector<double> yval;
static std::vector<double> yerr;

// function to minimize
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // sum negative log likelihood
  double val = 0;
  for (unsigned i = 0; i < ngen; i++)
  {
    double ymean = par[0] * xval[i] + par[1];
    double z = (yval[i] - ymean) / yerr[i];
    val += 0.5 * z * z;
    // printf("point %i x %f y %f par0  %f %f fcn val %f \n", i, xval[i], yval[i], par[0], par[1], val);
  }
  // printf("slope %f intercept %f fcn value %f \n", par[0], par[1], val);
  f = val;
}

void minuitTest(double slope = 1, double intercept = 0)
{
  TFile *fout = new TFile("minuitTest.root", "recreate");
  ngen = 100;
  double xlow = -1.;
  double xhigh = 1.;
  double sigma = 0.1;
  TNtuple *nt = new TNtuple("ntLine", "ntLine", "x:y:ymean");

  ran = new TRandom3();
  TF1 *f1 = new TF1("myLine", "[0]*x+[1]", xlow, xhigh);
  f1->SetParName(0, "slope");
  f1->SetParameter(0, slope);
  f1->SetParName(1, "intercept");
  f1->SetParameter(1, intercept);

  xval.clear();
  yval.clear();
  yerr.clear();

  // generate points about line
  for (int igen = 0; igen < ngen; ++igen)
  {
    double x = (xhigh - xlow) * ran->Rndm() + xlow;
    double ymean = f1->Eval(x);
    double y = ran->Gaus(ymean, sigma);
    xval.push_back(x);
    yval.push_back(y);
    yerr.push_back(sigma);
    nt->Fill(x, y, ymean);
  }
  TGraphErrors *gr = new TGraphErrors(ngen, &xval[0], &yval[0], nullptr, &yerr[0]);
  gr->SetName("ranLine");
  gr->SetTitle("random points about a Line");
  gr->Fit(f1, "myLine", "R+", xlow, xhigh);

  TCanvas *c = new TCanvas("minuitTest", "minuitTest");
  gr->Draw("ap");
  f1->Draw("same");
  // gStyle->SetOptStat(1001101);
  gStyle->SetOptFit(1011);
  c->Update();

  /****** now try with TMinuit ********/
  gMinuit = new TMinuit(2);
  printf("minimize with Minuit number of points %i \n", ngen);
  gMinuit->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 0.5; // UP for likelihood
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  // Set starting values and step sizes for parameters
  double step = 0.0001;
  /* the range must be specified unlike in bad example Ifit.C!! */
  gMinuit->mnparm(0, "slope", slope, step, -100. * slope, 100. * slope, ierflg);
  // gMinuit->mnparm(1, "intercept", intercept, step, -100. * intercept, 100. * intercept, ierflg); // because intercept value is zero cannot do this
  gMinuit->mnparm(1, "intercept", intercept, step, -100., 100., ierflg);
  // fix slope
  // gMinuit->FixParameter(1);

  double val, err;
  gMinuit->GetParameter(0, val, err);
  printf("starting par %i val %f +/- %f\n", 0, val, err);
  gMinuit->GetParameter(1, val, err);
  printf("starting par %i val %f +/- %f\n", 1, val, err);
  // Print results
  double amin, edm, errdef;
  int nvpar, nparx, icstat;
  gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  // Now ready for minimization step with MIGRAD
  gMinuit->mnexcm("MIGRAD", arglist, 0, ierflg);
  printf("\t\t ***** MIGRAD error code %i ******\n", ierflg);
  gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  double fval = 0;
  gMinuit->mnprin(0, fval);
  gMinuit->mnprin(1, fval);

  //  gMinuit->mnprin(3,amin);

  fout->Add(f1);
  fout->Add(gr);
  fout->Write();
}
