// lar_analytic.C
//
// CERN ROOT/C++ port of model_with_absorption_and_convolution7.py: closed-form
// (analytic) solution of the linear LAr/Xe kinetics system, the Neumeier
// cross-section fit, the emission-spectrum/QE/filter pipeline, time-domain
// smearing, and the full distance x ppm sweep with Neumeier 2.3mm replication.
//
// The ODE system (S,T,MS,MT,X,R) and its closed-form chain-of-exponentials
// solution are identical to the python script's lar_system / rates_from_params /
// lar_analytic_solution; see the header of lar_analyticVersion1.C (the
// previous, ODE-solver-only version of this file, kept as a reference) for a
// full derivation of the chain formula.
//
// Run:  root -l -b -q lar_analytic.C+        (ACLiC-compiled; recommended --
//                                              the full sweep is ~56 distance x
//                                              ppm combinations, each walking a
//                                              ~17000-point time grid)

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMinuit.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TStyle.h"
#include "TSystem.h"

// =========================================================
// 0. Generic helpers (CSV I/O, interpolation, trapezoid, Gaussian filtering)
// =========================================================

// Reads a 2-column CSV (wavelength, value). Matches python's
// np.genfromtxt(filename, delimiter=',', skip_header=1) -- the data files in
// this analysis have no header row, so skip_header=1 discards their first
// real data point; replicated here unchanged for fidelity with the existing
// (already-fit) analysis rather than "fixed".
bool ReadCSV(const char *filename, std::vector<double> &colA, std::vector<double> &colB)
{
   std::ifstream in(filename);
   if (!in.is_open())
      return false;
   std::string line;
   bool first = true;
   while (std::getline(in, line))
   {
      if (first)
      {
         first = false;
         continue;
      }
      if (line.empty())
         continue;
      std::replace(line.begin(), line.end(), ',', ' ');
      std::istringstream iss(line);
      double a, b;
      if (iss >> a >> b)
      {
         colA.push_back(a);
         colB.push_back(b);
      }
   }
   return true;
}

// Linear interpolation; extrapolate=true replicates kind='linear',
// fill_value='extrapolate'; extrapolate=false replicates bounds_error=False,
// fill_value=outOfRange (used by the polycarbonate/quartz transmission curves).
double LinearInterp(const std::vector<double> &x, const std::vector<double> &y, double xv,
                     bool extrapolate = true, double outOfRange = 0.0)
{
   const int n = (int)x.size();
   if (xv <= x[0])
      return extrapolate ? y[0] + (y[1] - y[0]) * (xv - x[0]) / (x[1] - x[0]) : outOfRange;
   if (xv >= x[n - 1])
      return extrapolate ? y[n - 2] + (y[n - 1] - y[n - 2]) * (xv - x[n - 2]) / (x[n - 1] - x[n - 2]) : outOfRange;
   for (int i = 0; i < n - 1; ++i)
      if (xv >= x[i] && xv <= x[i + 1])
         return y[i] + (y[i + 1] - y[i]) * (xv - x[i]) / (x[i + 1] - x[i]);
   return y[n - 1];
}

double Trapezoid(const std::vector<double> &y, const std::vector<double> &x)
{
   double sum = 0;
   for (size_t i = 1; i < x.size(); ++i)
      sum += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
   return sum;
}

std::vector<double> Linspace(double lo, double hi, int n)
{
   std::vector<double> v(n);
   for (int i = 0; i < n; ++i)
      v[i] = lo + (hi - lo) * i / (n - 1);
   return v;
}

// Matches scipy.ndimage.gaussian_filter1d(signal, sigma, mode='nearest'),
// default truncate=4.0: kernel radius = int(truncate*sigma + 0.5), boundary
// samples clamped to the nearest edge value.
std::vector<double> GaussianFilter1D(const std::vector<double> &input, double sigma, double truncate = 4.0)
{
   const int radius = (int)(truncate * sigma + 0.5);
   std::vector<double> kernel(2 * radius + 1);
   double sum = 0;
   for (int i = -radius; i <= radius; ++i)
   {
      kernel[i + radius] = std::exp(-0.5 * (i / sigma) * (i / sigma));
      sum += kernel[i + radius];
   }
   for (auto &k : kernel)
      k /= sum;
   const int n = (int)input.size();
   std::vector<double> out(n, 0.0);
   for (int idx = 0; idx < n; ++idx)
   {
      double acc = 0;
      for (int k = -radius; k <= radius; ++k)
      {
         int j = idx + k;
         if (j < 0)
            j = 0;
         else if (j >= n)
            j = n - 1;
         acc += input[j] * kernel[k + radius];
      }
      out[idx] = acc;
   }
   return out;
}

// Matches scipy.signal.convolve(signal, kernel, mode='same') with kernel =
// norm.pdf(t_kernel, 0, sigma) normalized to sum 1 (zero-padded boundary, NOT
// clamped -- this is a genuinely different edge behavior than GaussianFilter1D
// above, matching the python script's use of two different smoothing calls).
std::vector<double> ConvoluteTimeSmear(const std::vector<double> &signal, double dt, double sigma = 10.0)
{
   const int halfPts = (int)std::ceil(5.0 * sigma / dt);
   std::vector<double> kernel(2 * halfPts + 1);
   double sum = 0;
   for (int i = -halfPts; i <= halfPts; ++i)
   {
      const double t = i * dt;
      const double w = std::exp(-0.5 * (t / sigma) * (t / sigma)) / (sigma * std::sqrt(2 * M_PI));
      kernel[i + halfPts] = w;
      sum += w;
   }
   for (auto &k : kernel)
      k /= sum;
   const int n = (int)signal.size();
   std::vector<double> out(n, 0.0);
   for (int idx = 0; idx < n; ++idx)
   {
      double acc = 0;
      for (int k = -halfPts; k <= halfPts; ++k)
      {
         const int j = idx + k;
         if (j >= 0 && j < n)
            acc += signal[j] * kernel[k + halfPts];
      }
      out[idx] = acc;
   }
   return out;
}

// =========================================================
// 1. Load Absorption Data & Perform Global Cross-Section Fit
// =========================================================
const double x_ref = 11.6;
const double C_DF = -0.025;

std::vector<double> gLamHr;                                   // shared 110-140nm, 3000-pt grid for the fit model
std::map<double, std::pair<std::vector<double>, std::vector<double>>> gAbsDataByC; // c -> (wl, trans_norm), masked [118,132]
std::vector<double> gFitLam, gFitC, gFitTrans;                // flattened fit dataset

double A1_cs, E1_cs, sigma1_cs, A2_cs, E2_cs, sigma2_cs, slit_fit; // fit results, set by FitAbsorptionCrossSection()

// Model: for each unique concentration in cData, builds the true cross-section
// transmission on gLamHr, slit-convolves it, then linearly interpolates back
// onto the (possibly scattered) lamData points at that concentration -- an
// exact port of python's convolved_gaussian_transmission().
std::vector<double> ConvolvedGaussianTransmission(const std::vector<double> &lamData, const std::vector<double> &cData,
                                                   const double *par)
{
   const double A1 = par[0], E1 = par[1], s1 = par[2], A2 = par[3], E2 = par[4], s2 = par[5], slitSigma = par[6];
   const double dlam = gLamHr[1] - gLamHr[0];
   const double slitPixels = slitSigma / dlam;

   std::vector<double> uniqueC;
   for (double c : cData)
      if (std::find(uniqueC.begin(), uniqueC.end(), c) == uniqueC.end())
         uniqueC.push_back(c);

   std::vector<double> result(lamData.size(), 0.0);
   for (double cval : uniqueC)
   {
      std::vector<double> tTrue(gLamHr.size());
      for (size_t i = 0; i < gLamHr.size(); ++i)
      {
         const double E = 1239.84 / gLamHr[i];
         const double cs1 = A1 * std::exp(-((E - E1) * (E - E1)) / (2 * s1 * s1));
         const double cs2 = A2 * std::exp(-((E - E2) * (E - E2)) / (2 * s2 * s2));
         tTrue[i] = std::exp(-(cs1 + cs2) * cval * x_ref);
      }
      std::vector<double> tConv = GaussianFilter1D(tTrue, slitPixels);
      for (size_t i = 0; i < lamData.size(); ++i)
         if (cData[i] == cval)
            result[i] = LinearInterp(gLamHr, tConv, lamData[i], true);
   }
   return result;
}

void AbsFitFCN(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
   std::vector<double> model = ConvolvedGaussianTransmission(gFitLam, gFitC, par);
   double chi2 = 0;
   for (size_t i = 0; i < model.size(); ++i)
   {
      const double d = model[i] - gFitTrans[i];
      chi2 += d * d;
   }
   f = chi2;
}

void FitAbsorptionCrossSection()
{
   printf("Loading data and fitting Gaussian cross-section to Neumeier data...\n");

   std::vector<double> wlRayleigh, tRayleigh;
   if (!ReadCSV("Neumeier_Rayleigh.csv", wlRayleigh, tRayleigh))
   {
      printf("Warning: Neumeier_Rayleigh.csv not found. Using dummy 0.90 flat baseline.\n");
      wlRayleigh = Linspace(110, 150, 100);
      tRayleigh.assign(100, 0.90);
   }

   std::vector<std::pair<double, std::string>> files = {
       {0.1, "Neumeier_absorption.csv"}, {1.0, "1ppmabsorption.csv"}, {10.0, "10ppmabsorption.csv"}, {1000.0, "1000ppmabsorption.csv"}};

   for (auto &pr : files)
   {
      const double cVal = pr.first;
      std::vector<double> wlRaw, transRaw;
      if (!ReadCSV(pr.second.c_str(), wlRaw, transRaw))
      {
         printf("Warning: %s not found, skipping.\n", pr.second.c_str());
         continue;
      }
      std::vector<double> wlMasked, transNormMasked;
      for (size_t i = 0; i < wlRaw.size(); ++i)
      {
         const double rayleighAtWl = LinearInterp(wlRayleigh, tRayleigh, wlRaw[i], true);
         const double transNorm = transRaw[i] / (rayleighAtWl + C_DF);
         if (wlRaw[i] >= 118 && wlRaw[i] <= 132)
         {
            wlMasked.push_back(wlRaw[i]);
            transNormMasked.push_back(transNorm);
            gFitLam.push_back(wlRaw[i]);
            gFitC.push_back(cVal);
            gFitTrans.push_back(transNorm);
         }
      }
      gAbsDataByC[cVal] = {wlMasked, transNormMasked};
   }

   gLamHr = Linspace(110, 140, 3000);

   const double initPar[7] = {50.0, 9.785, 0.05, 20.0, 10.16, 0.08, 0.6};
   const double lowPar[7] = {0.0, 9.60, 0.001, 0.0, 10.00, 0.001, 0.01};
   const double highPar[7] = {10000.0, 9.90, 0.500, 10000.0, 10.40, 0.500, 2.0};
   const char *parNames[7] = {"A1_cs", "E1_cs", "sigma1_cs", "A2_cs", "E2_cs", "sigma2_cs", "slit_sigma"};

   TMinuit minuit(7);
   minuit.SetFCN(AbsFitFCN);
   minuit.SetPrintLevel(-1);
   int ierflg = 0;
   for (int i = 0; i < 7; ++i)
      minuit.mnparm(i, parNames[i], initPar[i], 0.01 * (highPar[i] - lowPar[i]), lowPar[i], highPar[i], ierflg);
   double arglist[2] = {5000, 1.0};
   minuit.mnexcm("MIGRAD", arglist, 2, ierflg);

   double err;
   minuit.GetParameter(0, A1_cs, err);
   minuit.GetParameter(1, E1_cs, err);
   minuit.GetParameter(2, sigma1_cs, err);
   minuit.GetParameter(3, A2_cs, err);
   minuit.GetParameter(4, E2_cs, err);
   minuit.GetParameter(5, sigma2_cs, err);
   minuit.GetParameter(6, slit_fit, err);

   printf("Fit Complete. Slit Sigma: %.2f nm\n", slit_fit);
   printf("A1_cs = %.2f, E1_cs = %.2f, sigma1_cs = %.2f\n", A1_cs, E1_cs, sigma1_cs);
   printf("A2_cs = %.2f, E2_cs = %.2f, sigma2_cs = %.2f\n", A2_cs, E2_cs, sigma2_cs);
}

void PlotAbsorptionFit()
{
   TCanvas *can = new TCanvas("canAbsFit", "Slit-Convolved Gaussian Fit of Xenon Transmission", 900, 600);
   TMultiGraph *mg = new TMultiGraph();
   mg->SetTitle("Slit-Convolved Gaussian Fit of Xenon Transmission;Wavelength (nm);Normalized Transmission");
   TLegend *leg = new TLegend(0.70, 0.55, 0.98, 0.9);
   const std::map<double, int> colors = {{0.1, kBlack}, {1.0, kBlue}, {10.0, kGreen + 2}, {1000.0, kRed}};
   std::vector<double> lamSmooth = Linspace(118, 132, 800);
   const double popt[7] = {A1_cs, E1_cs, sigma1_cs, A2_cs, E2_cs, sigma2_cs, slit_fit};

   for (auto &kv : gAbsDataByC)
   {
      const double cVal = kv.first;
      const int color = colors.count(cVal) ? colors.at(cVal) : kBlack;
      TGraph *gData = new TGraph((int)kv.second.first.size(), kv.second.first.data(), kv.second.second.data());
      gData->SetMarkerStyle(20);
      gData->SetMarkerSize(0.6);
      gData->SetMarkerColor(color);
      mg->Add(gData, "P");
      leg->AddEntry(gData, Form("%.1f ppm Data", cVal), "p");

      std::vector<double> cSmooth(lamSmooth.size(), cVal);
      std::vector<double> transFit = ConvolvedGaussianTransmission(lamSmooth, cSmooth, popt);
      TGraph *gFit = new TGraph((int)lamSmooth.size(), lamSmooth.data(), transFit.data());
      gFit->SetLineColor(color);
      gFit->SetLineWidth(2);
      mg->Add(gFit, "L");
      leg->AddEntry(gFit, Form("%.1f ppm Convolved Fit", cVal), "l");
   }
   mg->Draw("A");
   mg->GetXaxis()->SetLimits(118, 132);
   mg->SetMinimum(-0.05);
   mg->SetMaximum(1.05);
   leg->Draw();
   can->Update();
   gSystem->mkdir("BTB", kTRUE);
   can->SaveAs("BTB/convolved_gaussian_transmission_fit.pdf");
}

// =========================================================
// 2. Constants (matches section 2 of the python script)
// =========================================================
const std::vector<double> ppmx = {0.0001, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0};
const std::vector<double> distanceArr = {1.0, 10.0, 20.0, 30.0};
const double ppmn = 0.0;

const double taus = 4, taut = 1184, taumix_T = 4700, taumix_S = 35, taumix_R = 250;
const double taux = 20, tauquench = 2250, tauquenchmix = 3200;

const double Dab = 2.2e-5, RArAr = 5.21e-8, RArArsinglet = 52.0e-8, RArXe = 6.18e-8, RArN = 4.41e-8;
const double Avogadro = 6.022e23, rhoAr = 1.3954, molmassAr = 39.948;
const double mrhoAr = rhoAr / molmassAr;
const double NTotal = 36000;

const double b_rad = 0.70;

const double N_nbrs_prompt = NTotal * 0.0005;
const double N_rec = NTotal * 0.0020;
const double f_vis_slow = 0.0025;

const double k_rec = 1.0 / 36.0;
const double tau_nbrs = 0.5;

const double R_singlet_triplet = 0.313;
const double Sfrac = 1.0 / (1.0 + (1.0 / R_singlet_triplet));

const double t_min = 1300.0, t_start = 1450.0, tmax = 10000.0;
const double window_start = 1400.0, window_end = 1510.0;
const double dtGrid = 0.5;

const double d_target = 10.0;
const double inf = std::numeric_limits<double>::infinity();

// =========================================================
// 3. Closed-form (analytic) solution of the lar_system ODE
// =========================================================
// S and T decay autonomously (no absorption factor on their own decay rate --
// Absorb only redirects already-decaying photons into MS/MT via the shared
// absorbed-photon rate P_abs = Absorb*((1/taus)*S + (1/taut)*T), split
// between MS/MT by the radiative branching fraction b_rad; it does not slow
// S/T's own decay). MS is driven by S and T, MT by S, T and MS, X by MS+MT,
// and R decays autonomously. Every state is therefore a sum of "chain" terms
//   chain(l1,...,ln; t) = sum_i exp(-li*t) / prod_{j!=i} (lj - li)
// (the n-fold time-convolution of pure exponentials).
struct LArRates
{
   double lS, lT, lMS, lMT, lX;
   double aSS, aTS, aST, aTT;
   double b, g, k;
};
struct LArState
{
   double S, T, MS, MT, X, R;
};

double LArInv(double x) { return std::isinf(x) ? 0.0 : 1.0 / x; }

LArRates RatesFromTaus(double taudx_S, double taudx_T, double taudx_X, double taudn,
                        double Absorb)
{
   LArRates r;
   r.lS = LArInv(taus) + LArInv(taudx_S) + LArInv(taudn) + LArInv(tauquench);
   r.lT = LArInv(taut) + LArInv(taudx_T) + LArInv(taudn) + LArInv(tauquench);
   r.lMS = LArInv(taumix_S) + LArInv(taumix_R) + LArInv(taudx_X) + LArInv(taudn) + LArInv(tauquenchmix);
   r.lMT = LArInv(taumix_T) + LArInv(taudx_X) + LArInv(taudn) + LArInv(tauquenchmix);
   r.lX = LArInv(taux) + LArInv(taudn) + LArInv(tauquench);
   r.aSS = b_rad * Absorb * LArInv(taus) + LArInv(taudx_S);
   r.aTS = b_rad * Absorb * LArInv(taut);
   r.aST = (1.0 - b_rad) * Absorb * LArInv(taus);
   r.aTT = (1.0 - b_rad) * Absorb * LArInv(taut) + LArInv(taudx_T);
   r.b = LArInv(taumix_R);
   r.g = LArInv(taudx_X);
   r.k = k_rec;
   return r;
}

double Chain(const std::vector<double> &rates, double t)
{
   double out = 0.0;
   for (size_t i = 0; i < rates.size(); ++i)
   {
      double denom = 1.0;
      for (size_t j = 0; j < rates.size(); ++j)
         if (j != i)
            denom *= (rates[j] - rates[i]);
      out += std::exp(-rates[i] * t) / denom;
   }
   return out;
}

LArState AnalyticSolution(double t, const LArState &y0, const LArRates &r)
{
   const double S0 = y0.S, T0 = y0.T, MS0 = y0.MS, MT0 = y0.MT, X0 = y0.X, R0 = y0.R;
   LArState y;
   y.S = S0 * Chain({r.lS}, t);
   y.T = T0 * Chain({r.lT}, t);
   y.R = R0 * Chain({r.k}, t);
   y.MS = MS0 * Chain({r.lMS}, t) + r.aSS * S0 * Chain({r.lS, r.lMS}, t) + r.aTS * T0 * Chain({r.lT, r.lMS}, t);
   y.MT = MT0 * Chain({r.lMT}, t) + r.aST * S0 * Chain({r.lS, r.lMT}, t) + r.aTT * T0 * Chain({r.lT, r.lMT}, t) +
          r.b * MS0 * Chain({r.lMS, r.lMT}, t) + r.b * r.aSS * S0 * Chain({r.lS, r.lMS, r.lMT}, t) +
          r.b * r.aTS * T0 * Chain({r.lT, r.lMS, r.lMT}, t);
   y.X = X0 * Chain({r.lX}, t) +
         r.g * (MS0 * Chain({r.lMS, r.lX}, t) + r.aSS * S0 * Chain({r.lS, r.lMS, r.lX}, t) +
                r.aTS * T0 * Chain({r.lT, r.lMS, r.lX}, t)) +
         r.g * (MT0 * Chain({r.lMT, r.lX}, t) + r.aST * S0 * Chain({r.lS, r.lMT, r.lX}, t) +
                r.aTT * T0 * Chain({r.lT, r.lMT, r.lX}, t) + r.b * MS0 * Chain({r.lMS, r.lMT, r.lX}, t) +
                r.b * r.aSS * S0 * Chain({r.lS, r.lMS, r.lMT, r.lX}, t) +
                r.b * r.aTS * T0 * Chain({r.lT, r.lMS, r.lMT, r.lX}, t));
   return y;
}

struct LArLight
{
   double S, T, M, X, R;
};

// Converts populations at a given time into detected light per unit time.
// transFactor is the fraction of S/T's own radiative decay that escapes
// detection (1-Absorb in the main sweep, Tr_detector in the Neumeier
// section); tSinceStart = t - t_start, used to gate the prompt NBrS term
// (zero before the excitation time).
LArLight LightComponents(const LArState &y, double transFactor, double tSinceStart)
{
   LArLight light;
   light.S = transFactor * (1 / taus) * y.S;
   light.T = transFactor * (1 / taut) * y.T;
   light.M = (1 / taumix_S) * y.MS + (1 / taumix_T) * y.MT;
   light.X = (1 / taux) * y.X;
   const double promptNbrs = (tSinceStart >= 0) ? (N_nbrs_prompt / tau_nbrs) * std::exp(-tSinceStart / tau_nbrs) : 0.0;
   light.R = promptNbrs + k_rec * y.R + f_vis_slow * (1 / taut) * y.T;
   return light;
}

// Lightweight RK4 cross-check of the closed form above (no plots), run once
// at a representative operating point as a build-time sanity check.
void CheckAnalyticSolutionAgainstRK4()
{
   const double taudx_S = 50.0, taudx_T = 3000.0, taudx_X = 800.0, taudn = inf;
   const LArRates r = RatesFromTaus(taudx_S, taudx_T, taudx_X, taudn, /*Absorb*/ 0.3);
   const LArState y0 = {Sfrac * NTotal, (1.0 - Sfrac) * NTotal, 0.0, 0.0, 0.0, N_rec};

   auto derivs = [&](const LArState &y) -> LArState {
      return {-r.lS * y.S, -r.lT * y.T,
              r.aSS * y.S + r.aTS * y.T - r.lMS * y.MS,
              r.aST * y.S + r.aTT * y.T + r.b * y.MS - r.lMT * y.MT,
              r.g * (y.MS + y.MT) - r.lX * y.X, -r.k * y.R};
   };
   auto addScaled = [](const LArState &a, const LArState &d, double h) -> LArState {
      return {a.S + h * d.S, a.T + h * d.T, a.MS + h * d.MS, a.MT + h * d.MT, a.X + h * d.X, a.R + h * d.R};
   };

   const double h = 0.05;
   const int nSteps = 50000; // 2500 ns
   LArState yNum = y0;
   double maxDiff[6] = {0}, maxVal[6] = {0};
   for (int n = 0; n <= nSteps; ++n)
   {
      const double t = n * h;
      const LArState yA = AnalyticSolution(t, y0, r);
      const double a[6] = {yA.S, yA.T, yA.MS, yA.MT, yA.X, yA.R};
      const double v[6] = {yNum.S, yNum.T, yNum.MS, yNum.MT, yNum.X, yNum.R};
      for (int i = 0; i < 6; ++i)
      {
         maxDiff[i] = std::fmax(maxDiff[i], std::fabs(a[i] - v[i]));
         maxVal[i] = std::fmax(maxVal[i], std::fabs(v[i]));
      }
      if (n < nSteps)
      {
         const LArState k1 = derivs(yNum);
         const LArState k2 = derivs(addScaled(yNum, k1, 0.5 * h));
         const LArState k3 = derivs(addScaled(yNum, k2, 0.5 * h));
         const LArState k4 = derivs(addScaled(yNum, k3, h));
         yNum = {yNum.S + h / 6.0 * (k1.S + 2 * k2.S + 2 * k3.S + k4.S),
                 yNum.T + h / 6.0 * (k1.T + 2 * k2.T + 2 * k3.T + k4.T),
                 yNum.MS + h / 6.0 * (k1.MS + 2 * k2.MS + 2 * k3.MS + k4.MS),
                 yNum.MT + h / 6.0 * (k1.MT + 2 * k2.MT + 2 * k3.MT + k4.MT),
                 yNum.X + h / 6.0 * (k1.X + 2 * k2.X + 2 * k3.X + k4.X),
                 yNum.R + h / 6.0 * (k1.R + 2 * k2.R + 2 * k3.R + k4.R)};
      }
   }
   const char *names[6] = {"S", "T", "MS", "MT", "X", "R"};
   printf("Analytic-vs-RK4 sanity check (h=%.2f ns, %d steps):\n", h, nSteps);
   for (int i = 0; i < 6; ++i)
      printf("  %-3s max |analytic-numeric| / max = %.1e\n", names[i], maxDiff[i] / maxVal[i]);
}

// =========================================================
// 4. Emission profile, QE/filter transmission, and smearing functions
// =========================================================
const double A_E = 8132.0, mu_E = 128.2, sigma_E = 3.006;
const double Rayleigh_baseline = 1.00;

inline double Gaussian(double xv, double amp, double mu, double sigma)
{
   return amp * std::exp(-((xv - mu) * (xv - mu)) / (2 * sigma * sigma));
}

std::vector<double> GenerateEmissionSpectrum(const std::vector<double> &wavelengths, double y128, double y150,
                                              double y175, double y500, double x, double c)
{
   const size_t n = wavelengths.size();
   std::vector<double> shape128(n, 0.0);
   for (size_t i = 0; i < n; ++i)
   {
      const double wl = wavelengths[i];
      if (wl >= 110 && wl <= 900)
      {
         const double em128 = Gaussian(wl, A_E, mu_E, sigma_E);
         const double E = 1239.84 / wl;
         const double cs1 = A1_cs * std::exp(-((E - E1_cs) * (E - E1_cs)) / (2 * sigma1_cs * sigma1_cs));
         const double cs2 = A2_cs * std::exp(-((E - E2_cs) * (E - E2_cs)) / (2 * sigma2_cs * sigma2_cs));
         const double tScaled = std::exp(-(cs1 + cs2) * c * x) * Rayleigh_baseline;
         shape128[i] = em128 * tScaled;
      }
   }
   const double area128 = Trapezoid(shape128, wavelengths);
   if (area128 > 0)
      for (auto &v : shape128)
         v = v / area128 * y128;

   const double sig150 = 2.8, sig175 = 4.5;
   const double amp150 = y150 / (sig150 * std::sqrt(2 * M_PI));
   const double amp175 = y175 / (sig175 * std::sqrt(2 * M_PI));

   std::vector<double> shapeVis(n);
   for (size_t i = 0; i < n; ++i)
   {
      const double wl = wavelengths[i];
      const double nbrsContinuum = (wl >= 250) ? std::exp(-std::pow((wl - 560.0) / 220.0, 4)) : 0.0;
      const double nirHump = 0.6 * Gaussian(wl, 1.0, 760.0, 55.0);
      shapeVis[i] = nbrsContinuum + nirHump;
   }
   const double areaVis = Trapezoid(shapeVis, wavelengths);
   if (areaVis > 0)
      for (auto &v : shapeVis)
         v = v / areaVis * y500;

   std::vector<double> result(n);
   for (size_t i = 0; i < n; ++i)
      result[i] = shape128[i] + Gaussian(wavelengths[i], amp150, 149, sig150) +
                  Gaussian(wavelengths[i], amp175, 175, sig175) + shapeVis[i];
   return result;
}

// Cubic interpolation via TSpline3 (natural-spline boundary) -- an
// approximation of scipy's "not-a-knot" cubic interp1d; the QE curve is smooth
// enough that the boundary-condition difference is negligible for this use.
std::vector<double> SipmDetectorQE(const std::vector<double> &wavelengths)
{
   static const double wl[] = {120, 128, 150, 175, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900};
   static const double pde[] = {0.12, 0.15, 0.20, 0.24, 0.28, 0.33, 0.38, 0.44, 0.48, 0.50, 0.46, 0.38, 0.28, 0.20, 0.14, 0.09, 0.05, 0.03, 0.01};
   const int n = 19;
   TSpline3 spline("sipmQE", const_cast<double *>(wl), const_cast<double *>(pde), n);
   std::vector<double> out(wavelengths.size());
   for (size_t i = 0; i < wavelengths.size(); ++i)
   {
      const double w = wavelengths[i];
      out[i] = (w < wl[0] || w > wl[n - 1]) ? 0.0 : spline.Eval(w);
   }
   return out;
}

std::vector<double> PolycarbonateTrans(const std::vector<double> &wavelengths)
{
   static const std::vector<double> wl = {120, 128, 150, 175, 200, 250, 300, 380, 410, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900};
   static const std::vector<double> trans = {0.00, 0.00, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85};
   std::vector<double> out(wavelengths.size());
   for (size_t i = 0; i < wavelengths.size(); ++i)
      out[i] = std::clamp(LinearInterp(wl, trans, wavelengths[i], false, 0.0), 0.0, 1.0);
   return out;
}

std::vector<double> QuartzTrans(const std::vector<double> &wavelengths)
{
   static const std::vector<double> wl = {120, 130, 140, 150, 155, 160, 165, 170, 175, 180, 200, 300, 400, 500, 600, 700, 800, 850, 900};
   static const std::vector<double> trans = {0.0036, 0.0022, 0.0014, 0.0035, .1344, .3840, .7008, .8605, .8862, .8946, .9084, .9247, .9289, .9307, .9318, .9324, .9329, .9331, .9331};
   std::vector<double> out(wavelengths.size());
   for (size_t i = 0; i < wavelengths.size(); ++i)
      out[i] = std::clamp(LinearInterp(wl, trans, wavelengths[i], false, 0.0), 0.0, 1.0);
   return out;
}

// =========================================================
// 5. Exact integral: total 128nm light transmitted through distance x at
//    concentration c, relative to I_0 (x=0, c=0) -- fine-grid trapezoid
//    substituting for scipy.integrate.quad (smooth Gaussian-shaped integrand).
// =========================================================
double Integrand128(double lam, double x, double c)
{
   const double em = Gaussian(lam, A_E, mu_E, sigma_E);
   const double E = 1239.84 / lam;
   const double cs1 = A1_cs * std::exp(-((E - E1_cs) * (E - E1_cs)) / (2 * sigma1_cs * sigma1_cs));
   const double cs2 = A2_cs * std::exp(-((E - E2_cs) * (E - E2_cs)) / (2 * sigma2_cs * sigma2_cs));
   return em * std::exp(-(cs1 + cs2) * c * x);
}

double CalculateTotalLight(double x, double c)
{
   const int N = 3000;
   const double lo = 110, hi = 140, h = (hi - lo) / N;
   double sum = 0.5 * (Integrand128(lo, x, c) + Integrand128(hi, x, c));
   for (int i = 1; i < N; ++i)
      sum += Integrand128(lo + i * h, x, c);
   return sum * h;
}

// =========================================================
// Plotting helpers
// =========================================================
void addLabel(double ppm, double davg)
{
   TLatex lbl;
   lbl.SetNDC();
   lbl.SetTextSize(0.035);
   lbl.DrawLatex(0.15, 0.85, Form("PPM = %g, d = %g cm", ppm, davg));
}

TGraph *MakeGraph(const std::vector<double> &x, const std::vector<double> &y, int color, int markerStyle = 0)
{
   TGraph *g = new TGraph((int)x.size(), x.data(), y.data());
   g->SetLineColor(color);
   g->SetLineWidth(2);
   if (markerStyle > 0)
   {
      g->SetMarkerStyle(markerStyle);
      g->SetMarkerColor(color);
      g->SetMarkerSize(1.0);
   }
   return g;
}

// =========================================================
// 6 & 7. Main distance x ppm sweep, and post-loop Neumeier replication
// =========================================================
void lar_analytic()
{
   gStyle->SetOptStat(0);
   gSystem->mkdir("BTB", kTRUE);

   const bool PLOT = true;
   const bool PLOT_Neu = true;
   const bool PLOT_ABS = true;

   FitAbsorptionCrossSection();
   if (PLOT_ABS)
      PlotAbsorptionFit();

   CheckAnalyticSolutionAgainstRK4();

   const double I_0 = CalculateTotalLight(0.0, 0.0);

   // t_eval: 1300 to 10000 ns in 0.5 ns steps; t_ode: the >= t_start=1450 subset
   std::vector<double> t_eval;
   for (double t = t_min; t <= tmax + 1e-9; t += dtGrid)
      t_eval.push_back(t);
   const int nEval = (int)t_eval.size();
   int postStart = 0;
   while (postStart < nEval && t_eval[postStart] < t_start)
      ++postStart;
   const int nOde = nEval - postStart;

   const size_t nD = distanceArr.size(), nP = ppmx.size();
   std::vector<std::vector<double>> window_integral_all(nD, std::vector<double>(nP, 0.0));
   std::vector<std::vector<double>> singlet_yield_all(nD, std::vector<double>(nP, 0.0));

   for (size_t dIdx = 0; dIdx < nD; ++dIdx)
   {
      const double davg = distanceArr[dIdx];
      std::vector<double> Total(nP), Singlet(nP), Triplet(nP), Mixed(nP), Xenon(nP), Recomb(nP), ConvolutedTotal(nP);

      for (size_t i = 0; i < nP; ++i)
      {
         const double p = ppmx[i];
         const double cX = 1e-6 * mrhoAr * Avogadro * p;
         const double cN = 1e-6 * mrhoAr * Avogadro * ppmn;

         const double taudx_T = 1.0 / (4 * M_PI * Dab * RArAr * cX * 1e-9);
         const double taudx_X = 1.0 / (4 * M_PI * Dab * RArXe * cX * 1e-9);
         const double taudx_S = 1.0 / (4 * M_PI * Dab * RArArsinglet * cX * 1e-9);
         const double taudn = (cN > 0) ? 1.0 / (4 * M_PI * Dab * RArN * cN * 1e-9) : inf;

         const double Tr128 = CalculateTotalLight(davg, p) / I_0;
         const double Absorb = 1.0 - Tr128;

         const LArState y0 = {Sfrac * NTotal, (1.0 - Sfrac) * NTotal, 0.0, 0.0, 0.0, N_rec};
         const LArRates r = RatesFromTaus(taudx_S, taudx_T, taudx_X, taudn, Absorb);

         std::vector<double> S_sol(nEval, 0.0), T_sol(nEval, 0.0), MS_sol(nEval, 0.0), MT_sol(nEval, 0.0),
             X_sol(nEval, 0.0), R_sol(nEval, 0.0);
         for (int j = 0; j < nOde; ++j)
         {
            const LArState y = AnalyticSolution(t_eval[postStart + j] - t_start, y0, r);
            S_sol[postStart + j] = y.S;
            T_sol[postStart + j] = y.T;
            MS_sol[postStart + j] = y.MS;
            MT_sol[postStart + j] = y.MT;
            X_sol[postStart + j] = y.X;
            R_sol[postStart + j] = y.R;
         }

         std::vector<double> SLight(nEval), TLight(nEval), MLight(nEval), XLight(nEval), RLight(nEval);
         for (int j = 0; j < nEval; ++j)
         {
            const LArState y = {S_sol[j], T_sol[j], MS_sol[j], MT_sol[j], X_sol[j], R_sol[j]};
            const LArLight light = LightComponents(y, 1 - Absorb, t_eval[j] - t_start);
            SLight[j] = light.S;
            TLight[j] = light.T;
            MLight[j] = light.M;
            XLight[j] = light.X;
            RLight[j] = light.R;
         }

         Singlet[i] = Trapezoid(SLight, t_eval);
         Triplet[i] = Trapezoid(TLight, t_eval);
         Mixed[i] = Trapezoid(MLight, t_eval);
         Xenon[i] = Trapezoid(XLight, t_eval);
         Recomb[i] = Trapezoid(RLight, t_eval);
         Total[i] = Singlet[i] + Triplet[i] + Mixed[i] + Xenon[i] + Recomb[i];

         std::vector<double> wavelengths = Linspace(100, 900, 1000);
         std::vector<double> emissionCurve =
             GenerateEmissionSpectrum(wavelengths, Singlet[i] + Triplet[i], Mixed[i], Xenon[i], Recomb[i], davg, p);

         if (PLOT)
         {
            TCanvas can(Form("canEmission_%zu_%zu", dIdx, i), "emission", 800, 500);
            can.SetLogy();
            TGraph g((int)wavelengths.size(), wavelengths.data(), emissionCurve.data());
            g.SetLineColor(kBlue);
            g.SetTitle(Form("emission vs. wavelength, d = %g cm, ppmx = %g;Wavelength (nm);Yield", davg, p));
            g.GetXaxis()->SetLimits(100, 600);
            g.SetMinimum(0.005);
            g.SetMaximum(30000);
            g.Draw("AL");
            can.SaveAs(Form("BTB/emission_ppmx_%g_davg_%gcm.pdf", p, davg));
         }

         std::vector<double> qeCurve = SipmDetectorQE(wavelengths);
         std::vector<double> detectedSpectrum(wavelengths.size());
         for (size_t k = 0; k < wavelengths.size(); ++k)
            detectedSpectrum[k] = emissionCurve[k] * qeCurve[k];
         ConvolutedTotal[i] = Trapezoid(detectedSpectrum, wavelengths);

         auto bandQE = [&](double lo, double hi) {
            std::vector<double> wlM, emM, detM;
            for (size_t k = 0; k < wavelengths.size(); ++k)
               if (wavelengths[k] > lo && wavelengths[k] < hi)
               {
                  wlM.push_back(wavelengths[k]);
                  emM.push_back(emissionCurve[k]);
                  detM.push_back(detectedSpectrum[k]);
               }
            return Trapezoid(detM, wlM) / Trapezoid(emM, wlM);
         };
         const double SiPMQE128 = bandQE(119, 140);
         const double SiPMQE150 = bandQE(140, 160);
         const double SiPMQE175 = bandQE(160, 180);
         const double SiPMQE500 = bandQE(280, 900);

         std::vector<double> totTime(nEval);
         for (int j = 0; j < nEval; ++j)
            totTime[j] = SiPMQE128 * SLight[j] + SiPMQE128 * TLight[j] + SiPMQE175 * XLight[j] +
                         SiPMQE150 * MLight[j] + SiPMQE500 * RLight[j];

         std::vector<double> totTimeConv = ConvoluteTimeSmear(totTime, dtGrid, 10.0);
         std::vector<double> sLightScaled(nEval);
         for (int j = 0; j < nEval; ++j)
            sLightScaled[j] = SiPMQE128 * SLight[j];
         std::vector<double> slightConv = ConvoluteTimeSmear(sLightScaled, dtGrid, 10.0);

         std::vector<double> tWin, totConvWin, slightConvWin;
         for (int j = 0; j < nEval; ++j)
            if (t_eval[j] >= window_start && t_eval[j] <= window_end)
            {
               tWin.push_back(t_eval[j]);
               totConvWin.push_back(totTimeConv[j]);
               slightConvWin.push_back(slightConv[j]);
            }
         window_integral_all[dIdx][i] = Trapezoid(totConvWin, tWin);
         singlet_yield_all[dIdx][i] = Trapezoid(slightConvWin, tWin);

         if (PLOT)
         {
            const double duration = tmax - t_start;
            const double xlo = t_start - 0.1 * duration, xhi = t_start + 0.9 * duration;

            TCanvas can1(Form("canYieldTime_%zu_%zu", dIdx, i), "SiPM yield vs time", 800, 500);
            can1.SetLogy();
            TMultiGraph mg;
            std::vector<double> sS(nEval), sT(nEval), sX(nEval), sM(nEval), sR(nEval);
            for (int j = 0; j < nEval; ++j)
            {
               sS[j] = SiPMQE128 * SLight[j];
               sT[j] = SiPMQE128 * TLight[j];
               sX[j] = SiPMQE175 * XLight[j];
               sM[j] = SiPMQE150 * MLight[j];
               sR[j] = SiPMQE500 * RLight[j];
            }
            mg.Add(MakeGraph(t_eval, sS, kBlue));
            mg.Add(MakeGraph(t_eval, sT, kRed));
            mg.Add(MakeGraph(t_eval, sX, kGreen + 2));
            mg.Add(MakeGraph(t_eval, sM, kOrange + 1));
            mg.Add(MakeGraph(t_eval, sR, kMagenta));
            mg.Add(MakeGraph(t_eval, totTime, kBlack));
            mg.SetTitle(Form("SiPM yield vs. time, d = %g cm, ppmx = %g;Time (ns);Yield", davg, p));
            mg.Draw("AL");
            mg.GetXaxis()->SetLimits(xlo, xhi);
            mg.SetMinimum(0.005);
            mg.SetMaximum(1000.0);
            can1.SaveAs(Form("BTB/yield_time_SiPM_ppmx_%g_davg_%gcm.pdf", p, davg));

            TCanvas can2(Form("canConvYieldTime_%zu_%zu", dIdx, i), "Convoluted SiPM yield vs time", 800, 500);
            can2.SetLogy();
            TGraph gConv = *MakeGraph(t_eval, totTimeConv, kBlack);
            gConv.SetTitle(Form("Convoluted SiPM yield vs. time, d = %g cm, ppmx = %g;Time (ns);Yield", davg, p));
            gConv.Draw("AL");
            gConv.GetXaxis()->SetLimits(xlo, xhi);
            gConv.SetMinimum(0.005);
            gConv.SetMaximum(1000.0);
            can2.SaveAs(Form("BTB/conv_yield_time_SiPM_ppmx_%g_davg_%gcm.pdf", p, davg));

            if (std::fabs(davg - 10.0) < 1e-9)
            {
               std::vector<double> qCurve = QuartzTrans(wavelengths);
               std::vector<double> detSpecQ(wavelengths.size());
               for (size_t k = 0; k < wavelengths.size(); ++k)
                  detSpecQ[k] = detectedSpectrum[k] * qCurve[k];
               auto bandQEq = [&](double lo, double hi) {
                  std::vector<double> wlM, emM, dqM;
                  for (size_t k = 0; k < wavelengths.size(); ++k)
                     if (wavelengths[k] > lo && wavelengths[k] < hi)
                     {
                        wlM.push_back(wavelengths[k]);
                        emM.push_back(emissionCurve[k]);
                        dqM.push_back(detSpecQ[k]);
                     }
                  return Trapezoid(dqM, wlM) / Trapezoid(emM, wlM);
               };
               const double QE128q = bandQEq(119, 140), QE150q = bandQEq(140, 160), QE175q = bandQEq(160, 180), QE500q = bandQEq(280, 900);
               std::vector<double> totTimeQ(nEval);
               for (int j = 0; j < nEval; ++j)
                  totTimeQ[j] = QE128q * SLight[j] + QE128q * TLight[j] + QE150q * MLight[j] + QE175q * XLight[j] + QE500q * RLight[j];
               std::vector<double> convQ = ConvoluteTimeSmear(totTimeQ, dtGrid, 10.0);
               TCanvas can3(Form("canQuartz_%zu", i), "Quartz-filtered SiPM yield", 800, 500);
               can3.SetLogy();
               TGraph gq = *MakeGraph(t_eval, convQ, kGreen + 2);
               gq.SetTitle(Form("Convoluted Quartz-Filtered SiPM yield vs. time, d = %g cm, ppmx = %g;Time (ns);Yield", davg, p));
               gq.Draw("AL");
               gq.GetXaxis()->SetLimits(xlo, xhi);
               gq.SetMinimum(1e-4);
               gq.SetMaximum(1000.0);
               can3.SaveAs(Form("BTB/conv_yield_time_Quartz_SiPM_ppmx_%g_davg_%gcm.pdf", p, davg));
            }

            if (std::fabs(davg - 30.0) < 1e-9)
            {
               std::vector<double> pCurve = PolycarbonateTrans(wavelengths);
               std::vector<double> detSpecP(wavelengths.size());
               for (size_t k = 0; k < wavelengths.size(); ++k)
                  detSpecP[k] = detectedSpectrum[k] * pCurve[k];
               auto bandQEp = [&](double lo, double hi) {
                  std::vector<double> wlM, emM, dpM;
                  for (size_t k = 0; k < wavelengths.size(); ++k)
                     if (wavelengths[k] > lo && wavelengths[k] < hi)
                     {
                        wlM.push_back(wavelengths[k]);
                        emM.push_back(emissionCurve[k]);
                        dpM.push_back(detSpecP[k]);
                     }
                  return Trapezoid(dpM, wlM) / Trapezoid(emM, wlM);
               };
               const double QE128p = bandQEp(119, 140), QE150p = bandQEp(140, 160), QE175p = bandQEp(160, 180), QE500p = bandQEp(280, 900);
               std::vector<double> totTimeP(nEval);
               for (int j = 0; j < nEval; ++j)
                  totTimeP[j] = QE128p * SLight[j] + QE128p * TLight[j] + QE150p * MLight[j] + QE175p * XLight[j] + QE500p * RLight[j];
               std::vector<double> convP = ConvoluteTimeSmear(totTimeP, dtGrid, 10.0);
               TCanvas can4(Form("canPolycarb_%zu", i), "Polycarbonate-filtered SiPM yield", 800, 500);
               can4.SetLogy();
               TGraph gp = *MakeGraph(t_eval, convP, kMagenta);
               gp.SetTitle(Form("Convoluted Polycarbonate-Filtered SiPM yield vs. time, d = %g cm, ppmx = %g;Time (ns);Yield", davg, p));
               gp.Draw("AL");
               gp.GetXaxis()->SetLimits(xlo, xhi);
               gp.SetMinimum(1e-5);
               gp.SetMaximum(10.0);
               can4.SaveAs(Form("BTB/conv_yield_time_Polycarb_SiPM_ppmx_%g_davg_%gcm.pdf", p, davg));
            }
         }
      }

      if (PLOT)
      {
         TCanvas can(Form("canYieldPpmx_%zu", dIdx), "Light yield vs Xenon conc.", 800, 500);
         can.SetLogx();
         TMultiGraph mg;
         mg.Add(MakeGraph(ppmx, Total, kBlack, 29));
         mg.Add(MakeGraph(ppmx, ConvolutedTotal, kMagenta + 2, 3));
         mg.Add(MakeGraph(ppmx, Singlet, kBlue, 20));
         mg.Add(MakeGraph(ppmx, Triplet, kRed, 21));
         mg.Add(MakeGraph(ppmx, Mixed, kOrange + 1, 33));
         mg.Add(MakeGraph(ppmx, Xenon, kGreen + 2, 22));
         for (int k = 0; k < 6; ++k)
            ((TGraph *)mg.GetListOfGraphs()->At(k))->SetLineStyle(0);
         mg.SetTitle(Form("Light yield vs. Xenon conc., d = %g cm;Xenon Concentration (ppm);Light Yield", davg));
         mg.Draw("AP");
         mg.GetXaxis()->SetLimits(0.00005, 50.0);
         TLegend leg(0.7, 0.55, 0.98, 0.9);
         leg.AddEntry(mg.GetListOfGraphs()->At(0), "Total (Static QE)", "p");
         leg.AddEntry(mg.GetListOfGraphs()->At(1), "Total (Convoluted QE)", "p");
         leg.AddEntry(mg.GetListOfGraphs()->At(2), "Singlet", "p");
         leg.AddEntry(mg.GetListOfGraphs()->At(3), "Triplet", "p");
         leg.AddEntry(mg.GetListOfGraphs()->At(4), "Mixed", "p");
         leg.AddEntry(mg.GetListOfGraphs()->At(5), "Xenon", "p");
         leg.Draw();
         can.SaveAs(Form("BTB/yield_ppmx_davg_%gcm.pdf", davg));
      }
   }

   // =========================================================
   // 7. Post-Loop Plotting & Neumeier Replication
   // =========================================================
   if (PLOT_Neu)
   {
      for (size_t dIdx = 0; dIdx < nD; ++dIdx)
      {
         const double davg = distanceArr[dIdx];
         const double normSinglet0 = (singlet_yield_all[dIdx][0] != 0) ? singlet_yield_all[dIdx][0] : 1.0;
         const double normWindow0 = (window_integral_all[dIdx][0] != 0) ? window_integral_all[dIdx][0] : 1.0;
         std::vector<double> normSinglet(nP), normWindow(nP);
         for (size_t i = 0; i < nP; ++i)
         {
            normSinglet[i] = singlet_yield_all[dIdx][i] / normSinglet0;
            normWindow[i] = window_integral_all[dIdx][i] / normWindow0;
         }
         TCanvas can(Form("canNorm_%zu", dIdx), "Normalized yield comparison", 800, 500);
         can.SetLogx();
         TMultiGraph mg;
         mg.Add(MakeGraph(ppmx, normSinglet, kBlue, 20));
         mg.Add(MakeGraph(ppmx, normWindow, kRed, 29));
         mg.SetTitle(Form("Normalized Yield Comparison vs Xenon Concentration (d = %g cm);Xenon Concentration (ppm);Normalized Light Yield", davg));
         mg.Draw("APL");
         TLegend leg(0.55, 0.7, 0.9, 0.9);
         leg.AddEntry(mg.GetListOfGraphs()->At(0), "Normalized Total Singlet Yield (Model)", "lp");
         leg.AddEntry(mg.GetListOfGraphs()->At(1), Form("Normalized Convoluted Window (%g-%g ns)", window_start, window_end), "lp");
         leg.Draw();
         can.SaveAs(Form("BTB/normalized_windowed_yield_%gcm.pdf", davg));
      }

      // Surviving 128nm light profile vs absorption
      {
         std::vector<double> lamWindow = Linspace(110, 145, 500);
         std::vector<double> em128(lamWindow.size());
         double maxEm = 0;
         for (size_t i = 0; i < lamWindow.size(); ++i)
         {
            em128[i] = A_E * std::exp(-((lamWindow[i] - mu_E) * (lamWindow[i] - mu_E)) / (2 * sigma_E * sigma_E));
            maxEm = std::fmax(maxEm, em128[i]);
         }
         for (auto &v : em128)
            v /= maxEm;

         TCanvas can("canSurviving128", "Surviving 128nm profile", 800, 500);
         TMultiGraph mg;
         TGraph *gOrig = MakeGraph(lamWindow, em128, kBlack);
         gOrig->SetLineStyle(2);
         mg.Add(gOrig);
         TLegend leg(0.68, 0.55, 0.98, 0.9);
         leg.AddEntry(gOrig, "Original 128 nm Emission", "l");

         const std::vector<double> plotConcentrations = {0.0000001, 0.0001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0};
         const int plotColors[8] = {kBlack, kGray + 2, kBlue, kCyan + 1, kGreen + 2, kOrange + 1, kRed, kMagenta + 2};
         for (size_t ci = 0; ci < plotConcentrations.size(); ++ci)
         {
            const double cPlot = plotConcentrations[ci];
            std::vector<double> surviving(lamWindow.size());
            for (size_t i = 0; i < lamWindow.size(); ++i)
            {
               const double E = 1239.84 / lamWindow[i];
               const double cs1 = A1_cs * std::exp(-((E - E1_cs) * (E - E1_cs)) / (2 * sigma1_cs * sigma1_cs));
               const double cs2 = A2_cs * std::exp(-((E - E2_cs) * (E - E2_cs)) / (2 * sigma2_cs * sigma2_cs));
               const double tScaled = std::exp(-(cs1 + cs2) * cPlot * d_target) * Rayleigh_baseline;
               surviving[i] = em128[i] * tScaled;
            }
            TGraph *g = MakeGraph(lamWindow, surviving, plotColors[ci]);
            mg.Add(g);
            leg.AddEntry(g, Form("%g ppm Xenon", cPlot), "l");
         }
         mg.SetTitle(Form("Surviving 128 nm Light Profile after Absorption (d = %g cm);Wavelength (nm);Relative Intensity", d_target));
         mg.Draw("AL");
         leg.Draw();
         can.SaveAs("BTB/surviving_128nm_shape.pdf");
      }

      // Neumeier 2.3mm full-spectrum replication
      printf("Generating Neumeier 2.3 mm Replication Plot...\n");
      {
         std::vector<double> neuLam = Linspace(115, 200, 2000);
         const double dlamNeu = neuLam[1] - neuLam[0];
         const double slitPixelsNeu = slit_fit / dlamNeu;
         const std::vector<double> neuConcs = {0.0001, 0.1, 1.0, 10.0, 100.0, 1000.0};
         const std::map<double, int> neuColors = {{0.0001, kBlack}, {0.1, kBlack}, {1.0, kRed}, {10.0, kGreen + 2}, {100.0, kBlue}, {1000.0, kMagenta}};
         const double d_neu = 0.23;
         const double L_bulk = 4.0;

         // t_eval becomes the post-trigger-only grid for this section, matching python's `t_eval = t_ode`
         std::vector<double> tEvalNeu(t_eval.begin() + postStart, t_eval.end());
         const int nEvalNeu = (int)tEvalNeu.size();

         TCanvas can("canNeumeier", "Neumeier 2.3mm Replication", 1000, 600);
         TMultiGraph mg;
         TLegend leg(0.75, 0.55, 0.98, 0.9);

         for (double p : neuConcs)
         {
            const double cX = 1e-6 * mrhoAr * Avogadro * p;
            const double cN = 1e-6 * mrhoAr * Avogadro * ppmn;
            const double taudx_T = 1.0 / (4 * M_PI * Dab * RArAr * cX * 1e-9);
            const double taudx_X = 1.0 / (4 * M_PI * Dab * RArXe * cX * 1e-9);
            const double taudx_S = 1.0 / (4 * M_PI * Dab * RArArsinglet * cX * 1e-9);
            const double taudn = (cN > 0) ? 1.0 / (4 * M_PI * Dab * RArN * cN * 1e-9) : inf;

            const double Tr128Global = CalculateTotalLight(L_bulk, p) / I_0;
            const double AbsorbGlobal = 1.0 - Tr128Global;
            const double TrDetector = CalculateTotalLight(d_neu, p) / I_0;

            const LArState y0 = {Sfrac * NTotal, (1.0 - Sfrac) * NTotal, 0.0, 0.0, 0.0, N_rec};
            const LArRates r = RatesFromTaus(taudx_S, taudx_T, taudx_X, taudn, AbsorbGlobal);

            std::vector<double> SLight(nEvalNeu), TLight(nEvalNeu), MLight(nEvalNeu), XLight(nEvalNeu), RLight(nEvalNeu);
            for (int j = 0; j < nEvalNeu; ++j)
            {
               const LArState y = AnalyticSolution(tEvalNeu[j] - t_start, y0, r);
               const LArLight light = LightComponents(y, TrDetector, tEvalNeu[j] - t_start);
               SLight[j] = light.S;
               TLight[j] = light.T;
               MLight[j] = light.M;
               XLight[j] = light.X;
               RLight[j] = light.R;
            }

            const double sYield = Trapezoid(SLight, tEvalNeu);
            const double tYield = Trapezoid(TLight, tEvalNeu);
            const double mYield = Trapezoid(MLight, tEvalNeu) * 0.4;
            const double xYield = Trapezoid(XLight, tEvalNeu);
            const double rYield = Trapezoid(RLight, tEvalNeu);

            std::vector<double> rawSpectrum = GenerateEmissionSpectrum(neuLam, sYield + tYield, mYield, xYield, rYield, d_neu, p);
            std::vector<double> convolvedSpectrum = GaussianFilter1D(rawSpectrum, slitPixelsNeu);

            const int color = neuColors.count(p) ? neuColors.at(p) : kBlack;
            TGraph *g = MakeGraph(neuLam, convolvedSpectrum, color);
            mg.Add(g);
            leg.AddEntry(g, Form("LAr %g ppm Xe", p), "l");
         }
         mg.SetTitle("Simulated Raw Emission Spectrum at 2.3 mm (Neumeier Replication);Wavelength (nm);Signal (Raw Simulated Photons)");
         mg.Draw("AL");
         mg.GetXaxis()->SetLimits(115, 200);
         leg.Draw();
         can.SaveAs("BTB/Neumeier_2.3mm_Replication.pdf");
      }
      printf("Done!\n");
   }
}
