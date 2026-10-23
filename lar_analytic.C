// lar_analytic.C
//
// CERN ROOT macro: closed-form solution of the linear LAr/Xe kinetics system
//
//   dS  = -lS*S
//   dT  = -lT*T
//   dMS = aS*S - lMS*MS
//   dMT = aT*T + b*MS - lMT*MT
//   dX  = g*(MS + MT) - lX*X
//   dR  = -k_rec*R
//
// with
//   lS  = 1/taus + 1/taudx_S + 1/taudn + 1/tauquench
//   lT  = 1/taut + 1/taudx_T + 1/taudn + 1/tauquench
//   lMS = 1/taumix_S + 1/taumix_R + 1/taudx_X + 1/taudn + 1/tauquenchmix
//   lMT = 1/taumix_T + 1/taudx_X + 1/taudn + 1/tauquenchmix
//   lX  = 1/taux + 1/taudn + 1/tauquench
//   aS  = Absorb/taus + 1/taudx_S
//   aT  = Absorb/taut + 1/taudx_T
//   b   = 1/taumix_R
//   g   = 1/taudx_X
//
// Every component is a sum of chain terms
//   C(l1,...,ln; t) = sum_i exp(-li t) / prod_{j!=i} (lj - li)
// (the time convolution of exp(-l1 t), ..., exp(-ln t)), with
//   integral_0^inf C dt = 1 / (l1*...*ln).
//
// Run:  root -l lar_analytic.C        (interpreted)
//       root -l lar_analytic.C+       (compiled with ACLiC)

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#ifndef LAR_NO_ROOT
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMultiGraph.h"
#endif
TFile *fout;
std::vector<double> startIdx;
std::vector<double> startVal;
double RsingeltTriplet = 0.25;                           // fraction of singlet/triplet light that is prompt (1/4 for LAr, 1/3 for LXe)
double Sfrac = RsingeltTriplet / (1. + RsingeltTriplet); // fraction of prompt light that is singlet (1/5 for LAr, 1/4 for LXe)
static double LY = 41.;                                  // Doke, April 2009 https://arxiv.org/abs/0910.4956v1
static double nPhotons = 60. * LY;
std::vector<TH1D *> hAnalytic; // analytic solution histograms for each component (S, T, MS, MT, X, R)

struct LArRates
{
   double lS, lT, lMS, lMT, lX; // total loss rates (1/ns)
   double aS, aT, b, g;         // feed rates (1/ns)
   double k;                    // recombination rate (1/ns)
};

struct LArState
{
   double S, T, MS, MT, X, R;
};

// 1/x, treating an infinite lifetime as a zero rate
double LArInv(double x) { return std::isinf(x) ? 0.0 : 1.0 / x; }

LArRates RatesFromTaus(double taus, double taut, double taudx_S, double taudx_T,
                       double taudx_X, double taudn, double tauquench, double tauquenchmix,
                       double taumix_S, double taumix_T, double taumix_R, double taux,
                       double k_rec, double Absorb)
{
   printf("RatesFromTaus inputs:\n");
   printf("  taus %.3E  taut %.3E taudx_S %.3E  taudx_T %.3E\n", taus, taut, taudx_S, taudx_T);
   printf("  taudx_X %.3E  taudn %.3E tauquench %.3E  tauquenchmix %.3E\n", taudx_X, taudn, tauquench, tauquenchmix);
   printf("  taumix_S %.3E  taumix_T %.3E  taumix_R %.3E  taux %.3E\n", taumix_S, taumix_T, taumix_R, taux);
   printf("  k_rec %.3E  Absorb %.3f \n", k_rec, Absorb);
   printf("  singlet to triplet %.3f singlet fraction %.3f\n", RsingeltTriplet, Sfrac);
   LArRates r;
   r.lS = LArInv(taus) + LArInv(taudx_S) + LArInv(taudn) + LArInv(tauquench);
   r.lS *= (1.0 - Absorb); // apply absorption factor to singlet loss rate
   r.lT = LArInv(taut) + LArInv(taudx_T) + LArInv(taudn) + LArInv(tauquench);
   r.lT *= (1.0 - Absorb); // apply absorption factor to triplet loss rate
   r.lMS = LArInv(taumix_S) + LArInv(taumix_R) + LArInv(taudx_X) + LArInv(taudn) + LArInv(tauquenchmix);
   r.lMT = LArInv(taumix_T) + LArInv(taudx_X) + LArInv(taudn) + LArInv(tauquenchmix);
   r.lX = LArInv(taux) + LArInv(taudn) + LArInv(tauquench);
   r.aS = Absorb * (LArInv(taus) + LArInv(taudx_S));
   r.aT = Absorb * (LArInv(taut) + LArInv(taudx_T));
   r.b = LArInv(taumix_R);
   r.g = LArInv(taudx_X);
   r.k = k_rec;
   startIdx.clear();
   startVal.clear();
   startIdx.push_back(0.0);
   startVal.push_back(r.lS);
   startIdx.push_back(1.);
   startVal.push_back(r.lT);
   startIdx.push_back(2.);
   startVal.push_back(r.lMS);
   startIdx.push_back(3.);
   startVal.push_back(r.lMT);
   startIdx.push_back(4.);
   startVal.push_back(r.lX);
   startIdx.push_back(5.);
   startVal.push_back(r.aS);
   startIdx.push_back(6.);
   startVal.push_back(r.aT);

   printf("Starting rates:\n");
   printf("  S %.3E  T %.3E  MS %.3E  MT %.3E  X %.3E  R %.3E\n",
          startVal[0], startVal[1], startVal[2], startVal[3], startVal[4], startVal[5]);

   return r;
}

// Convolution of exponentials exp(-l t) for l in rates (rates must be distinct)
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

// Populations at time t, measured from the start of integration
/* This is standard C++11 list-initialization: any place a std::vector<T> is expected, a brace-enclosed list of T-convertible values is automatically turned into one. */
LArState AnalyticSolution(double t, const LArState &y0, const LArRates &r)
{
   const double S0 = y0.S, T0 = y0.T, MS0 = y0.MS, MT0 = y0.MT, X0 = y0.X, R0 = y0.R;
   LArState y;
   y.S = S0 * Chain({r.lS}, t);
   y.T = T0 * Chain({r.lT}, t);
   y.R = R0 * Chain({r.k}, t);
   y.MS = MS0 * Chain({r.lMS}, t) + r.aS * S0 * Chain({r.lS, r.lMS}, t);
   y.MT = MT0 * Chain({r.lMT}, t) + r.aT * T0 * Chain({r.lT, r.lMT}, t) + r.b * MS0 * Chain({r.lMS, r.lMT}, t) + r.b * r.aS * S0 * Chain({r.lS, r.lMS, r.lMT}, t);
   y.X = X0 * Chain({r.lX}, t) + r.g * (MS0 * Chain({r.lMS, r.lX}, t) + r.aS * S0 * Chain({r.lS, r.lMS, r.lX}, t) + MT0 * Chain({r.lMT, r.lX}, t) + r.aT * T0 * Chain({r.lT, r.lMT, r.lX}, t) + r.b * MS0 * Chain({r.lMS, r.lMT, r.lX}, t) + r.b * r.aS * S0 * Chain({r.lS, r.lMS, r.lMT, r.lX}, t));
   // finish by mutiplying by mulitplying by rates double lS, lT, lMS, lMT, lX; // total loss rates (1/ns)
   y.S *= r.lS;
   y.T *= r.lT;
   y.MS *= r.lMS;
   y.MT *= r.lMT;
   y.X *= r.lX;
   y.R *= r.k;

   return y;
}

// Integrals from 0 to infinity of each population (multiply by 1/tau for light yield)
LArState TimeIntegrals(const LArState &y0, const LArRates &r)
{
   LArState I;
   I.S = y0.S / r.lS;
   I.T = y0.T / r.lT;
   I.MS = (y0.MS + r.aS * I.S) / r.lMS;
   I.MT = (y0.MT + r.aT * I.T + r.b * I.MS) / r.lMT;
   I.X = (y0.X + r.g * (I.MS + I.MT)) / r.lX;
   I.R = y0.R / r.k;
   return I;
}

// Right-hand side of the ODE system (used only for the numerical check)
LArState Derivs(const LArState &y, const LArRates &r)
{
   LArState d;
   d.S = -r.lS * y.S;
   d.T = -r.lT * y.T;
   d.MS = r.aS * y.S - r.lMS * y.MS;
   d.MT = r.aT * y.T + r.b * y.MS - r.lMT * y.MT;
   d.X = r.g * (y.MS + y.MT) - r.lX * y.X;
   d.R = -r.k * y.R;
   return d;
}

LArState AddScaled(const LArState &a, const LArState &d, double h)
{
   return {a.S + h * d.S, a.T + h * d.T, a.MS + h * d.MS,
           a.MT + h * d.MT, a.X + h * d.X, a.R + h * d.R};
}

// One classical Runge-Kutta 4 step
LArState RK4Step(const LArState &y, const LArRates &r, double h)
{
   LArState k1 = Derivs(y, r);
   LArState k2 = Derivs(AddScaled(y, k1, 0.5 * h), r);
   LArState k3 = Derivs(AddScaled(y, k2, 0.5 * h), r);
   LArState k4 = Derivs(AddScaled(y, k3, h), r);
   return {y.S + h / 6.0 * (k1.S + 2 * k2.S + 2 * k3.S + k4.S),
           y.T + h / 6.0 * (k1.T + 2 * k2.T + 2 * k3.T + k4.T),
           y.MS + h / 6.0 * (k1.MS + 2 * k2.MS + 2 * k3.MS + k4.MS),
           y.MT + h / 6.0 * (k1.MT + 2 * k2.MT + 2 * k3.MT + k4.MT),
           y.X + h / 6.0 * (k1.X + 2 * k2.X + 2 * k3.X + k4.X),
           y.R + h / 6.0 * (k1.R + 2 * k2.R + 2 * k3.R + k4.R)};
}

void lar_analytic()
{
   fout = new TFile("lar_analytic.root", "RECREATE");
   // Representative values (ns)
   const double inf = std::numeric_limits<double>::infinity();
   LArRates r = RatesFromTaus(/*taus*/ 4, /*taut*/ 900, /*taudx_S*/ 50.0, /*taudx_T*/ 3000.0,
                              /*taudx_X*/ 800.0, /*taudn*/ inf, /*tauquench*/ 7700,
                              /*tauquenchmix*/ 7700, /*taumix_S*/ 6.0, /*taumix_T*/ 4700.0,
                              /*taumix_R*/ 30.0, /*taux*/ 20, /*k_rec*/ 1.0 / 36.0,
                              /*Absorb*/ 0.3);

   // --- Check the closed form against RK4 integration ---
   const double tEnd = 8550.0;
   const double h = 0.05;
   const int nSteps = static_cast<int>(tEnd / h + 0.5);
   const int nPlot = 2000;
   const int plotEvery = nSteps / nPlot;

   const char *names7[7] = {"S", "T", "MS", "MT", "X", "RS", "RT"};
   const char *names[6] = {"S", "T", "MS", "MT", "X", "R"};
   double maxDiff[7] = {0}, maxVal[7] = {0};
   std::vector<double> tp;
   std::vector<double> yp[7];

   // const double S0 = y0.S, T0 = y0.T, MS0 = y0.MS, MT0 = y0.MT, X0 = y0.X, R0 = y0.R;
   const LArState y0 = {Sfrac * nPhotons, (1.0 - Sfrac) * nPhotons, 0.0, 0.0, 0.0, 36.0};
   printf("Initial population values in struct y0 \n");
   printf("  S %.3E  T %.3E  MS %.3E  MT %.3E  X %.3E  R %.3E\n", y0.S, y0.T, y0.MS, y0.MT, y0.X, y0.R);

   LArState yNum = y0;
   for (int n = 0; n <= nSteps; ++n)
   {
      const double t = n * h;
      const LArState yA = AnalyticSolution(t, y0, r);
      const double a[6] = {yA.S, yA.T, yA.MS, yA.MT, yA.X, yA.R};
      const double v[6] = {yNum.S, yNum.T, yNum.MS, yNum.MT, yNum.X, yNum.R};
      if (n < 6)
      {
         printf("t %.2f ns  S %.3E  T %.3E  MS %.3E  MT %.3E  X %.3E  R %.3E\n", t, a[0], a[1], a[2], a[3], a[4], a[5]);
      }
      for (int i = 0; i < 6; ++i)
      {
         maxDiff[i] = std::fmax(maxDiff[i], std::fabs(a[i] - v[i]));
         maxVal[i] = std::fmax(maxVal[i], std::fabs(v[i]));
      }
      if (n % plotEvery == 0)
      {
         tp.push_back(t);
         for (int i = 0; i < 6; ++i)
            yp[i].push_back(a[i]);
      }
      if (n < nSteps)
         yNum = RK4Step(yNum, r, h);
   }

   printf("Closed form vs RK4 (h = %.2f ns):\n", h);
   for (int i = 0; i < 6; ++i)
      printf("  %-3s max |analytic - numeric| / max = %.1e\n", names[i], maxDiff[i] / maxVal[i]);

   // --- Total time integrals ---
   const LArState I = TimeIntegrals(y0, r);
   printf("Integrals from 0 to infinity:\n");
   printf("  S %.6g  T %.6g  MS %.6g  MT %.6g  X %.6g  R %.6g\n",
          I.S, I.T, I.MS, I.MT, I.X, I.R);

#ifndef LAR_NO_ROOT

   // --- Graph of starting values, one point per component ---
   {
      TCanvas *cStart = new TCanvas("cLArStart", "LAr/Xe kinetics (start values)", 900, 600);
      cStart->SetLogy();
      cStart->SetGrid();
      TGraph *gInt = new TGraph(startVal.size(), &startIdx[0], &startVal[0]);
      gInt->SetName("gTimeIntegrals");
      gInt->SetTitle("Starting rates ;component;rate [1/ns]");
      gInt->SetMarkerStyle(21);
      gInt->SetMarkerSize(1.2);
      gInt->Draw("ap");
      TLatex lbl;
      lbl.SetTextAlign(21);
      lbl.SetTextSize(0.035);
      for (int i = 0; i < 7; ++i)
         lbl.DrawLatex(startIdx[i], startVal[i] * 1.3, names7[i]);
      cStart->Update();
      cStart->SaveAs("lar_analytic_start.pdf");
      fout->Append(gInt);
      fout->Append(cStart);
   }
   // --- Graph of the time integrals, one point per component ---
   {
      TCanvas *cInt = new TCanvas("cLArIntegrals", "LAr/Xe kinetics (time integrals)", 900, 600);
      cInt->SetLogy();
      cInt->SetGrid();
      double idx[6] = {0, 1, 2, 3, 4, 5};
      double ival[6] = {I.S, I.T, I.MS, I.MT, I.X, I.R};
      TGraph *gInt = new TGraph(6, idx, ival);
      gInt->SetName("gTimeIntegrals");
      gInt->SetTitle("Time integrals (0 to #infty);component;integral");
      gInt->SetMarkerStyle(21);
      gInt->SetMarkerSize(1.2);
      gInt->Draw("ap");
      TLatex lbl;
      lbl.SetTextAlign(21);
      lbl.SetTextSize(0.035);
      for (int i = 0; i < 6; ++i)
         lbl.DrawLatex(idx[i], ival[i] * 1.3, names[i]);
      cInt->Update();
      cInt->SaveAs("lar_analytic_integrals.pdf");
      fout->Append(gInt);
      fout->Append(cInt);
   }

   // --- Plot the populations ---
   TCanvas *c = new TCanvas("cLAr", "LAr/Xe kinetics (closed form)", 900, 600);
   c->SetLogy();
   TMultiGraph *mg = new TMultiGraph();
   mg->SetTitle("Populations vs. time (closed form);Time (ns);Population");
   TLegend *leg = new TLegend(0.75, 0.6, 0.9, 0.9);
   const int colors[6] = {kBlue, kRed, kOrange + 1, kMagenta, kGreen + 2, kBlack};
   // log scale: keep only values at or above this plotted floor (mg->SetMinimum below).
   // Values from 0 down through subnormal doubles are physically irrelevant tail decay,
   // but their enormous dynamic range (>100 orders of magnitude) overflows the PDF
   // renderer's log-axis transform if included, corrupting the whole canvas.
   const double plotFloor = 1e-3;
   for (int i = 0; i < 6; ++i)
   {
      std::vector<double> tPos, yPos;
      for (size_t j = 0; j < tp.size(); ++j)
         if (yp[i][j] > plotFloor && std::isnormal(yp[i][j]))
         {
            tPos.push_back(tp[j]);
            yPos.push_back(yp[i][j]);
         }
      if (tPos.empty())
         continue;
      TGraph *gr = new TGraph(tPos.size(), tPos.data(), yPos.data());
      gr->SetName(names[i]);
      fout->Append(gr);
      gr->SetLineColor(colors[i]);
      gr->SetLineWidth(2);
      mg->Add(gr, "L");
      leg->AddEntry(gr, names[i], "l");
   }
   mg->Draw("A");
   mg->SetMinimum(1e-3);
   leg->Draw();
   c->Update();
   c->SaveAs("lar_analytic.pdf");
   fout->Append(c);

   // --- Fill hAnalytic histograms (7500 bins, 2 ns/bin, 0-15000 ns) from the
   //     closed-form solution evaluated at each bin center ---
   const int nBins = 7500;
   const double binWidth = 2.0; // ns
   hAnalytic.resize(6);
   for (int i = 0; i < 6; ++i)
   {
      hAnalytic[i] = new TH1D(Form("hAnalytic%s", names[i]), Form("hAnalytic%s;time [ns];population", names[i]),
                              nBins, 0, nBins * binWidth);
      fout->Append(hAnalytic[i]);
   }
   for (int ibin = 1; ibin <= nBins; ++ibin)
   {
      const double t = hAnalytic[0]->GetBinCenter(ibin);
      const LArState yAhist = AnalyticSolution(t, y0, r);
      const double a[6] = {yAhist.S, yAhist.T, yAhist.MS, yAhist.MT, yAhist.X, yAhist.R};
      for (int i = 0; i < 6; ++i)
         hAnalytic[i]->SetBinContent(ibin, a[i]);
   }

   // --- Multiplot of the hAnalytic histograms themselves (overlaid, not converted to TGraph) ---
   gStyle->SetOptStat(0);
   TCanvas *cHist = new TCanvas("cLArHist", "LAr/Xe kinetics (from hAnalytic histograms)", 900, 600);
   cHist->SetLogy();
   TLegend *legHist = new TLegend(0.75, 0.6, 0.9, 0.9);
   double histMax = 0;
   for (int i = 0; i < 6; ++i)
      histMax = std::fmax(histMax, hAnalytic[i]->GetMaximum());
   for (int i = 0; i < 6; ++i)
   {
      hAnalytic[i]->SetLineColor(colors[i]);
      hAnalytic[i]->SetLineWidth(2);
      hAnalytic[i]->SetTitle("Populations vs. time (hAnalytic);Time (ns);Population");
      hAnalytic[i]->GetYaxis()->SetRangeUser(plotFloor, histMax * 1.2);
      fout->Append(hAnalytic[i]);
      hAnalytic[i]->Draw(i == 0 ? "hist" : "hist same");
      legHist->AddEntry(hAnalytic[i], names[i], "l");
   }
   legHist->Draw();
   gPad->RangeAxis(hAnalytic[0]->GetXaxis()->GetXmin(), plotFloor,
                   hAnalytic[0]->GetXaxis()->GetXmax(), hAnalytic[0]->GetMaximum() * 1.2);
   cHist->Update();
   cHist->SaveAs("lar_analytic_hist.pdf");
   fout->Append(cHist);
#endif

   fout->Write();
   // Disown all objects from fout before Close(): fout->Write() has already
   // persisted everything, and if objects are left in fout's list, Close()
   // will delete them a second time after ROOT's own canvas/object cleanup
   // already did so at macro exit, producing a segfault. Clear() removes the
   // list entries without deleting the objects themselves.
   fout->GetList()->Clear();
   fout->Close();
}
