/*
  compare the full multi-exponential Tr128 (Tr128All) against the
  single-exponential approximation (Tr128Exp) used in tbDraw.cc
*/
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "distanceLevels.hh"
#include "modelAllFit.hh"

using namespace TMath;

double Tr128All(double ppm, double dist, double absorb1Const = absorb1ConstDefault)
{
    ppm = max(1.0E-9, ppm);
    double A = Aconstant;
    double C = Cconstant;
    double lambda1 = absorb1Const * 0.1 / ppm;
    double lambda2 = 3.38 * 0.1 / ppm;
    double lambda3 = 50.5 * 0.1 / ppm;
    return A * exp(-dist / lambda1) + C * exp(-dist / lambda2) + (1 - A - C) * exp(-dist / lambda3);
}

double Tr128Exp(double ppm, double dist, double absorb1Const = absorb1ConstDefault)
{
    ppm = max(1.0E-9, ppm);
    double A = Aconstant;
    double lambda1 = absorb1Const * 0.1 / ppm;
    // printf("Tr128Exp: ppm %.3f dist %.3f A %.3E lambda1 %.3E tr128 %.3f\n", ppm, dist, A, A * exp(-dist / lambda1));
    return A * exp(-dist / lambda1);
}

void drawAbsorption(double ppm = 1., double absorb1Const = absorb1ConstDefault)
{
    int npoints = 10000;
    std::vector<double> distance(npoints);
    std::vector<double> tr128All(npoints);
    std::vector<double> tr128Exp(npoints);
    for (int i = 0; i < npoints; ++i)
    {
        distance[i] = double(i) * .005;
        tr128All[i] = Tr128All(ppm, distance[i], absorb1ConstDefault);
        tr128Exp[i] = Tr128Exp(ppm, distance[i], absorb1ConstDefault);
        // printf("ppm %.3f dist %.3f Tr128All %.3E Tr128Exp %.3E\n", ppm, distance[i], tr128All[i], tr128Exp[i]);
    }

    TGraph *gAll = new TGraph(npoints, &distance[0], &tr128All[0]);
    gAll->SetName("Tr128All");
    gAll->SetLineColor(kBlue);
    gAll->SetLineWidth(2);

    TGraph *gExp = new TGraph(npoints, &distance[0], &tr128Exp[0]);
    gExp->SetName("Tr128Exp");
    gExp->SetLineColor(kRed);
    gExp->SetLineWidth(2);
    gExp->SetLineStyle(2);

    TMultiGraph *mg = new TMultiGraph("mgTr128Compare", Form("Tr128All vs Tr128Exp  PPM=%.3f absorb1Const=%.3E", ppm, absorb1Const));
    mg->Add(gAll, "l");
    mg->Add(gExp, "l");

    TCanvas *can = new TCanvas("canTr128Compare", "Tr128All vs Tr128Exp");
    can->SetGrid();
    can->SetLogx();
    mg->GetXaxis()->SetTitle("distance [cm]");
    mg->GetYaxis()->SetTitle("Tr128");
    mg->Draw("a");

    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(gAll, "Tr128All (3-exp)", "l");
    leg->AddEntry(gExp, "Tr128Exp (1-exp)", "l");
    leg->Draw();

    can->Print("Tr128Compare.pdf");
}
