
/*
look at  postAna root file
*/

// pass bit failures hex
#include "modelAllFit.hh"
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <numeric>
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TKey.h"
#include <TROOT.h>

using namespace std;
enum FAILURECODES
{
    PASS = 0,
    BASEFAIL = 0x1,
    EARLYCUT = 0x2,
    FIRSTTIME = 0x4,
    COSMIC = 0x8,
    GAMMA = 0x10,
    TRIGFAIL = 0x20,
    TRIANGLE = 0x40, // 2^6
    TOTALCODES = 2 * TRIANGLE
};

enum
{
    FAILBITS = 8
};

std::string sdate;
TFile *fin0;
TFile *fin1;
TFile *fout;
std::vector<TString> codeNames;
std::vector<int> failCodes;
std::vector<TString> bitNames;
std::vector<std::vector<double>> bitCount; // [bit][file]
std::vector<double> bitFile;               // file number
std::vector<std::vector<TH1D *>> hnorm;
std::vector<std::vector<TH1D *>> hcurve;
std::vector<std::vector<TH1D *>> hEffNorm;

/*
 * Fit histogram to exponential plus constant model
 * Model: f(x) = A * exp(-x/tau) + C
 *
 * Parameters:
 *   h - TH1* histogram to fit
 *   xmin - minimum x value for fit range
 *   xmax - maximum x value for fit range
 *   amp_init - initial amplitude guess (default 1.0)
 *   tau_init - initial decay constant guess (default 1.0)
 *   const_init - initial constant background guess (default 0.1)
 *
 * Returns:
 *   TF1* - fitted function (or nullptr if fit fails)
 */
static TF1 *fitExpPlusConst(TH1 *h, double xmin, double xmax,
                            double amp_init = 1.0, double tau_init = 1000.0,
                            double const_init = 0.0)
{
    if (!h)
    {
        printf("fitExpPlusConst: invalid histogram pointer\n");
        return nullptr;
    }

    // Create fit function: A*exp(-x/tau) + C
    TF1 *fitFunc = new TF1("fitExpPlusConst", "[0]*exp(-x/[1])+[2]", xmin, xmax);
    fitFunc->SetLineColor(kBlack);
    fitFunc->SetLineWidth(4);

    // Set parameter names for clarity
    fitFunc->SetParName(0, "Amplitude");
    fitFunc->SetParName(1, "Tau");
    fitFunc->SetParName(2, "Const");

    // Set initial parameter values
    fitFunc->SetParameter(0, amp_init);
    fitFunc->SetParameter(1, tau_init);
    fitFunc->SetParameter(2, const_init);

    // Set reasonable ranges to help convergence
    fitFunc->SetParLimits(0, 1.E-10, 1.0E4); // amplitude must be positive
    fitFunc->SetParLimits(1, 1.0E01, 1.0E6); // tau must be positive
    fitFunc->SetParLimits(2, 0, 1.0E-5);     // constant must be non-negative
    // fitFunc->FixParameter(2, 0);

    // Perform the fit // weighted likelihood
    // minos errors  IMPROVE algorithm cannor be computed in ML fit
    // M improved fit MIGRAD
    h->Fit(fitFunc, "SWLM+", "", xmin, xmax);

    // Print fit results
    printf("fitExpPlusConst: fit results for histogram %s\n", h->GetName());
    printf("  Amplitude = %.4E +/- %.4E\n", fitFunc->GetParameter(0), fitFunc->GetParError(0));
    printf("  Tau       = %.4E +/- %.4E\n", fitFunc->GetParameter(1), fitFunc->GetParError(1));
    printf("  Const     = %.4E +/- %.4E\n", fitFunc->GetParameter(2), fitFunc->GetParError(2));
    printf("  Chi-square = %.4f\n", fitFunc->GetChisquare());
    printf("  NDF = %d\n", fitFunc->GetNDF());

    return fitFunc;
}

/*
 * Fit histogram to double exponential plus constant model
 * Model: f(x) = A1 * exp(-x/tau1) + A2 * exp(-x/tau2) + C
 *
 * Parameters:
 *   h - TH1* histogram to fit
 *   xmin - minimum x value for fit range
 *   xmax - maximum x value for fit range
 *   amp1_init - initial amplitude 1 guess (default 1.0)
 *   tau1_init - initial decay constant 1 guess (default 100.0)
 *   amp2_init - initial amplitude 2 guess (default 1.0)
 *   tau2_init - initial decay constant 2 guess (default 1000.0)
 *   const_init - initial constant background guess (default 0.0)
 *
 * Returns:
 *   TF1* - fitted function (or nullptr if fit fails)
 */
static TF1 *fitDoubleExp(TH1 *h, double xmin, double xmax,
                         double amp1_init = 1.0, double tau1_init = 100.0,
                         double amp2_init = 1.0, double tau2_init = 1000.0)
{
    if (!h)
    {
        printf("fitDoubleExpPlusConst: invalid histogram pointer\n");
        return nullptr;
    }

    // Create fit function: A1*exp(-x/tau1) + A2*exp(-x/tau2) + C
    TF1 *fitFunc = new TF1("fitDoubleExpPlusConst", "[0]*exp(-x/[1])+[2]*exp(-x/[3])", xmin, xmax);
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(4);

    // Set parameter names for clarity
    fitFunc->SetParName(0, "Amplitude1");
    fitFunc->SetParName(1, "Tau1");
    fitFunc->SetParName(2, "Amplitude2");
    fitFunc->SetParName(3, "Tau2");

    // Set initial parameter values
    fitFunc->SetParameter(0, amp1_init);
    fitFunc->SetParameter(1, tau1_init);
    fitFunc->SetParameter(2, amp2_init);
    fitFunc->SetParameter(3, tau2_init);

    // Set reasonable ranges to help convergence
    fitFunc->SetParLimits(0, 1.E-10, 1.0);   // amplitude1 must be positive
    fitFunc->SetParLimits(1, 1.0E01, 1.0E6); // tau1 must be positive
    fitFunc->SetParLimits(2, 1.E-10, 1.0);   // amplitude2 must be positive
    fitFunc->SetParLimits(3, 1.0E01, 1.0E6); // tau2 must be positive

    // Perform the fit // weighted likelihood
    h->Fit(fitFunc, "SWLM+", "", xmin, xmax);

    // Print fit results
    printf("fitDoubleExpPlusConst: fit results for histogram %s\n", h->GetName());
    printf("  Amplitude1 = %.4E +/- %.4E\n", fitFunc->GetParameter(0), fitFunc->GetParError(0));
    printf("  Tau1       = %.4E +/- %.4E\n", fitFunc->GetParameter(1), fitFunc->GetParError(1));
    printf("  Amplitude2 = %.4E +/- %.4E\n", fitFunc->GetParameter(2), fitFunc->GetParError(2));
    printf("  Tau2       = %.4E +/- %.4E\n", fitFunc->GetParameter(3), fitFunc->GetParError(3));
    printf("  Chi-square = %.4f\n", fitFunc->GetChisquare());
    printf("  NDF = %d\n", fitFunc->GetNDF());

    return fitFunc;
}

void zeroErrors(TH1D *h)
{
    for (int ibin = 0; ibin < h->GetNbinsX(); ++ibin)
        h->SetBinError(ibin, 0);
}

string currentDate()
{
    time_t rawtime;
    struct tm *timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    char output[30];
    strftime(output, 30, "%Y-%m-%d-%H-%M", timeinfo);
    return string(output);
}

TCanvas *makeCanCompare(int ichan, TString canName)
{
    printf(" makeCanCompare chan %i name %s  \n", ichan, canName.Data());
    TCanvas *can = new TCanvas(canName, canName);
    // zeroErrors(hnorm[0][ichan]);
    // zeroErrors(hnorm[1][ichan]);
    hnorm[0][ichan]->GetYaxis()->SetRangeUser(1.E-9, 1.E-3);
    hnorm[1][ichan]->GetYaxis()->SetRangeUser(1.E-9, 1.E-3);
    hnorm[0][ichan]->Rebin(50);
    hnorm[1][ichan]->Rebin(50);
    hnorm[0][ichan]->SetLineColor(kRed);
    hnorm[1][ichan]->SetLineColor(kBlue);
    hnorm[0][ichan]->Draw("HIST");
    hnorm[1][ichan]->Draw("HISTSAME");
    can->BuildLegend();
    can->SetLogy();
    return can;
}

TCanvas *makeCanCutNorm(int i1, int i2, TString canName)
{
    printf(" makeCanCutfNorm %s from %i to %i size %lu \n", canName.Data(), i1, 12, hnorm[0].size());
    bool firstPlot = true;
    TCanvas *can = new TCanvas(canName, canName);
    for (int i = i1; i >= i2; --i)
    {
        // printf("%i %s \n", i, hnorm[i]->GetName());
        //  hnorm[i]->GetYaxis()->SetRangeUser(tiny, ymax);
        // hnorm[i]->GetYaxis()->SetRangeUser(0.1, 30);
        // hnorm[i]->GetXaxis()->SetRangeUser(1000, 4000);
        // printf("eff draw %i \n", i);
        if (firstPlot)
        {
            printf("%i %s \n", i, hnorm[0][i]->GetName());
            hnorm[0][i]->Draw("HIST");
            hnorm[1][i]->Draw("HISTSAME");
            firstPlot = false;
        }
        else if (!isBadChannel(i))
        {
            hnorm[0][i]->Draw("HISTSAME");
            hnorm[1][i]->Draw("HISTSAME");
        }
        else
            printf("skip bad channel %i %s \n", i, hnorm[0][i]->GetName());
    }
    can->BuildLegend();
    can->SetLogy();
    return can;
}

TCanvas *makeCanEffNorm(int i1, int i2, TString canName)
{
    printf(" makeCanEffNorm %s from %i to %i size %lu \n", canName.Data(), i1, 12, hEffNorm.size());
    bool firstPlot = true;
    TCanvas *can = new TCanvas(canName, canName);
    double tiny = 1.E-4;
    double ymax = .1;
    for (int i = i1; i >= i2; --i)
    {
        printf("%i %s %f \n", i, hEffNorm[0][i]->GetName(), hEffNorm[0][i]->Integral());
        printf("%i %s %f \n", i, hEffNorm[1][i]->GetName(), hEffNorm[1][i]->Integral());
        //  hEffNorm[i]->GetYaxis()->SetRangeUser(tiny, ymax);
        // hEffNorm[i]->GetYaxis()->SetRangeUser(tiny, ymax);
        // hEffNorm[i]->GetXaxis()->SetRangeUser(1000, 4000);
        // printf("eff draw %i \n", i);
        if (firstPlot)
        {
            printf("%i %s \n", i, hEffNorm[0][i]->GetName());
            hEffNorm[0][i]->Draw("HIST");
            hEffNorm[1][i]->Draw("HISTSAME");
            firstPlot = false;
        }
        else if (!isBadChannel(i))
        {
            hEffNorm[0][i]->Draw("HISTSAME");
            hEffNorm[1][i]->Draw("HISTSAME");
        }
    }
    can->BuildLegend();
    can->SetLogy();
    return can;
}

void getCurves(TFile *fin = nullptr, int fnum = 0)
{
    printf("get curves from file %s number %i \n", fin->GetName(), fnum);
    TKey *key;
    TIter next(fin->GetListOfKeys());
    int ifile = 0;
    while (TKey *key = (TKey *)next())
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());

        if (!cl->InheritsFrom("TH1D"))
            continue;
        TH1D *h = (TH1D *)key->ReadObj();

        if (TString(h->GetName()).Contains("LightCurve"))
        {
            hcurve[fnum].push_back(h);
            hcurve[fnum][hcurve[fnum].size() - 1]->SetTitle(Form("%s-file%i", h->GetName(), fnum));
        }

        if (TString(h->GetName()).Contains("LightNorm"))
        {
            hnorm[fnum].push_back(h);
            hnorm[fnum][hnorm[fnum].size() - 1]->SetTitle(Form("%s-file%i", h->GetName(), fnum));
        }
        printf("curves from file %s num %i %lu\n", fin->GetName(), fnum, hcurve[fnum].size());
    }
}

// collect trigger bits by file
void getTriggerBits(TFile *fin)
{
    bitFile.clear();
    bitCount.resize(TOTALCODES / 2);

    TIter next(fin->GetListOfKeys());
    TKey *key;
    int ifile = 0;
    while (TKey *key = (TKey *)next())
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());

        if (!cl->InheritsFrom("TH1D"))
            continue;

        TH1D *h = (TH1D *)key->ReadObj();

        if (!TString(h->GetName()).Contains("EventPassFile"))
            continue;

        cout << ifile << " name " << h->GetName() << endl;

        for (int ibin = 0; ibin < h->GetNbinsX(); ++ibin)
            bitCount[ibin].push_back(h->GetBinContent(ibin + 1));

        bitFile.push_back(++ifile);
        printf("files %lu bits %lu  %lu  %lu \n", bitCount[0].size(), bitCount.size(), bitCount[bitCount.size() - 1].size(), bitFile.size());
    }

    // cout << " #bits   " << bitCount.size() << "  # files " << bitCount[0].size() << " " << bitFile.size() << endl;
    for (unsigned ibit = 0; ibit < bitCount.size(); ++ibit)
        printf("\t bit %u files %lu \n", ibit, bitCount[ibit].size());

    // normalize to number of files
    for (unsigned ibit = 0; ibit < bitCount.size(); ++ibit)
    {
        printf("bit %i files %lu \n", ibit, bitCount[ibit].size());
        for (unsigned ifile = 0; ifile < bitCount[ibit].size(); ++ifile)
            bitCount[ibit][ifile] = bitCount[ibit][ifile] / double(bitFile.size());
    }

    cout << " make graphs" << endl;

    // make graph for this bit
    for (unsigned ibit = 0; ibit < bitCount.size(); ++ibit)
    {
        int sum = std::accumulate(bitCount[ibit].begin(), bitCount[ibit].end(), 0.0);
        printf(" sum bit %i # name %s  sum %i  \n", ibit, codeNames[ibit].Data(), sum);
        if (sum > 1000)
        {
            printf(" graph for bit %i # name %s files %lu sum %i  \n", ibit, codeNames[ibit].Data(), bitCount[ibit].size(), sum);
            TGraph *gBit = new TGraph(bitFile.size(), &bitFile[0], &bitCount[ibit][0]);
            gBit->SetName(Form("TrigFailuresBit%i%s", ibit, codeNames[ibit].Data()));
            gBit->SetTitle(Form("TrigFailures for bit %i %s per file ", ibit, codeNames[ibit].Data()));
            gBit->SetMarkerStyle(21);
            gBit->GetYaxis()->SetTitle("falures per file ");
            gBit->GetXaxis()->SetTitle("file number ");
            fout->Append(gBit);
        }
    }
}

void postFit()
{
    /* set bad channels */

    /* set bad channels */
    std::vector<unsigned> badList;
    // badList.push_back(0);
    // badList.push_back(1);
    // badList.push_back(8);
    setBadChannels(badList);

    sdate = currentDate();
    printf(" postMacro on %s \n", sdate.c_str());

    failCodes.resize(FAILBITS);
    failCodes[0] = PASS;
    failCodes[1] = BASEFAIL;
    failCodes[2] = EARLYCUT;
    failCodes[3] = FIRSTTIME;
    failCodes[4] = COSMIC;
    failCodes[5] = GAMMA;
    failCodes[6] = TRIGFAIL;
    failCodes[7] = TRIANGLE;

    bitNames.resize(FAILBITS);
    bitNames[0] = TString("Pass");
    bitNames[1] = TString("Baseline");
    bitNames[2] = TString("Earlycut");
    bitNames[3] = TString("Firsttime");
    bitNames[4] = TString("Cosmic");
    bitNames[5] = TString("Gamma");
    bitNames[6] = TString("Trigger");
    bitNames[7] = TString("Triangle");

    codeNames.resize(TOTALCODES);

    // build trigger bit pattern names
    for (int ic = 0; ic < TOTALCODES; ++ic)
    {
        for (int ibit = 0; ibit < FAILBITS; ++ibit)
            if (ic & failCodes[ibit])
                codeNames[ic] += bitNames[ibit];
    }
    for (unsigned ic = 0; ic < TOTALCODES; ++ic)
        printf("code %i %x name %s \n", ic, ic, codeNames[ic].Data());

    // put in explicit file name and get tag
    // TString fileName("post-10_06_2025-10_06_2025-2392606.root");
    // TString fileName("post-10_16_2025-10_16_2025-2371051.root");
    // TString fileName("post-10_16_2025-10_16_2025-2371051.root");
    TString fileName1;
    TString fileName2;
    fileName1 = TString("post-04_16_2026-04_16_2026-10281297.root");
    fileName2 = TString("post-04_24_2026-04_24_2026-8404446.root");
    TString tag1 = TString(fileName1(fileName1.First("-") + 1, 21));
    TString tag2 = TString(fileName2(fileName2.First("-") + 1, 21));
    cout << " gains from file " << fileName1 << " with date tag " << tag1 << "   " << fileName2 << " with date tag " << tag2 << endl;

    // open file
    fin0 = new TFile(fileName1, "readonly");
    if (fin0->IsZombie())
    {
        printf("no file %s \n", fileName1.Data());
        return;
    }

    // open file
    fin1 = new TFile(fileName2, "readonly");
    if (fin1->IsZombie())
    {
        printf("no file %s \n", fileName2.Data());
        return;
    }

    fout = new TFile(Form("postMacro-%s-%s.root", tag1.Data(), tag2.Data()), "recreate");

    // printf("open file %s date %s \n", fileName.Data(), sdate.c_str());
    // getTriggerBits();

    // geometric eff
    bool geoVersionOld = false;
    setDistanceLevels(geoVersionOld);
    for (unsigned i = 0; i < 12; ++i)
    {
        printf("chan %i distance %f geo eff %.2E\n", i, distanceLevel[getLevel(i)], effGeoFunc(i));
    }

    // light curves
    printf("call getCurves\n");
    hnorm.resize(2);
    hcurve.resize(2);
    hEffNorm.resize(2);
    getCurves(fin0, 0);
    for (unsigned ic = 0; ic < 12; ++ic)
        printf("chan%i %s int %.3E \n", ic, hcurve[0][ic]->GetName(), hcurve[0][ic]->Integral());

    getCurves(fin1, 1);

    for (unsigned ic = 0; ic < 12; ++ic)
        printf("chan%i %s int %.3E \n", ic, hcurve[1][ic]->GetName(), hcurve[1][ic]->Integral());

    int color[12] = {kGreenRedViolet, 2, 3, 4, 5, 6, 7, 8, 9, kTeal, kOrange, kAzure};

    /*
    double ymax = 0;
    double ymin = 1.E9;
    double tiny = 1.E-5;
    for (unsigned i = 0; i < 9; ++i)
    {
    printf("hist %i ymax %f \n", i, hnorm[i]->GetBinContent(hnorm[i]->GetMaximumBin()));
    if (hnorm[i]->GetBinContent(hnorm[i]->GetMaximumBin()) > ymax)
        ymax = hnorm[i]->GetBinContent(hnorm[i]->GetMaximumBin());
    if (hnorm[i]->GetBinContent(hnorm[i]->GetMinimumBin()) < ymin && hnorm[i]->GetBinContent(hnorm[i]->GetMinimumBin()) > tiny)
        ymin = hnorm[i]->GetBinContent(hnorm[i]->GetMinimumBin());
    }
    printf("ymin %f ymax = %f \n ", ymin, ymax);
    */

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    // for (int i = 0; i < hnorm.size(); ++i)
    //     fout->Add(hnorm[0][i]);

    // clone and normalize to geometric efficiency
    hEffNorm[0].resize(12);
    hEffNorm[1].resize(12);

    vector<vector<double>> noiseSum(2);

    for (int fnum = 0; fnum < 2; ++fnum)
        for (unsigned i = 0; i < 12; ++i)
            noiseSum[fnum].push_back(hnorm[fnum][i]->Integral(0, 600));

    // subtact noise
    for (int fnum = 0; fnum < 2; ++fnum)
    {
        for (unsigned i = 0; i < 12; ++i)
        {
            for (int ibin = 0; ibin < hnorm[fnum][i]->GetNbinsX(); ++ibin)
            {
                double val = abs(hnorm[fnum][i]->GetBinContent(ibin) - noiseSum[fnum][i]);
                // hnorm[fnum][i]->SetBinContent(ibin, val);
            }
        }
    }

    // only fit first fiile
    /*   0 amp1_init - initial amplitude 1 guess (default 1.0)
     *   1 tau1_init - initial decay constant 1 guess (default 100.0)
     *   2 amp2_init - initial amplitude 2 guess (default 1.0)
     *   3 tau2_init - initial decay constant 2 guess (default 1000.0)
     *   4 const_init - initial constant background guess (default 0.0)
     */
    vector<vector<double>> amplitude(2);
    vector<vector<double>> tau(2);
    vector<vector<double>> tauErr(2);
    vector<double> constant;
    int fnum = 0;
    for (unsigned i = 0; i < 12; ++i)
    {
        printf(" eff for chan%i %s int %.3E \n", i, hnorm[fnum][i]->GetName(), hnorm[fnum][i]->Integral());
        hEffNorm[fnum][i] = (TH1D *)hnorm[fnum][i]->Clone(Form("effNormChan%ifile %i,", i, fnum));
        for (int ibin = 0; ibin < hnorm[fnum][i]->GetNbinsX(); ++ibin)
        {
            hEffNorm[fnum][i]->SetBinContent(ibin, hnorm[fnum][i]->GetBinContent(ibin) / effGeoFunc(i));
            hEffNorm[fnum][i]->SetBinError(ibin, hnorm[fnum][i]->GetBinError(ibin) / effGeoFunc(i));
        }
        // fout->Append(hEffNorm[i]);
    }

    for (unsigned i = 0; i < 12; ++i)
    {
        hEffNorm[fnum][i]->SetLineColor(color[i]);
        hEffNorm[fnum][i]->SetLineWidth(1);
        hnorm[fnum][i]->SetLineColor(color[i]);
        hnorm[fnum][i]->SetMarkerColor(color[i]);
        hnorm[fnum][i]->SetMarkerStyle(1);
        hnorm[fnum][i]->SetLineWidth(1);
        // hnorm[fnum][i]->Rebin(50);
        fout->Add(hnorm[fnum][i]);
        TF1 *result = fitDoubleExp(hnorm[fnum][i], 1600, 15000);
        // TF1 *result = fitExpPlus(hnorm[fnum][i], 1600, 15000);
        amplitude[0].push_back(result->GetParameter(0));
        tau[0].push_back(result->GetParameter(1));
        tauErr[0].push_back(result->GetParError(1));
        amplitude[1].push_back(result->GetParameter(2));
        tau[1].push_back(result->GetParameter(3));
        tauErr[1].push_back(result->GetParError(3));
    }
    /*
    for (int fnum = 0; fnum < 2; ++fnum)
        for (unsigned i = 0; i < 12; ++i)
            // noise per bin
            noiseSum[fnum][i] /= double(600.);
            */

    TCanvas *can = new TCanvas("fit", "Comparison", 800, 600);
    for (unsigned i = 0; i < 12; ++i)
    {
        can = new TCanvas(Form("can_%d_%d", fnum, i), Form("Channel %d - File %d", i, fnum), 800, 600);
        gStyle->SetOptStat(1);
        gStyle->SetOptFit(1);
        hnorm[fnum][i]->Draw("hist");
        hnorm[fnum][i]->Draw("funcsame");
        can->SetLogy(true);
        can->Print(Form("ExpFit2File%dChannel%d.pdf", fnum, i));
    }

    for (unsigned i = 0; i < 12; ++i)
    {
        printf(" Channel %d: Amplitude1  = %.4E, Tau1  = %.2f +/- %.2f, Amplitude2  = %.4E, Tau2  = %.2f +/- %.2f \n",

               i, amplitude[0][i], tau[0][i], tauErr[0][i],
               amplitude[1][i], tau[1][i], tauErr[1][i]);
    }

    // fout->ls();
    fout->Write();
    return;

    /* make canvases */
    TCanvas *cChan[12];

    for (int i = 0; i < 12; ++i)
        cChan[i] = makeCanCompare(i, TString(Form("channel%i-%s-%s", i, tag1.Data(), tag2.Data())));

    return;

    TCanvas *canCutNormLevelAll = makeCanCutNorm(11, 0, TString("CutNormedLevelAll"));
    TCanvas *canEffNormLevelTrig = makeCanEffNorm(11, 9, TString("EffNormedLevelTrig"));
    TCanvas *canEffNormLevel0 = makeCanEffNorm(2, 0, TString("EffNormedLevel0"));
    TCanvas *canEffNormLevel1 = makeCanEffNorm(5, 3, TString("EffNormedLevel1"));
    TCanvas *canEffNormLevel2 = makeCanEffNorm(8, 6, TString("EffNormedLevel2"));
    TCanvas *canEffNormLevelAll = makeCanEffNorm(11, 0, TString("EffNormedAll"));
}