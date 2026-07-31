/*
look at  postAna root file
*/

// C++ standard library includes
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <numeric>
#include <algorithm>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TString.h>
#include <TObject.h>

// pass bit failures hex
#include "modelAllFit.hh"
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

vector<double> ppmFile;

unsigned ifile;

std::string sdate;
TFile *fin;
TFile *fout;
TGraph *gSingletIntegrals;
TGraph *gLateIntegrals;
std::vector<double> singletIntegral;
std::vector<double> lateIntegral;
std::vector<TString> codeNames;
std::vector<int> failCodes;
std::vector<TString> bitNames;
std::vector<std::vector<double>> bitCount; // [bit][file]
std::vector<double> bitFile;               // file number
std::vector<TH1D *> hnorm;
std::vector<TH1D *> hcurve;
double rayleighLength = 99.;
double ppm = 0.15;

TReadGains *readGains;
std::vector<double> absorptionFactor;
std::vector<double> rayleighAtten;
std::vector<double> relativeEff;
std::vector<TH1D *> hGeoNorm;
std::vector<TH1D *> hEffNorm;
std::vector<TH1D *> hRayleighNorm;
std::vector<TString> fileList;

void getOtherEffCorrections(bool isSimulation = false)
{
    bool doCorrections = true;
    rayleighAtten.resize(13);
    for (unsigned i = 0; i < 13; ++i)
        rayleighAtten[i] = 1.;
    /*
if (!isSimulation)
{
    for (unsigned i = 0; i < 13; ++i)
    {
        if (i == 9 || i == 10 || i == 11) // no correction for trigger
            continue;
        // half backward scattered 1+cos^2
        rayleighAtten[i] = 1. - (1. - TMath::Exp(-distanceLevel[getLevel(i)] / rayleighLength)) / 2.;
        printf("chan %i distance %f geo eff %.2E rayleigh %.2E\n", i, distanceLevel[getLevel(i)], effGeoFunc(i), rayleighAtten[i]);
    }
}
    */

    relativeEff.resize(13);
    for (int ichan = 0; ichan < relativeEff.size(); ++ichan)
        relativeEff[ichan] = 1.;
    // set relative efficiencies
    /*
    if (!isSimulation && doCorrections)
    {
        for (int ichan = 0; ichan < relativeEff.size(); ++ichan)
            relativeEff[ichan] = readGains->relativeEff[ichan];
    }
            */

    absorptionFactor.resize(13);
    for (int ichan = 0; ichan < absorptionFactor.size(); ++ichan)
        absorptionFactor[ichan] = 1.;

    /*
    printf("absorption at %.3f\n", ppm);
    for (int i = 0; i < 13; ++i)
    {
        double dist = distanceLevel[getLevel(i)];
        double abs = Absorbtion(ppm, dist);
        absorptionFactor[i] = 1 - abs;
    }
        */
}

void makeEffNorm(int i)
{
    double baseline = hnorm[i]->Integral(0, 600) / 600.; // sum bin range
    // double corr = 1. / relativeEff[i] / effGeoFunc(i) / rayleighAtten[i] / absorptionFactor[i];
    double corr = 1. / effGeoFunc(i) / readGains->relativeNorm[i];
    printf("chan%i  integral %.3E base %.3E geo %.3E abs %.3E corr %.3E \n", i, hnorm[i]->Integral(), baseline, effGeoFunc(i), absorptionFactor[i], corr);
    hGeoNorm[i] = (TH1D *)hnorm[i]->Clone(Form("GeoNormChan%i", i));
    hRayleighNorm[i] = (TH1D *)hnorm[i]->Clone(Form("RayleighNormChan%i", i));
    hEffNorm[i] = (TH1D *)hnorm[i]->Clone(Form("effNormChan%i", i));
    for (int ibin = 0; ibin < hnorm[i]->GetNbinsX(); ++ibin)
    // correct for geometric efficiency and rayleigh
    {
        double val = hnorm[i]->GetBinContent(ibin) - baseline;
        val = max(val, 0.); // avoid negative values
        //      *absorptionFactor[i];
        hGeoNorm[i]->SetBinContent(ibin, val / effGeoFunc(i));
        hRayleighNorm[i]->SetBinContent(ibin, val / effGeoFunc(i) / rayleighAtten[i]);
        hEffNorm[i]->SetBinContent(ibin, val * corr);
        hEffNorm[i]->SetBinError(ibin, hnorm[i]->GetBinError(ibin) * corr);
    }
}

// Color palette with 12 distinct colors for plotting
int colorPalette[13] = {
    kBlack,      // 0 - black
    kRed - 2,    // 1 - red
    kBlue,       // 2 - blue
    kGreen,      // 3 - green
    kMagenta,    // 4 - magenta
    kCyan,       // 5 - cyan
    kYellow,     // 6 - yellow
    kOrange,     // 7 - orange
    kViolet,     // 8 - violet
    kSpring,     // 9 - spring (green-cyan)
    kRed + 2,    // 10 - teal (blue-green)
    kTeal,       // 11 - azure (light blue)
    kMagenta + 3 // 11 - azure (light blue)
};

// Helper function to get color from palette (wraps around if index >= 12)
int getColor(int index)
{
    return colorPalette[index % 12];
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

TCanvas *canEffNormFile()
{
    printf(" makeCanEffNorm file %u \n", ifile);
    bool firstPlot = true;
    TString canName(Form("canEffNormAllFile%u", ifile));
    TString canTitle(Form("canEffNormAllFile%u", ifile));
    TCanvas *can = new TCanvas(canName, canTitle);
    double tiny = 1.E-4;
    double ymin = 1.E5;
    double ymax = 1.E-4;
    for (int i = hEffNorm.size() - 1; i >= 0; --i)
    {
        // be careful here!
        // hEffNorm[i]->Rebin(50);
        int maxBin = hEffNorm[i]->GetMaximumBin();
        double maxValue = hEffNorm[i]->GetBinContent(maxBin);
        if (maxValue > ymax)
            ymax = maxValue;
        int minBin = hEffNorm[i]->GetMaximumBin();
        double minValue = hEffNorm[i]->GetBinContent(minBin);
        if (minValue < ymin)
            ymin = minValue;
    }
    ymin *= 0.1;
    ymax *= 1.2;
    printf("ymin %f ymax %f \n", ymin, ymax);
    for (int i = hEffNorm.size() - 1; i >= 0; --i)
    {
        int icolor = colorPalette[i];
        hEffNorm[i]->SetLineColor(icolor);
        hEffNorm[i]->SetMarkerColor(icolor);
        hEffNorm[i]->SetLineWidth(1);
        hEffNorm[i]->SetName(Form("effNormChan%iFile%i", i, ifile));
        hEffNorm[i]->SetTitle(Form("effNormChan%iPPM%.3fFile%i", i, ppmFile[ifile], ifile));
        // hEffNorm[i]->GetXaxis()->SetRangeUser(1000, 75000); // ramge im bins
        hEffNorm[i]->GetYaxis()->SetRangeUser(ymin, ymax);
        printf("%i %s %f \n", i, hEffNorm[i]->GetName(), hEffNorm[i]->Integral());
        if (firstPlot)
        {
            printf("%i %s \n", i, hEffNorm[i]->GetName());
            hEffNorm[i]->Draw("HISTILINE");
            firstPlot = false;
        }
        else if (!isBadChannel(i))
            hEffNorm[i]->Draw("HISTLINESAME");
    }

    // Create legend in upper right corner
    TLegend *leg = new TLegend(0.6, 0.7, .99, .99);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    for (unsigned i = 0; i < hEffNorm.size(); ++i)
    {
        leg->AddEntry(hEffNorm[i], Form("chan %u", i), "l");
    }
    leg->Draw();

    can->SetLogy();
    can->SetGrid();
    can->Print(".pdf");
    return can;
}

// Get histogram by name from input file
TH1D *getHistFromFile(TString histName)
{
    if (!fin || fin->IsZombie())
    {
        printf("ERROR: fin is not open or is zombie\n");
        return 0;
    }

    TObject *obj = fin->Get(histName);
    if (!obj)
    {
        printf("ERROR: Histogram '%s' not found in file\n", histName.Data());
        return 0;
    }

    TH1D *hist = (TH1D *)obj;
    if (!hist)
    {
        printf("ERROR: Object '%s' is not a TH1D\n", histName.Data());
        return 0;
    }

    printf("Retrieved histogram '%s' peak bin %i from file %s\n", histName.Data(), hist->GetMaximumBin(), fin->GetName());
    fout->Append(hist);
    return hist;
}

// collect trigger bits by file
void getTriggerBits()
{
    bitFile.clear();
    bitCount.resize(TOTALCODES / 2);

    std::set<TString> seen;
    TIter next(fin->GetListOfKeys());
    int ifile = 0;
    while (TKey *key = (TKey *)next())
    {
        // last key only — skip earlier cycles
        if (fin->GetKey(key->GetName()) != key)
            continue;

        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1D"))
            continue;

        TH1D *h = (TH1D *)fin->Get(key->GetName());
        if (!h || !TString(h->GetName()).Contains("EventPassFile"))
            continue;

        cout << ifile << " name " << h->GetName() << " nbins " << h->GetNbinsX() << endl;

        for (int ibin = 0; ibin < h->GetNbinsX() && ibin < (int)bitCount.size(); ++ibin)
            bitCount[ibin].push_back(h->GetBinContent(ibin + 1));

        bitFile.push_back(++ifile);
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
        // printf(" sum bit %i # name %s  sum %i  \n", ibit, codeNames[ibit].Data(), sum);
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

// read all post files and make summary plots
void postMacroAllFile(unsigned theFile = 0)
{
    readGains = new TReadGains();
    ifile = theFile;
    fout = new TFile(Form("postMacroAllFile%i.root", ifile), "recreate");

    fileList.push_back("post-04_16_2026-04_16_2026-10281297.root");
    fileList.push_back("post-04_24_2026-04_24_2026-10341507.root");
    fileList.push_back("post-04_28_2026-04_28_2026-11344902.root");
    fileList.push_back("post-04_30_2026-04_30_2026-6910781.root");
    fileList.push_back("post-05_01_2026-05_01_2026-10587911.root");
    fileList.push_back("post-05_03_2026-05_03_2026-10135057.root");
    fileList.push_back("post-05_05_2026-05_05_2026-12500971.root");
    fileList.push_back("post-05_07_2026-05_07_2026-10063119.root");
    fileList.push_back("post-05_09_2026-05_09_2026-11006186.root");
    fileList.push_back("post-05_11_2026-05_11_2026-12249461.root");
    fileList.push_back("post-05_14_2026-05_14_2026-10233851.root");
    // fileList.push_back("post-05_18_2026-05_18_2026-7107472.root"); 10 random trigger
    fileList.push_back("post-05_26_2026-05_26_2026-8083870.root");
    fileList.push_back("post-06-09-2026-06-09-2026-9598597.root");

    ppmFile.push_back(0);
    ppmFile.push_back(0.01);
    ppmFile.push_back(0.03);
    ppmFile.push_back(0.05);
    ppmFile.push_back(0.1);
    ppmFile.push_back(0.3);
    ppmFile.push_back(0.5);
    ppmFile.push_back(1.);
    ppmFile.push_back(2.);
    ppmFile.push_back(5.);
    ppmFile.push_back(10.);
    ppmFile.push_back(1.);
    ppmFile.push_back(30.);

    printf("postMacroAll read %lu files \n", fileList.size());
    if (ifile > fileList.size() - 1)
    {
        printf("no such file %u \n", ifile);
    }

    for (unsigned i = 0; i < fileList.size(); ++i)
        printf("file %u %s PPM %f \n", i, fileList[i].Data(), ppmFile[i]);

    // open file
    printf("read file %u %s \n", ifile, fileList[ifile].Data());
    fin = new TFile(fileList[ifile], "readonly");
    if (fin->IsZombie())
    {
        printf("no file %s \n", fileList[ifile].Data());
        return;
    }
    // collect histograms from this file
    for (unsigned ichan = 0; ichan < 13; ++ichan)
    {
        TString histName = Form("LightNormChan%i", ichan);
        hnorm.push_back(getHistFromFile(histName));
        // fin->Close();
    }
    TString tag = TString(fileList[ifile](fileList[ifile].First("-") + 1, 21));

    printf("read %lu histograms \n", hnorm.size());
    for (unsigned ichan = 0; ichan < hnorm.size(); ++ichan)
    {
        if (!hnorm[ichan])
        {
            printf("WARNING: hnorm[%u] is null\n", ichan);
            continue;
        }
        printf("hist %s \n", hnorm[ichan]->GetName());
    }

    // geometric eff
    bool geoVersionOld = false;
    setDistanceLevels(geoVersionOld);
    for (unsigned i = 0; i < 12; ++i)
    {
        printf("chan %i distance %f geo eff %.2E\n", i, distanceLevel[getLevel(i)], effGeoFunc(i));
    }

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
    // for (unsigned ic = 0; ic < TOTALCODES; ++ic)
    //    printf("code %i %x name %s \n", ic, ic, codeNames[ic].Data());

    getOtherEffCorrections();
    for (unsigned i = 0; i < 13; ++i)
        printf("chan %u Rayleigh %.3E relative eff %.3E   \n", i, rayleighAtten[i], relativeEff[i]);
    printf("file %u makeEffHistos \n", ifile);

    // clone and normalize to geometric efficiency
    // cd into fout so Clone() places new histograms in fout, not fin
    fout->cd();
    hGeoNorm.resize(13);
    hRayleighNorm.resize(13);
    hEffNorm.resize(13);
    for (unsigned i = 0; i < hEffNorm.size(); ++i)
    {
        if (!hnorm[i])
        {
            printf("WARNING: hnorm[%u] null, skipping makeEffNorm\n", i);
            continue;
        }
        makeEffNorm(i);
        fout->Append(hGeoNorm[i]);
        fout->Append(hRayleighNorm[i]);
        fout->Append(hEffNorm[i]);
    }
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    printf("file %i make canvas \n", ifile);
    /* make canvases */
    TCanvas *canEffNormCanvas = canEffNormFile();
    // peak plot
    // graph of relative peaks
    TH1D *hPeak = new TH1D(Form("PeakFile%i", ifile), Form("Geometric normalized Peak values by channel file %i", ifile), 13, 0, 13);
    hPeak->GetYaxis()->SetTitle("geometric corrected singlet peak value (SPE) ");
    hPeak->GetXaxis()->SetTitle("channel number (+1) ");
    TH1D *hPeakEff = new TH1D(Form("PeakFileEff%i", ifile), Form("Corrected Peak values by channel file %i", ifile), 13, 0, 13);
    hPeak->GetYaxis()->SetTitle("corrected singlet peak value (SPE) ");
    hPeak->GetXaxis()->SetTitle("channel number (+1) ");
    hPeakEff->GetYaxis()->SetTitle("singlet peak value (SPE) ");
    hPeakEff->GetXaxis()->SetTitle("channel number (+1) ");
    fout->Append(hPeak);
    fout->Append(hPeakEff);

    TCanvas *canPeak = new TCanvas(Form("SingletPeakFile%i", ifile), Form("SintletPeakFile%i", ifile));
    hPeak->Draw();
    // canPeak->SetLogy();
    canPeak->Print(".pdf");

    TCanvas *canPeakEff = new TCanvas(Form("SingletEffPeakFile%i", ifile), Form("SintletEffPeakFile%i", ifile));
    hPeakEff->Draw();
    canPeakEff->Print(".pdf");

    double peakValue[13];
    double peakRayleighCorr[13];
    double peakValueCorr[13];
    double relativeValue[13];
    double xchan[13];
    double peakAve = 0;
    for (unsigned i = 0; i < hEffNorm.size(); ++i)
    {
        xchan[i] = double(i);
        // double baseline = hnorm[i]->Integral(0, 600) / 600.; // sum bin range
        // peakValue[i] = (hnorm[i]->GetBinContent(hnorm[i]->GetMaximumBin() - baseline) / effGeoFunc(i));
        peakValue[i] = hGeoNorm[i]->GetBinContent(hGeoNorm[i]->GetMaximumBin());
        peakRayleighCorr[i] = hRayleighNorm[i]->GetBinContent(hRayleighNorm[i]->GetMaximumBin());
        peakValueCorr[i] = hEffNorm[i]->GetBinContent(hEffNorm[i]->GetMaximumBin());
        hPeak->SetBinContent(i + 1, peakValue[i]);
        hPeakEff->SetBinContent(i + 1, peakValueCorr[i]);
        printf("MESSAGE chan %u max bin %i geo %.3E peakVal %.3E peakRayleigh %.3E peakEffVal %.3E \n", i, hnorm[i]->GetMaximumBin(), effGeoFunc(i), peakValue[i], peakRayleighCorr[i], peakValueCorr[i]);
        // peak ave is of non trigger
        if (i > 0 && i < 8)
            peakAve += peakValueCorr[i];
    }
    peakAve /= 7.;
    printf("ave %f \n", peakAve);

    printf("rayleigh scattering length is %.3f\n", rayleighLength);

    for (unsigned i = 0; i < 13; ++i)
    {
        relativeValue[i] = (peakValueCorr[i] - peakAve) / peakAve;
        printf("RELATIVEVALUES chan %i peak %.3E relative %E (%.3E) distance %.3f rayleigh %.3f  peak %.3E corr %.3E; \n", i, peakValue[i], relativeValue[i], readGains->relativeEff[i] - 1., distanceLevel[getLevel(i)], rayleighAtten[i], peakValue[i], peakValueCorr[i]);
    }
    // fill integral arrays
    singletIntegral.resize(NCHAN);
    lateIntegral.resize(NCHAN);
    for (int i = 0; i < NCHAN; ++i)
    {
        // integral is over bins
        int earlyBin = 1400 / 2;
        int lateBin = 1460 / 2;
        // Integral(Int_t binx1, Int_t binx2, Option_t *option="") const
        singletIntegral[i] = hEffNorm[i]->Integral(earlyBin, lateBin);
        lateIntegral[i] = hEffNorm[i]->Integral(lateBin, 3000 / 2);
        printf("line546 chan %i singlet %E \n", i, singletIntegral[i]);
        //  singletIntegral[i] = hGeoNorm[i]->Integral(0., 1400.);
        //  lateIntegral[i] = hGeoNorm[i]->Integral(1400., 3500.);

        /*
        for (int j = 0; j < hEffNorm[i]->GetNbinsX(); ++j)
        {
            double val = hEffNorm[i]->GetBinContent(j);
            if (j < 1400 / 2)
                singletIntegral[i] += val;
            else if (j > 1400 / 2 && j < 3000 / 2)
                lateIntegral[i] += val;
        }
        */
    }

    gSingletIntegrals = new TGraph(NCHAN, &xchan[0], &singletIntegral[0]);
    gLateIntegrals = new TGraph(NCHAN, &xchan[0], &lateIntegral[0]);

    gSingletIntegrals->SetName(Form("gSingletIntegralsFile%i", theFile));
    gSingletIntegrals->SetTitle(Form("gSingletIntegralsFile%i", theFile));
    fout->Append(gSingletIntegrals);

    gLateIntegrals->SetName(Form("gLateIntegralsFile%i", theFile));
    gLateIntegrals->SetTitle(Form("gLateIntegralsFile%i", theFile));
    fout->Append(gLateIntegrals);

    // to fix the scale for the double Draw, use TMultiGraph

    TCanvas *cintegral = new TCanvas(Form("integrals%s", tag.Data()), Form("integral%s", tag.Data()));
    gSingletIntegrals->GetHistogram()->GetXaxis()->SetTitle("channel");
    gSingletIntegrals->GetHistogram()->GetYaxis()->SetTitle("integral value");
    gLateIntegrals->GetHistogram()->GetXaxis()->SetTitle("channel");
    gLateIntegrals->GetHistogram()->GetYaxis()->SetTitle("integral value");
    gSingletIntegrals->SetMarkerStyle(21);
    gLateIntegrals->SetMarkerStyle(22);
    gLateIntegrals->Draw("ap");
    gSingletIntegrals->Draw("psame");
    cintegral->SetLogx();
    cintegral->BuildLegend();

    TGraph *gPeak = new TGraph(12, &xchan[0], &peakValue[0]);
    gPeak->SetName("gPeak");
    TGraph *gRelative = new TGraph(12, &xchan[0], &relativeValue[0]);
    gRelative->SetName("gRelative");
    gPeak->SetMarkerStyle(21);
    gRelative->SetMarkerStyle(21);

    TCanvas *cpeak = new TCanvas(Form("peak%s", tag.Data()), Form("peak%s", tag.Data()));
    gPeak->GetHistogram()->GetXaxis()->SetTitle("channel");
    gPeak->GetHistogram()->GetYaxis()->SetTitle("peak value");
    gPeak->Draw("ap");
    fout->Append(gPeak);

    TCanvas *crelative = new TCanvas(Form("relative%s", tag.Data()), Form("relative%s", tag.Data()));
    gRelative->GetHistogram()->GetXaxis()->SetTitle("channel");
    gRelative->GetHistogram()->GetYaxis()->SetTitle("relative efficiency");
    gRelative->Draw("ap");
    gRelative->Print("all");
    fout->Append(gRelative);

    // corrected peaks
    TGraph *gPeakCorr = new TGraph(12, &xchan[0], &peakValueCorr[0]);
    gPeakCorr->SetName("gPeakCorr");
    gPeakCorr->SetMarkerStyle(21);
    TCanvas *cpeakCorr = new TCanvas(Form("peakCorr%s", tag.Data()), Form("corrected peak%s", tag.Data()));
    cpeakCorr->SetGrid();
    gPeakCorr->GetHistogram()->GetXaxis()->SetTitle("channel");
    gPeakCorr->GetHistogram()->GetYaxis()->SetTitle("corrected peak value");
    gPeakCorr->Draw("ap");
    gPeakCorr->Print("all");
    fout->Append(gPeakCorr);

    // gRelativeEff->Print("all");
    for (unsigned ichan = 0; ichan < 13; ++ichan)
        printf("relativeEff[%u]=1.+%f;\n", ichan, relativeValue[ichan]);

    for (unsigned ichan = 0; ichan < 13; ++ichan)
        printf("chan %i singlet %f late %f \n", ichan, singletIntegral[ichan], lateIntegral[ichan]);
    printf("MESSAGE end of file %s \n", fout->GetName());

    // fout->ls();
    fout->Write();
    // Disown all objects from fout before Close().  fout->Write() has already
    // persisted everything; if we leave objects in fout's list, Close() will
    // delete them a second time after ROOT's canvas cleanup already did so
    // at macro exit, producing the segfault.  Clear() removes list entries
    // without deleting the objects; the canvases become the sole owners.
    fout->GetList()->Clear();
    fin->Close();
    fout->Close();
}
