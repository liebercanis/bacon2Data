
/*
look at  postAna root file
*/

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

std::string sdate;
TFile *fin;
TFile *fout;
std::vector<TString> codeNames;
std::vector<int> failCodes;
std::vector<TString> bitNames;
std::vector<std::vector<double>> bitCount; // [bit][file]
std::vector<double> bitFile;               // file number
std::vector<TH1D *> hnorm;
std::vector<TH1D *> hcurve;
std::vector<TH1D *> hEffNorm;

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

TCanvas *makeCanCutNorm(int i1, int i2, TString canName)
{
    printf(" makeCanCutfNorm %s from %i to %i size %lu \n", canName.Data(), i1, 12, hnorm.size());
    bool firstPlot = true;
    TCanvas *can = new TCanvas(canName, canName);
    for (int i = i1; i >= i2; --i)
    {
        // printf("%i %s \n", i, hnorm[i]->GetName());
        //  hnorm[i]->GetYaxis()->SetRangeUser(tiny, ymax);
        // hnorm[i]->GetYaxis()->SetRangeUser(0.1, 30);
        hnorm[i]->GetXaxis()->SetRangeUser(1000, 4000);
        // printf("eff draw %i \n", i);
        if (firstPlot)
        {
            printf("%i %s \n", i, hnorm[i]->GetName());
            hnorm[i]->Draw("HIST");
            firstPlot = false;
        }
        else if (!isBadChannel(i))
            hnorm[i]->Draw("HISTSAME");
    }
    can->BuildLegend();
    return can;
}

TCanvas *makeCanEffNorm(int i1, int i2, TString canName)
{
    printf(" makeCanEffNorm %s from %i to %i size %lu \n", canName.Data(), i1, 12, hEffNorm.size());
    bool firstPlot = true;
    TCanvas *can = new TCanvas(canName, canName);
    double tiny = 1.E-2;
    double ymax = .1;
    for (int i = i1; i >= i2; --i)
    {
        // printf("%i %s \n", i, hEffNorm[i]->GetName());
        //  hEffNorm[i]->GetYaxis()->SetRangeUser(tiny, ymax);
        hEffNorm[i]->GetYaxis()->SetRangeUser(0.1, 30);
        hEffNorm[i]->GetXaxis()->SetRangeUser(1000, 4000);
        // printf("eff draw %i \n", i);
        if (firstPlot)
        {
            printf("%i %s \n", i, hEffNorm[i]->GetName());
            hEffNorm[i]->Draw("HIST");
            firstPlot = false;
        }
        else if (!isBadChannel(i))
            hEffNorm[i]->Draw("HISTSAME");
    }
    can->BuildLegend();
    return can;
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
            hnorm.push_back(h);
    }
}

// collect trigger bits by file
void getTriggerBits()
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

void postMacro()
{

    /* set bad channels */

    /* set bad channels */
    std::vector<unsigned> badList;
    badList.push_back(0);
    badList.push_back(1);
    badList.push_back(8);
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
    TString fileName("post-anaCRun-btbSimNEW-2026-02-13-100000-7857.root");
    fileName = TString("post-11_19_2025-11_19_2025-1371746.root");
    TString tag = TString(fileName(fileName.First("-") + 1, 21));
    cout << " gains from file " << fileName << " with date tag " << tag << endl;

    // open file
    fin = new TFile(fileName, "readonly");
    if (fin->IsZombie())
    {
        printf("no file %s \n", fileName.Data());
        return;
    }
    fout = new TFile("postMacro.root", "recreate");

    printf("open file %s date %s \n", fileName.Data(), sdate.c_str());
    getTriggerBits();

    // geometric eff
    bool geoVersionOld = false;
    setDistanceLevels(geoVersionOld);
    for (unsigned i = 0; i < 12; ++i)
    {
        printf("chan %i distance %f geo eff %.2E\n", i, distanceLevel[getLevel(i)], effGeoFunc(i));
    }

    // light curves
    getCurves();
    for (unsigned ic = 0; ic < hcurve.size(); ++ic)
        printf("chan%i %s int %.3E \n", ic, hcurve[ic]->GetName(), hcurve[ic]->Integral());

    int color[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, kTeal, kOrange, kAzure};

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

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    // clone and normalize to geometric efficiency
    hEffNorm.resize(12);
    for (unsigned i = 0; i < 12; ++i)
    {
        hEffNorm[i] = (TH1D *)hnorm[i]->Clone(Form("effNormChan%i", i));
        for (int ibin = 1; ibin < hnorm[i]->GetNbinsX(); ++ibin)
        {
            hEffNorm[i]->SetBinContent(ibin, hEffNorm[i]->GetBinContent(ibin) / effGeoFunc(i));
            hEffNorm[i]->SetBinError(ibin, hEffNorm[i]->GetBinError(ibin) / effGeoFunc(i));
        }
        fout->Append(hEffNorm[i]);
    }

    fout->Write();

    for (unsigned i = 0; i < 12; ++i)
    {
        hEffNorm[i]->SetLineColor(color[i]);
        hEffNorm[i]->SetLineWidth(1);
        hnorm[i]->SetLineColor(color[i]);
        hnorm[i]->SetLineWidth(1);
    }

    /* make canvases */
    TCanvas *canCutNormLevelAll = makeCanCutNorm(11, 0, TString("CutNormedLevelAll"));
    TCanvas *canEffNormLevelTrig = makeCanEffNorm(11, 9, TString("EffNormedLevelTrig"));
    TCanvas *canEffNormLevel0 = makeCanEffNorm(2, 0, TString("EffNormedLevel0"));
    TCanvas *canEffNormLevel1 = makeCanEffNorm(5, 3, TString("EffNormedLevel1"));
    TCanvas *canEffNormLevel2 = makeCanEffNorm(8, 6, TString("EffNormedLevel2"));
    TCanvas *canEffNormLevelAll = makeCanEffNorm(11, 0, TString("EffNormedAll"));
}