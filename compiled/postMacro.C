
/*
look at  postAna root file
*/
#include "TReadGains.hxx"
#include "distanceLevels.hh"
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
std::vector<double> rayleighAtten;
std::vector<TString> codeNames;
std::vector<int> failCodes;
std::vector<TString> bitNames;
std::vector<std::vector<double>> bitCount; // [bit][file]
std::vector<double> bitFile;               // file number
std::vector<TH1D *> hnorm;
std::vector<TH1D *> hcurve;
std::vector<TH1D *> hEffNorm;
TH1D *hPeak;
TReadGains *readGains;
std::vector<double> relativeEff;

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
        // hnorm[i]->GetXaxis()->SetRangeUser(1000, 4000);
        // printf("eff draw %i \n", i);
        if (firstPlot)
        {
            printf("%i %s \n", i, hnorm[i]->GetName());
            hnorm[i]->Draw("HISTLINE ");
            firstPlot = false;
        }
        else if (!isBadChannel(i))
            hnorm[i]->Draw("HISTLINESAME");
        else
            printf("skip bad channel %i %s \n", i, hnorm[i]->GetName());
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
        hEffNorm[i]->Rebin(75 * 2);
        printf("%i %s %f \n", i, hEffNorm[i]->GetName(), hEffNorm[i]->Integral());
        //  hEffNorm[i]->GetYaxis()->SetRangeUser(tiny, ymax);
        // hEffNorm[i]->GetYaxis()->SetRangeUser(tiny, ymax);
        // hEffNorm[i]->GetXaxis()->SetRangeUser(1000, 4000);
        // printf("eff draw %i \n", i);
        if (firstPlot)
        {
            printf("%i %s \n", i, hEffNorm[i]->GetName());
            hEffNorm[i]->Draw("HISTLINE");
            firstPlot = false;
        }
        else if (!isBadChannel(i))
            hEffNorm[i]->Draw("HISTLINESAME");
    }
    can->BuildLegend();
    can->SetLogy();
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

        if (TString(h->GetName()).Contains("LightCurve"))
            hcurve.push_back(h);

        if (TString(h->GetName()).Contains("LightNorm"))
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

void postMacro(TString fileName = "post-btbSim-06_17_2026-06_17_2026-822.root")
//"post-04_16_2026-04_16_2026-10281297.root")
{
    bool isSimulation = false;
    if (fileName.Contains("btb"))
        isSimulation = true;
    relativeEff.resize(13);
    for (int ichan = 0; ichan < relativeEff.size(); ++ichan)
        relativeEff[ichan] = 1.;
    // set relative efficiencies

    if (!isSimulation)
    {
        relativeEff[1] = 1. - 0.0656607;
        relativeEff[2] = 1. - 0.0180333;
        relativeEff[3] = 1. - 0.00212032;
        relativeEff[4] = 1. + 0.024446;
        relativeEff[5] = 1. - 0.135226;
        relativeEff[6] = 1. + 0.102685;
        relativeEff[7] = 1. + 0.093910;
    }

    /* set bad channels */
    if (fileName.Sizeof() == 0)
    {
        printf("ERROR: fileName argument is empty. Usage: postMacro(\"post-XX_XX_XXXX-XX_XX_XXXX-XXXXXXX.root\")\n");
        return;
    }

    readGains = new TReadGains();

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

    // fileName is passed as macro argument
    // Example filenames:
    // "post-10_06_2025-10_06_2025-2392606.root"
    // "post-04_24_2026-04_24_2026-2815814.root"
    if (fileName.Length() == 0)
    {
        printf("ERROR: fileName argument is empty. Usage: postMacro(\"post-XX_XX_XXXX-XX_XX_XXXX-XXXXXXX.root\")\n");
        return;
    }
    TString tag = TString(fileName(fileName.First("-") + 1, 21));
    cout << " gains from file " << fileName << " with date tag " << tag << endl;

    // open file
    fin = new TFile(fileName, "readonly");
    if (fin->IsZombie())
    {
        printf("no file %s \n", fileName.Data());
        return;
    }
    fout = new TFile(Form("postMacro-%s.root", tag.Data()), "recreate");

    printf("open file %s date %s \n", fileName.Data(), sdate.c_str());
    getTriggerBits();

    // geometric eff
    bool geoVersionOld = false;
    setDistanceLevels(geoVersionOld);

    // geometric eff
    rayleighAtten.resize(13);
    for (unsigned i = 0; i < 13; ++i)
        rayleighAtten[i] = 1.;
    if (!isSimulation)
        for (unsigned i = 0; i < 13; ++i)
        {
            rayleighAtten[i] = TMath::Exp(distanceLevel[getLevel(i)] / 66.);
            printf("chan %i distance %f geo eff %.2E rayleigh %.2E\n", i, distanceLevel[getLevel(i)], effGeoFunc(i), rayleighAtten[i]);
        }

    // light curves
    printf("call getCurves\n");
    getCurves();
    for (unsigned ic = 0; ic < hcurve.size(); ++ic)
        printf("chan%i %s int %.3E \n", ic, hcurve[ic]->GetName(), hcurve[ic]->Integral());

    int color[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, kTeal, kOrange, kAzure};

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
    for (int i = 0; i < hnorm.size(); ++i)
        fout->Add(hnorm[i]);

    // clone and normalize to geometric efficiency
    hEffNorm.resize(12);
    for (unsigned i = 0; i < 12; ++i)
    {
        printf("chan%i %s int %.3E \n", i, hnorm[i]->GetName(), hnorm[i]->Integral());
        hEffNorm[i] = (TH1D *)hnorm[i]->Clone(Form("effNormChan%i", i));
        for (int ibin = 0; ibin < hnorm[i]->GetNbinsX(); ++ibin)
        // correct for geometric efficiency and rayleigh
        {
            hEffNorm[i]->SetBinContent(ibin, hnorm[i]->GetBinContent(ibin) / relativeEff[i] / effGeoFunc(i) * rayleighAtten[i]);
            hEffNorm[i]->SetBinError(ibin, hnorm[i]->GetBinError(ibin) / relativeEff[i] / effGeoFunc(i) * rayleighAtten[i]);
        }
        // fout->Append(hEffNorm[i]);
    }

    // fout->ls();
    fout->Write();

    for (unsigned i = 0; i < 12; ++i)
    {
        hEffNorm[i]->SetLineColor(color[i]);
        hEffNorm[i]->SetLineWidth(1);
        hnorm[i]->SetLineColor(color[i]);
        hnorm[i]->SetLineWidth(1);
        hnorm[i]->Rebin(75 * 2);
    }

    /* make canvases */
    TCanvas *canCutNormLevelAll = makeCanCutNorm(11, 0, TString("CutNormedLevelAll"));
    TCanvas *canEffNormLevelTrig = makeCanEffNorm(11, 9, TString("EffNormedLevelTrig"));
    TCanvas *canEffNormLevel0 = makeCanEffNorm(2, 0, TString("EffNormedLevel0"));
    TCanvas *canEffNormLevel1 = makeCanEffNorm(5, 3, TString("EffNormedLevel1"));
    TCanvas *canEffNormLevel2 = makeCanEffNorm(8, 6, TString("EffNormedLevel2"));
    TCanvas *canEffNormLevelAll = makeCanEffNorm(11, 0, TString("EffNormedAll"));

    // graph of relative peaks
    hPeak = new TH1D("Peak", "Peak values by channel", 13, 0, 13);
    hPeak->GetYaxis()->SetTitle("singlet peak value (SPE) ");
    hPeak->GetXaxis()->SetTitle("channel number (+1) ");
    double peakValue[12];
    double relativeValue[12];
    double xchan[12];
    double peakAve = 0;
    for (unsigned i = 0; i < 12; ++i)
    {
        xchan[i] = double(i);
        peakValue[i] = hEffNorm[i]->GetBinContent(hEffNorm[i]->GetMaximumBin());
        hPeak->SetBinContent(i + 1, peakValue[i]);
        if (i > 0 && i < 8)
            peakAve += peakValue[i];
    }
    peakAve /= 7.;
    printf("ave %f \n", peakAve);

    for (unsigned i = 0; i < 9; ++i)
    {
        relativeValue[i] = (peakValue[i] - peakAve) / peakAve;
        printf("chan %i peak %f relative %f (%f) \n", i, peakValue[i], relativeValue[i], readGains->relativeGain[i] - 1.);
    }

    TGraph *gPeak = new TGraph(9, &xchan[0], &peakValue[0]);
    TGraph *gRelative = new TGraph(9, &xchan[0], &relativeValue[0]);
    gPeak->SetMarkerStyle(21);
    gRelative->SetMarkerStyle(21);

    TCanvas *cpeak = new TCanvas(Form("peak%s", tag.Data()), Form("peak%s", tag.Data()));
    gPeak->GetHistogram()->GetXaxis()->SetTitle("channel");
    gPeak->GetHistogram()->GetYaxis()->SetTitle("peak value");
    gPeak->Draw("ap");

    TCanvas *crelative = new TCanvas(Form("relative%s", tag.Data()), Form("relative%s", tag.Data()));
    gRelative->GetHistogram()->GetXaxis()->SetTitle("channel");
    gRelative->GetHistogram()->GetYaxis()->SetTitle("relative efficiency");
    gRelative->Draw("ap");
    gRelative->Print("all");

    //
}