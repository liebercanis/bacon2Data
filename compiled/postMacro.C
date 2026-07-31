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
    MULTI = 0x80,
    TOTALCODES = 2 * MULTI // must match postAna.cc
};

enum
{
    FAILBITS = 9
};

std::string sdate;
TFile *fin;
TFile *fout;
bool isSimulation;

std::vector<TString> codeNames;
std::vector<int> failCodes;
std::vector<TString> bitNames;
std::vector<std::vector<double>> bitCount; // [bit][file]
std::vector<double> bitFile;               // file number
std::vector<TH1D *> hnorm;
std::vector<TH1D *> hcurve;
std::vector<TH1D *> hGeoNorm;
std::vector<TH1D *> hEffNorm;
std::vector<TH1D *> hRayleighNorm;
TH1D *hPeak;
TReadGains *readGains;
double rayleighLength = 99.;
double ppm = 0.15;

std::vector<double> absorptionFactor;
std::vector<double> rayleighAtten;
std::vector<double> relativeEff;

void getOtherEffCorrections(bool isSimulation = false)
{
    bool doCorrections = true;
    rayleighAtten.resize(13);
    for (unsigned i = 0; i < 13; ++i)
        rayleighAtten[i] = 1.;
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

    relativeEff.resize(13);
    for (int ichan = 0; ichan < relativeEff.size(); ++ichan)
        relativeEff[ichan] = 1.;
    // set relative efficiencies
    if (!isSimulation && doCorrections)
    {
        for (int ichan = 0; ichan < relativeEff.size(); ++ichan)
            relativeEff[ichan] = readGains->relativeEff[ichan];
    }

    printf("absorption at %.3f\n", ppm);
    for (int i = 0; i < 13; ++i)
    {
        double dist = distanceLevel[getLevel(i)];
        double abs = Absorbtion(ppm, dist);
        absorptionFactor.push_back(1 - abs);
    }
}

void makeEffNorm(int i)
{
    double baseline = hnorm[i]->Integral(0, 600) / 600.; // sum bin range
    double corr = 1. / relativeEff[i] / effGeoFunc(i) / rayleighAtten[i] / absorptionFactor[i];
    printf("chan%i  integral %.3E base %.3E geo %.3E abs %.3E corr %.3E \n", i, hnorm[i]->Integral(), baseline, effGeoFunc(i), absorptionFactor[i], corr);
    hGeoNorm[i] = (TH1D *)hnorm[i]->Clone(Form("GeoNormChan%i", i));
    hRayleighNorm[i] = (TH1D *)hnorm[i]->Clone(Form("RayleighNormChan%i", i));
    hEffNorm[i] = (TH1D *)hnorm[i]->Clone(Form("effNormChan%i", i));
    for (int ibin = 0; ibin < hnorm[i]->GetNbinsX(); ++ibin)
    // correct for geometric efficiency and rayleigh
    {
        //      *absorptionFactor[i];
        hGeoNorm[i]->SetBinContent(ibin, (hnorm[i]->GetBinContent(ibin) - baseline) / effGeoFunc(i));
        hRayleighNorm[i]->SetBinContent(ibin, (hnorm[i]->GetBinContent(ibin) - baseline) / effGeoFunc(i) / rayleighAtten[i]);
        hEffNorm[i]->SetBinContent(ibin, (hnorm[i]->GetBinContent(ibin) - baseline) * corr);
        hEffNorm[i]->SetBinError(ibin, hnorm[i]->GetBinError(ibin) * corr);
    }
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
        // hEffNorm[i]->Rebin(50);
        printf("%i %s %f \n", i, hEffNorm[i]->GetName(), hEffNorm[i]->Integral());
        //  hEffNorm[i]->GetYaxis()->SetRangeUser(tiny, ymax);
        // hEffNorm[i]->GetYaxis()->SetRangeUser(tiny, ymax);
        // hEffNorm[i]->GetXaxis()->SetRangeUser(1300, 1450);
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
    while (TKey *key = (TKey *)next())
    {
        // skip earlier cycles — GetKey returns the highest-cycle key for this name
        if (fin->GetKey(key->GetName()) != key)
            continue;

        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1D"))
            continue;

        TH1D *h = (TH1D *)key->ReadObj();
        if (!h)
            continue;

        TString name(h->GetName());
        if (name.Contains("LightCurve"))
            hcurve.push_back(h);
        if (name.Contains("LightNorm"))
            hnorm.push_back(h);
    }
    printf("getCurves: %lu LightCurve  %lu LightNorm histograms\n", hcurve.size(), hnorm.size());
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
        TString name(key->GetName());
        if (seen.count(name))
            continue;
        seen.insert(name);

        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1D"))
            continue;

        if (!name.Contains("EventPassFile"))
            continue;

        TH1D *h = (TH1D *)fin->Get(name);
        if (!h)
            continue;

        printf("getTriggerBits: %s  nbins=%i\n", name.Data(), h->GetNbinsX());
        for (int ibin = 0; ibin < h->GetNbinsX() && ibin < (int)bitCount.size(); ++ibin)
            bitCount[ibin].push_back(h->GetBinContent(ibin + 1));

        bitFile.push_back(++ifile);
    }

    // cout << " #bits   " << bitCount.size() << "  # files " << bitCount[0].size() << " " << bitFile.size() << endl;
    // or (unsigned ibit = 0; ibit < bitCount.size(); ++ibit)
    //   printf("\t bit %u files %lu \n", ibit, bitCount[ibit].size());

    // normalize to number of files
    for (unsigned ibit = 0; ibit < bitCount.size(); ++ibit)
    {
        // printf("bit %i files %lu \n", ibit, bitCount[ibit].size());
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

void postMacro(TString fileName = "post-04_16_2026-04_16_2026-10281297.root")
{
    isSimulation = false;
    if (fileName.Contains("btb"))
        isSimulation = true;

    /* set bad channels */
    if (fileName.Sizeof() == 0)
    {
        printf("ERROR: fileName argument is empty. Usage: postMacro(\"post-XX_XX_XXXX-XX_XX_XXXX-XXXXXXX.root\")\n");
        return;
    }

    readGains = new TReadGains();

    // get absorption
    setupModelAllFit();

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
    failCodes[8] = MULTI;

    bitNames.resize(FAILBITS);
    bitNames[0] = TString("Pass");
    bitNames[1] = TString("Baseline");
    bitNames[2] = TString("Earlycut");
    bitNames[3] = TString("Firsttime");
    bitNames[4] = TString("Cosmic");
    bitNames[5] = TString("Gamma");
    bitNames[6] = TString("Trigfail");
    bitNames[7] = TString("Triangle");
    bitNames[8] = TString("Multi");

    codeNames.resize(TOTALCODES);

    // build trigger bit pattern names
    for (int ic = 0; ic < TOTALCODES; ++ic)
    {
        for (int ibit = 0; ibit < FAILBITS; ++ibit)
            if (ic & failCodes[ibit])
                codeNames[ic] += bitNames[ibit];
    }
    // for (unsigned ic = 0; ic < TOTALCODES; ++ic)
    //     printf("code %i %x name %s \n", ic, ic, codeNames[ic].Data());

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

    // light curves
    printf("call getCurves\n");
    getCurves();
    for (unsigned ic = 0; ic < hnorm.size(); ++ic)
        printf("chan%i %s int %.3E \n", ic, hnorm[ic]->GetName(), hnorm[ic]->Integral());

    int color[13] = {1, 2, 3, 4, 5, 6, 7, 8, 9, kTeal, kOrange, kAzure, kBlack};

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
    getOtherEffCorrections(isSimulation);
    for (unsigned i = 0; i < 13; ++i)
        printf("chan %u Rayleigh %.3E relative eff %.3E  \n", i, rayleighAtten[i], relativeEff[i]);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    for (int i = 0; i < hnorm.size(); ++i)
        fout->Add(hnorm[i]);

    // clone and normalize to geometric efficiency
    hGeoNorm.resize(13);
    hRayleighNorm.resize(13);
    hEffNorm.resize(13);
    for (unsigned i = 0; i < hEffNorm.size(); ++i)
    {
        makeEffNorm(i);
        fout->Append(hGeoNorm[i]);
        fout->Append(hRayleighNorm[i]);
        fout->Append(hEffNorm[i]);
    }

    for (unsigned i = 0; i < hEffNorm.size(); ++i)
    {
        hEffNorm[i]->SetLineColor(color[i]);
        hEffNorm[i]->SetLineWidth(1);
        hnorm[i]->SetLineColor(color[i]);
        hnorm[i]->SetLineWidth(1);
        // hnorm[i]->Rebin(75 * 2);
    }

    /* make canvases */
    // TCanvas *canCutNormLevelAll = makeCanCutNorm(11, 0, TString("CutNormedLevelAll"));
    // TCanvas *canEffNormLevelTrig = makeCanEffNorm(11, 9, TString("EffNormedLevelTrig"));
    TCanvas *canEffNormLevel0 = makeCanEffNorm(2, 0, TString("EffNormedLevel0"));
    TCanvas *canEffNormLevel1 = makeCanEffNorm(5, 3, TString("EffNormedLevel1"));
    TCanvas *canEffNormLevel2 = makeCanEffNorm(8, 6, TString("EffNormedLevel2"));
    TCanvas *canEffNormLevelAll = makeCanEffNorm(11, 0, TString("EffNormedAll"));

    // graph of relative peaks
    TH1D *hPeak = new TH1D("Peak", "Geometric normalized Peak values", 13, 0, 13);
    hPeak->GetXaxis()->SetTitle("channel number (+1) ");
    hPeak->GetYaxis()->SetTitle("geometric normalized singlet peak value (SPE) ");
    hPeak->GetXaxis()->SetTitle("channel number (+1) ");

    TH1D *hPeakEff = new TH1D("PeakFileEff", " correctd Peak values ", 13, 0, 13);
    hPeakEff->GetYaxis()->SetTitle("correctd singlet peak value (SPE) ");
    hPeakEff->GetXaxis()->SetTitle("channel number (+1) ");

    fout->Append(hPeak);
    fout->Append(hPeakEff);
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

    TGraph *gRay = new TGraph(12, &xchan[0], &rayleighAtten[0]);
    TCanvas *cray = new TCanvas(Form("rayleigh%s", tag.Data()), Form("rayleigh%s ", tag.Data()));
    gRay->SetName("gRay");
    gRay->SetMarkerStyle(21);
    gRay->SetMarkerSize(1.);
    cray->SetGrid();
    gRay->SetTitle(Form("absorption  ppm %.4f", ppm));
    gRay->GetHistogram()->GetXaxis()->SetTitle("channel");
    gRay->GetHistogram()->GetYaxis()->SetTitle("absorption factor");
    gRay->Draw("ap");
    // cray->Print(".pdf");
    fout->Append(gRay);
    gRay->Print("all");
    // gPeakRayleigh->Print("all");

    // peaks
    TGraph *gPeak = new TGraph(12, &xchan[0], &peakValue[0]);
    gPeak->SetName("gPeak");
    gPeak->SetMarkerStyle(21);
    gPeak->SetMarkerColor(kBlue);
    TCanvas *cpeak = new TCanvas(Form("peak%s", tag.Data()), Form("peak%s", tag.Data()));
    gPeak->GetHistogram()->GetXaxis()->SetTitle("channel");
    gPeak->GetHistogram()->GetYaxis()->SetTitle("peak value");
    cpeak->SetGrid();
    gPeak->Draw("ap");

    // rayleigh corrected peaks
    TGraph *gPeakRayleigh = new TGraph(12, &xchan[0], &peakRayleighCorr[0]);
    gPeakRayleigh->SetName("gRayleighPeak");
    gPeakRayleigh->SetMarkerStyle(21);
    gPeakRayleigh->SetMarkerColor(kRed);
    TCanvas *cpeakRayleigh = new TCanvas(Form("peakRayleigh%s", tag.Data()), Form("Rayleigh corrected peak%s", tag.Data()));
    cpeakRayleigh->SetGrid();
    gPeakRayleigh->GetHistogram()->GetXaxis()->SetTitle("channel");
    gPeakRayleigh->GetHistogram()->GetYaxis()->SetTitle("Rayleigh corrected peak value");
    gPeakRayleigh->Draw("ap");
    gPeak->Draw("psame");
    fout->Append(gPeakRayleigh);

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

    TGraph *gRelativeEff = new TGraph(12, &xchan[0], &relativeEff[0]);
    gRelativeEff->SetName("relativeEff");
    gRelativeEff->SetMarkerStyle(21);

    TCanvas *crelativeEff = new TCanvas(Form("relativeEff%s-%.3f", tag.Data(), ppm), Form("relativeEff%s %.3f", tag.Data(), ppm));
    gRelativeEff->GetHistogram()->GetXaxis()->SetTitle("channel");
    gRelativeEff->GetHistogram()->GetYaxis()->SetTitle("relative efficiency");
    crelativeEff->SetGrid();
    gRelativeEff->Draw("ap");
    fout->Append(gRelativeEff);
    // gRelativeEff->Print("all");
    for (unsigned ichan = 0; ichan < 13; ++ichan)
        printf("relativeEff[%u]=1.0+%f;\n", ichan, relativeValue[ichan]);

    /*
    TGraph *gAbs = new TGraph(12, &xchan[0], &absorptionFactor[0]);
    TCanvas *cabsorb = new TCanvas(Form("absorption%s-%.4f", tag.Data(), ppm), Form("absorption%s ppm %.4f", tag.Data(), ppm));
    gAbs->SetName("gAbs");
    gAbs->SetMarkerStyle(21);
    gAbs->SetMarkerSize(1.);
    gAbs->SetTitle(Form("absorption  ppm %.4f", ppm));
    gAbs->GetHistogram()->GetXaxis()->SetTitle("channel");
    gAbs->GetHistogram()->GetYaxis()->SetTitle("absorption factor");
    gAbs->Draw("ap");
    cabsorb->Print(".pdf");
    fout->Append(gAbs);
    */

    // fout->ls();:w
    vector<double> normFactor;
    normFactor.resize(13);
    for (unsigned i = 0; i < 13; ++i)
        normFactor[i] = 2.46000e+03; // default normalization factor
    /* with large bad fit region!
normFactor[1] = 1.58E3;
normFactor[2] = 1.378E3;
normFactor[3] = 9.1665E2;
normFactor[4] = 1.0E3;
normFactor[5] = 8.3176E2;
normFactor[6] = 9.531E2;
normFactor[7] = 9.96E2;
*/
    normFactor[1] = 9.15E2;
    normFactor[2] = 7.33E2;
    normFactor[3] = 6.86E2;
    normFactor[4] = 6.62E2;
    normFactor[5] = 6.11E2;
    normFactor[6] = 8.84E2;
    normFactor[7] = 8.848E2;

    double aveNormFactor = 0;
    for (unsigned i = 1; i < 8; ++i)
        aveNormFactor += normFactor[i];
    aveNormFactor /= double(7);

    printf("average normalization factor %.3E \n", aveNormFactor);

    vector<double> relativeNorm;
    relativeNorm.resize(13);
    for (unsigned ichan = 0; ichan < 13; ++ichan)
        relativeNorm[ichan] = 1.0 + (normFactor[ichan] - aveNormFactor) / aveNormFactor;

    TGraph *gNormFactor = new TGraph(12, &xchan[0], &relativeNorm[0]);
    gNormFactor->SetName("normFactor");
    gNormFactor->SetMarkerStyle(21);

    TCanvas *cnormFactor = new TCanvas(Form("relativeNormFactor%s-%.3f", tag.Data(), ppm), Form("relateiveNormFactor%s %.3f", tag.Data(), ppm));
    gNormFactor->GetHistogram()->GetXaxis()->SetTitle("channel");
    gNormFactor->GetHistogram()->GetYaxis()->SetTitle("fit normalization factor");
    cnormFactor->SetGrid();
    gNormFactor->Draw("ap");
    fout->Append(gNormFactor);

    for (unsigned ichan = 0; ichan < 13; ++ichan)
        printf("relativeNorm[%u]=1.0+%f;\n", ichan, relativeNorm[ichan]);

    fout->Write();
    //
}