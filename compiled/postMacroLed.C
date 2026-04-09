/*
look at  postAna root file
*/
#include "modelAllFit.hh"
enum
{
    NONSUMCHANNELS = 9
};

std::string sdate;
TFile *fin;
TFile *fout;
std::vector<TH1D *> hLateSumChan;

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

void getCurves()
{

    TDirectory *ledDir = fin->GetDirectory("ledDir");
    TIter next(ledDir->GetListOfKeys());
    TKey *key;
    int ifile = 0;
    while (TKey *key = (TKey *)next())
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());

        if (!cl->InheritsFrom("TH1D"))
            continue;
        TH1D *h = (TH1D *)key->ReadObj();

        if (TString(h->GetName()).Contains("LateSumChan"))
        {
            hLateSumChan.push_back(h);
            fout->Append(h);
        }
    }
}

void postMacroLed()
{

    // set file name and tag
    TString fileName = TString("post-02_26_2026-02_26_2026-265584.root");
    TString tag = TString(fileName(fileName.First("-") + 1, 21));
    cout << " led data file " << fileName << " with date tag " << tag << endl;

    // open file
    fin = new TFile(fileName, "readonly");
    if (fin->IsZombie())
    {
        printf("no file %s \n", fileName.Data());
        return;
    }
    fout = new TFile("postMacro.root", "recreate");

    // geometric eff
    bool geoVersionOld = false;
    setDistanceLevels(geoVersionOld);
    for (unsigned i = 0; i < 12; ++i)
    {
        printf("chan %i distance %f geo eff %.2E\n", i, distanceLevel[getLevel(i)], effGeoFunc(i));
    }
    getCurves();
    printf("Number of late sum channels: %lu\n", hLateSumChan.size());
    TH1D *hIntChannel = new TH1D("hIntChannnel", "int by channel", 9, 0, 9);
    TH1D *hGeoIntChannel = new TH1D("hGeoIntChannnel", "int by channel geo normalized", 9, 0, 9);
    hIntChannel->GetXaxis()->SetTitle("channel");
    hIntChannel->GetYaxis()->SetTitle("late sum average [SPE]");
    hGeoIntChannel->GetXaxis()->SetTitle("channel");
    hGeoIntChannel->GetYaxis()->SetTitle("late sum average [SPE] / geometric efficiency");

    int color[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, kTeal, kOrange, kAzure};
    for (unsigned i = 0; i < 9; ++i)
    {
        hLateSumChan[i]->SetLineColor(color[i]);
        hLateSumChan[i]->SetLineWidth(1);
        hLateSumChan[i]->SetStats(kFALSE);
    }

    TCanvas *can = new TCanvas("canLateSum", "canLateSum");
    can->Divide(3, 3);
    for (unsigned i = 0; i < 9; ++i)
    {
        can->cd(i + 1);
        hLateSumChan[i]->Draw("hist");
    }

    TCanvas *can0 = new TCanvas("canLateSumLevel0", "canLateSumLevel0");
    for (unsigned i = 0; i < 3; ++i)
    {
        if (i == 0)
            hLateSumChan[i]->Draw("hist");
        else
            hLateSumChan[i]->Draw("histsame");
    }
    can0->BuildLegend();

    TCanvas *can1 = new TCanvas("canLateSumLevel1", "canLateSumLevel1");
    for (unsigned i = 5; i >= 3; --i)
    {
        if (i == 5)
            hLateSumChan[i]->Draw("hist");
        else
            hLateSumChan[i]->Draw("histsame");
    }
    can1->BuildLegend();

    TCanvas *can2 = new TCanvas("canLateSumLevel2", "canLateSumLevel2");
    for (unsigned i = 6; i < 9; ++i)
    {
        if (i == 6)
            hLateSumChan[i]->Draw("hist");
        else
            hLateSumChan[i]->Draw("histsame");
    }
    can2->BuildLegend();

    // Open file and get Tree
    TTree *ntuple = (TTree *)fin->Get("ntLateSum");
    Long64_t nentries = ntuple->GetEntries();
    printf("Tree name: %s\n", ntuple->GetName());
    ntuple->Print();

    // Create variables to hold data
    float event;
    ntuple->SetBranchAddress("event", &event);
    float chan;
    ntuple->SetBranchAddress("chan", &chan);
    float geo;
    ntuple->SetBranchAddress("geo", &geo);
    float lateSum;
    ntuple->SetBranchAddress("lateSum", &lateSum);

    // Loop
    double lateSumAve[NONSUMCHANNELS] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    printf("number of entries in ntLateSum %lld \n", nentries);

    for (Long64_t i = 0; i < nentries; i++)
    {
        ntuple->GetEntry(i);
        if (i / 10000 * 10000 == i)
            std::cout << "entry " << i << " event " << event << " chan " << chan << " geo " << geo << "  " << lateSum << std::endl;

        if (chan < 9)
            lateSumAve[int(chan)] += lateSum;
    }

    for (int i = 0; i < NONSUMCHANNELS; i++)
    {
        double eff = effGeoFunc(i);
        printf("chan %i late sum average %f eff %f \n", i, lateSumAve[i], eff);
        hIntChannel->SetBinContent(i + 1, lateSumAve[i]);
        hIntChannel->SetBinError(i + 1, sqrt(lateSumAve[i]));
        hGeoIntChannel->SetBinContent(i + 1, lateSumAve[i] / eff);
        hGeoIntChannel->SetBinError(i + 1, sqrt(lateSumAve[i]) / eff);
    }
    printf("write file \n");
    fout->ls();
    fout->Write();
    printf("postMacroLed completed \n");
}