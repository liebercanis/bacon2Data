#include "TFile.h"
#include "TString.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TH1D.h"
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>

TFile *fin;
TFile *fout;

TString getTag(TString fileName)
{
    std::string sfileName = fileName.Data();

    // Find position of first dash
    size_t firstDash = sfileName.find("-"); // position 4

    // Find position of second dash
    size_t secondDash = sfileName.find("-", firstDash + 1); // position 15

    // Extract substring between the two dashes
    std::string extracted = sfileName.substr(firstDash + 1, secondDash - firstDash - 1);
    // Result: "04_16_2026"
    TString tag(extracted.c_str());
    return tag;
}

TCanvas *drawEarlyRatesCanvas(int minChan, int maxChan, std::vector<TGraph *> &graphs, TString canvasName = "EarlyRatesCanvas")
{
    if (graphs.empty())
    {
        printf("ERROR: graphs vector is empty. No graphs to draw.\n");
        return 0;
    }

    int numGraphs = graphs.size();
    printf("Drawing minChan=%d, maxChan=%d early rate graphs on canvas %s using TMultigraph\n", minChan, maxChan, canvasName.Data());

    // Create TMultiGraph
    TMultiGraph *mg = new TMultiGraph();
    mg->SetName(canvasName);
    mg->SetTitle("Early Hit Rates by Channel");

    // Add all graphs to multigraph with different colors
    int colorArray[] = {kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kYellow, kOrange, kViolet, kSpring, kTeal, kAzure};

    for (int i = minChan; i <= maxChan; ++i)
    {
        printf("Adding graph for channel %d %s to multigraph\n", i, graphs[i]->GetName());
        if (graphs[i])
        {
            graphs[i]->SetLineColor(colorArray[i % 12]);
            graphs[i]->SetMarkerColor(colorArray[i % 12]);
            graphs[i]->SetLineWidth(2);
            graphs[i]->SetMarkerStyle(20);
            mg->Add(graphs[i], "P");
        }
    }

    // Create canvas and draw multigraph
    TCanvas *canvas = new TCanvas(canvasName, canvasName, 900, 600);
    mg->Draw("APL");
    mg->SetTitle("Early Hit Rates by Channel");
    mg->GetXaxis()->SetTitle("File Number");
    mg->GetYaxis()->SetTitle("Early Hit Rate (count/sec)");

    // Create legend with channel labels
    TLegend *leg = new TLegend(0.7, 0.7, 0.99, 0.99);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    for (int i = minChan; i <= maxChan; ++i)
    {
        if (graphs[i])
        {
            leg->AddEntry(graphs[i], Form("Channel %d", i), "lp");
        }
    }
    leg->Draw();

    gPad->SetGrid();
    canvas->Update();

    return canvas;
}

void noiseMacro(TString fileName = "")
{

    /* set bad channels */
    if (fileName.Sizeof() == 0)
    {
        printf("ERROR: fileName argument is empty. Usage: noiseMacro(\"noiseMacro-XX_XX_XXXX-XX_XX_XXXX-XXXXXXX.root\")\n");
        return;
    }
    // open file
    fin = new TFile(fileName, "readonly");
    TString tag = getTag(fileName);
    if (fin->IsZombie())
    {
        printf("no file %s \n", fileName.Data());
        return;
    }

    TNtuple *ntHitCount = (TNtuple *)fin->Get("ntHitCount");
    if (!ntHitCount)
    {
        printf("ERROR: ntuple 'ntHitCount' not found in file %s with date tag %s\n", fileName.Data(), tag.Data());
        return;
    }
    printf("Retrieved ntuple 'ntHitCount' from file %s with date tag %s entries: %lld\n", fileName.Data(), tag.Data(), ntHitCount->GetEntries());

    // open output file
    fout = new TFile(Form("noiseMacro-%s.root", tag.Data()), "recreate");
    // add some histograms to output file
    vector<TH1D *> hNoiseEarly(12);
    for (unsigned i = 0; i < 9; ++i)
    {
        hNoiseEarly[i] = new TH1D(Form("hNoiseEarlyChan%u", i), Form("Early Hit Count rate Channel%u", i), 100, 1.E-7, 1.E-5);
    }
    for (unsigned i = 9; i < 12; ++i)
    {
        hNoiseEarly[i] = new TH1D(Form("hNoiseEarlyChan%u", i), Form("Early Hit Count rate Channel%u", i), 100, 1.E-4, 1.E-3);
    }

    // Declare variables to hold branch data
    Float_t file, nev, chan, earlyCount, lateCount;

    // Set branch addresses
    ntHitCount->SetBranchAddress("file", &file);
    ntHitCount->SetBranchAddress("nev", &nev);
    ntHitCount->SetBranchAddress("chan", &chan);
    ntHitCount->SetBranchAddress("early", &earlyCount);
    ntHitCount->SetBranchAddress("late", &lateCount);

    // Loop over entries and fill histograms
    Long64_t nEntries = ntHitCount->GetEntries();
    std::vector<float> files;
    std::vector<float> earlyRates[12];
    printf("Looping over %lld entries in ntHitCount\n", nEntries);
    int oldFile = -1;
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        ntHitCount->GetEntry(i);
        hNoiseEarly[chan]->Fill(earlyCount / 600. / nev);

        // new file
        if (file != oldFile)
        {
            oldFile = file;
            files.push_back(file);
        }
        float val = earlyCount / 600. / nev;
        // fill for each channel
        for (int ichan = 0; ichan < 12; ++ichan)
            if (chan == ichan)
                earlyRates[ichan].push_back(val);
    }

    // check early rates
    if (0)
    {
        for (int ichan = 0; ichan < 12; ++ichan)
        {
            printf("\t\tichan = %d\n", ichan);
            for (int k = 0; k < earlyRates[ichan].size(); ++k)
                printf("file %.0f  chan %i  entry %i rate %E ....... ", files[k], k, ichan, earlyRates[ichan][k]);
            printf("\n");
        }
    }
    // Create TGraphs for each channel
    std::vector<TGraph *> graphEarlyRates(12);
    for (int chan = 0; chan < 12; ++chan)
    {
        for (size_t j = 0; j < earlyRates[chan].size(); ++j)
        {
            if (earlyRates[chan][j] > 1)
            {
                printf("Channel %d, Point %lu: File %.0f, Rate %E\n", chan, j, files[j], earlyRates[chan][j]);
            }
        }
        printf("Creating TGraph for channel %d fiiles %lu with %lu points\n", chan, files.size(), earlyRates[chan].size());
        graphEarlyRates[chan] = new TGraph(files.size(), &files[0], &earlyRates[chan][0]);
        graphEarlyRates[chan]->SetName(Form("gEarlyRatesChan%d", chan));
        graphEarlyRates[chan]->SetTitle(Form("Early Hit Rate vs File - Channel %d", chan));
    }

    // Draw multigraph canvas
    TCanvas *canvas = drawEarlyRatesCanvas(0, 8, graphEarlyRates, Form("EarlyRatesNonTrig-%s", tag.Data()));
    TCanvas *canvas2 = drawEarlyRatesCanvas(9, 11, graphEarlyRates, Form("EarlyRatesTrig-%s", tag.Data()));

    // Write histograms and graphs to output file
    for (unsigned i = 0; i < 12; ++i)
    {
        hNoiseEarly[i]->Write();
        graphEarlyRates[i]->Write();
    }

    canvas->Write();
    canvas2->Write();

    vector<double> vmean;
    vector<double> vrms;
    vector<double> vchan;
    vector<double> vchanError;

    printf("summary for %s noise rates per event\n", tag.Data());
    for (int ichan = 0; ichan < 12; ++ichan)
    {
        double mean = hNoiseEarly[ichan]->GetMean();
        double rms = hNoiseEarly[ichan]->GetRMS();
        printf("chan %i mean early rate %.2E rms %.2E \n", ichan, mean, rms);
        vmean.push_back(mean);
        vrms.push_back(rms);
        vchan.push_back(ichan);
        vchanError.push_back(0);
    }
    printf("summary for %s noise rates per event %lu %lu\n", tag.Data(), vmean.size(), vrms.size());
    TGraph *gMean = new TGraphErrors(vchan.size(), &vchan[0], &vmean[0], &vchanError[0], &vrms[0]);
    gMean->SetName(Form("gNoiseMeanRms-%s", tag.Data()));
    gMean->SetTitle(Form("Noise Mean vs RMS - for %s", tag.Data()));
    gMean->SetMarkerStyle(21);
    gMean->SetMarkerColor(kRed);
    gMean->SetLineColor(kRed);
    gMean->GetXaxis()->SetTitle("channel number");
    gMean->GetYaxis()->SetTitle("mean noise rate per event ");
    TCanvas *canMeanRMS = new TCanvas(Form("canNoiseMeanRMS-%s", tag.Data()), Form("Noise Mean vs RMS - %s", tag.Data()));
    gMean->Draw("AP");
    gPad->SetGrid();
    gPad->SetLogy();
    // canMeanRMS->Update();
    canMeanRMS->Write();
    fout->Write();
}