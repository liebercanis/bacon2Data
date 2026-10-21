
/*
    pickup graphs from postMacroFile.C
*/
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "distanceLevels.hh"
#include "modelAllFit.hh"
TFile *fout;
TString fileName;
TGraph *gr;
vector<TGraph *> gpeakGraph;
vector<TGraph *> gsinglet;
vector<TGraph *> glate;
vector<TGraph *> gtotal;
vector<double> ppmFile;
vector<int> ppmFileNumber;
std::vector<double> singletAve;
std::vector<double> singletRelative;
vector<vector<double>> absorb;
vector<vector<double>> modelSingletVal;
vector<vector<double>> modelLateVal;
vector<vector<double>> modelTotalVal;
vector<vector<double>> peakVal;
vector<vector<double>> singletVal;
vector<vector<double>> lateVal;
vector<vector<double>> totalVal;
vector<vector<TH1D *>> hEffNormByFile;   // by [file][channel]
vector<vector<double>> singletValByFile; // by [file][channel]

void fitBreiWigner(int ichan, TH1D *h, double &normVal, double &normErr, double &meanVal, double &meanErr)
{
    double peakX = h->GetBinCenter(h->GetMaximumBin());
    double peakY = h->GetMaximum();
    double gammaGuess = 20.;
    double normGuess = peakY * TMath::Pi() * gammaGuess / 2.;

    TF1 *fbw = new TF1(Form("fbwChan%i", ichan), "breitwigner", peakX - 100, peakX + 100);
    fbw->SetParNames("norm", "mean", "gamma");
    fbw->SetParameters(normGuess, peakX, gammaGuess);
    h->Fit(fbw, "RQ");
    fbw->SetLineColor(kRed);
    normVal = fbw->GetParameter(0);
    normErr = fbw->GetParError(0);
    meanVal = fbw->GetParameter(1);
    meanErr = fbw->GetParError(1);
}

void makeIntegralGraphs()
{

    int skip = 0;
    ppmFile.push_back(0.001);
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
    ppmFile.push_back(15.);
    ppmFile.push_back(30.);

    /* need ppm file number */
    for (unsigned ifile = 0; ifile < ppmFile.size(); ++ifile)
        ppmFileNumber.push_back(skip + ifile);

    hEffNormByFile.resize(ppmFile.size());
    singletValByFile.resize(ppmFile.size());
    for (int i = 0; i < ppmFile.size(); ++i)
    {
        hEffNormByFile[i].resize(NCHAN);   // number of channels
        singletValByFile[i].resize(NCHAN); // number of channels
    }

    // Initialize model parameters and optical properties from modelAllFit.hh
    setupModelAllFit();

    // ppm of file
    absorb.resize(ppmFile.size());
    for (int ifile = 0; ifile < ppmFile.size(); ++ifile)
    {
        {
            for (int ilevel = 0; ilevel < NLEVELS; ++ilevel)
            {
                absorb[ifile].push_back(Absorption(ppmFile[ifile], distanceLevel[ilevel]));
                printf("line48  file %i level %i ab %f \n", ifile, ilevel, absorb[ifile].back());
            }
        }
    }

    /*
    for (unsigned ifile = 0; ifile < absorb.size(); ++ifile)
    {
        for (unsigned il = 0; il < absorb[ifile].size(); ++il)
            printf("file %u level %u absorb %E\n", ifile, il, absorb[ifile][il]);
    }
    */

    // for (unsigned ifile = 0; ifile < ppmFile.size(); ++ifile)
    //   tbDraw(-1, ppmFile[ifile]);

    fout = new TFile("integralGraphs.root", "recreate");
    gsinglet.resize(ppmFile.size(), nullptr);
    gtotal.resize(ppmFile.size(), nullptr);
    glate.resize(ppmFile.size(), nullptr);
    gpeakGraph.resize(ppmFile.size(), nullptr);
    for (int ifile = 0; ifile < ppmFileNumber.size(); ++ifile)
    {
        fileName.Form("postMacroAllFile%i.root", ppmFileNumber[ifile]);
        cout << " open file  " << fileName.Data() << endl;
        TFile *fin = new TFile(fileName, "readonly");
        TString graphName;
        graphName.Form("gSingletIntegralsFile%i", ppmFileNumber[ifile]);
        // cout << ifile << " find graph " << graphName.Data() << endl;
        TGraph *gr = nullptr;
        fin->GetObject(graphName, gr);
        if (gr != nullptr)
        {
            gsinglet[ifile] = gr;
            fout->Add(gr);
            cout << ifile << " got graph " << gr->GetName() << endl;
        }

        graphName.Form("gPeakCorr");
        // cout << ifile << " find graph " << graphName.Data() << endl;
        gr = nullptr;

        fin->GetObject(graphName, gr);
        if (gr != nullptr)
        {
            gpeakGraph[ifile] = gr;
            fout->Add(gr);
            cout << ifile << " got graph " << gr->GetName() << endl;
        }

        graphName.Form("gLateIntegralsFile%i", ppmFileNumber[ifile]);
        // cout << ifile << " find graph " << graphName.Data() << endl;
        gr = nullptr;
        fin->GetObject(graphName, gr);
        if (gr != nullptr)
        {
            cout << " got graph " << gr->GetName() << endl;
            fout->Add(gr);
            glate[ifile] = gr;
        }

        graphName.Form("gTotalIntegralsFile%i", ppmFileNumber[ifile]);
        // cout << ifile << " find graph " << graphName.Data() << endl;
        gr = nullptr;
        fin->GetObject(graphName, gr);
        if (gr != nullptr)
        {
            cout << " got graph " << gr->GetName() << endl;
            fout->Add(gr);
            gtotal[ifile] = gr;
        }
        fin->Close();
    }
    printf("line 81 ppm file %li gsinglet %li gtotal %li \n", ppmFile.size(), gsinglet.size(), gtotal.size());
    printf("glate size %lu \n", glate.size());
    // for (unsigned ilate = 0; ilate < glate.size(); ++ilate)
    //     printf("%i %s \n", ilate, glate[ilate]->GetName());
    //  fout->ls();
    modelSingletVal.resize(gsinglet[0]->GetN()); // number of channels
    modelLateVal.resize(glate[0]->GetN());       // number of channels
    modelTotalVal.resize(glate[0]->GetN());      // number of channels
    for (int ifile = 0; ifile < ppmFileNumber.size(); ++ifile)
    {
        fileName.Form("tbModelPPM%i.root", ppmFileNumber[ifile]);
        printf("line 110 open %s \n", fileName.Data());
        TFile *fin = new TFile(fileName, "readonly");
        for (int ichan = 0; ichan < NCHAN; ++ichan)
        {
            TString histName;
            TH1D *hist = nullptr;
            histName.Form("fitWaveFitChan%i", ichan);
            fin->GetObject(histName, hist);
            if (!hist)
            {
                printf("cannot get %s from %s \n", histName.Data(), fileName.Data());
                modelSingletVal[ichan].push_back(0);
                modelLateVal[ichan].push_back(0);
                modelTotalVal[ichan].push_back(0);
                continue;
            }
            // TH1D *hcopy = (TH1D *)hist->Clone(Form("fitWaveFitChan%iFile%i", ichan, ifile));
            // hcopy->SetTitle(Form("fitWaveChan%iFile%i %.3fPPM", ichan, ifile, ppmFile[ifile]));
            //  integral is over bins
            // printf("get %s \n", hist->GetName());
            double singletIntegral = hist->Integral(earlyBin, lateBin);
            double lateIntegral = hist->Integral(lateBin, 3000 / 2);
            double totalIntegral = hist->Integral(earlyBin, 15000 / 2);
            modelSingletVal[ichan].push_back(singletIntegral);
            modelLateVal[ichan].push_back(lateIntegral);
            modelTotalVal[ichan].push_back(totalIntegral);

            // printf("line128 !!!!!! chan %i file %i PPM %.3f singlet integral %.3E late integral %.3E total integral%.3E !!!!!!\n", ichan, ifile, ppmFile[ifile], singletIntegral, lateIntegral, totalIntegral);
            //  int jfile = modelSingletVal[ichan].size() - 1;
            if (ichan == 7)
                printf("\t line141  ******* modelSingletVal ich %i ifile %i val %f \n", ichan, ifile, modelSingletVal[ichan][ifile]);
            //   hist->SetLineColor(ifile + 1);
            //   fout->Add(hist);
        }
        fin->Close();
    }
    /* pick up normalized Light Curves*/
    for (int ifile = 0; ifile < ppmFile.size(); ++ifile)
    {
        fileName.Form("postMacroAllFile%i.root", ifile);
        TFile *fin = new TFile(fileName, "readonly");
        for (int ichan = 0; ichan < NCHAN; ++ichan)
        {
            TString histName;
            TH1D *hist = nullptr;
            histName.Form("effNormChan%iFile%i", ichan, ifile);
            fin->GetObject(histName, hist);
            if (hist != nullptr)
            {
                hist->SetDirectory(0); // detach so it survives fin->Close()
            }
            else
            {
                printf("line 166 file %i chan %i hist %s not found \n", ifile, ichan, histName.Data());
            }
            hEffNormByFile[ifile][ichan] = hist;
        }
        fin->Close();
    }
    for (int ifile = 0; ifile < ppmFile.size(); ++ifile)
    {
        printf("line 172 file %i \n", ifile);
        for (int ichan = 0; ichan < NCHAN; ++ichan)
        {
            printf("line 175 file %i chan %i hist %s integral %.3f \n", ifile, ichan, hEffNormByFile[ifile][ichan]->GetName(), hEffNormByFile[ifile][ichan]->Integral());
            fout->Add(hEffNormByFile[ifile][ichan]);
        }
    }

    // printf("modelSingletVal size %lu \n", modelSingletVal.size());
    for (int ichan = 0; ichan < modelSingletVal.size(); ++ichan)
    {
        printf("line 148 modelSingletVal ich %i size %lu \n", ichan, modelSingletVal[ichan].size());
        if (ichan >= 12)
            continue;
        for (int ifile = 0; ifile < modelSingletVal[ichan].size(); ++ifile)
        {
            printf("line 129 modelSingletVal ich %i ifile %i val %f \n", ichan, ifile, modelSingletVal[ichan][ifile]);
        }
    }

    // each one of these corresponds to a file
    // for (unsigned i = 0; i < gsinglet.size(); ++i)
    //    printf("%f %s %s \n", ppmFile[i], gsinglet[i]->GetName(), glate[i]->GetName())

    // get points from graphs
    singletVal.resize(gsinglet[0]->GetN()); // number of channels
    // by channel get points from each file singletVal[ich][ifile]
    for (unsigned ifile = 0; ifile < gsinglet.size(); ++ifile) // loop over all files
    {
        if (gsinglet[ifile] == nullptr)
            continue;
        Double_t x, y;
        for (int ic = 0; ic < gsinglet[ifile]->GetN(); ++ic) // loop over channels
        {
            if (ic == 0)
                continue;
            gsinglet[ifile]->GetPoint(ic, x, y);
            singletVal[ic].push_back(y);
            singletValByFile[ifile][ic] = y;
            printf("line150 SINGLETVAL file %i ic %i PPM %f point %lu singlet %f\n", ifile, ic, ppmFile[ifile], singletVal[ic].size() - 1, singletVal[ic].back());
        }
    }

    singletAve.resize(4, 0);
    std::vector<int> nlevel;
    nlevel.resize(4, 0);
    for (int ich = 0; ich < singletValByFile[0].size(); ++ich)
    {
        // get level
        if (ich == 0 || ich == 8 || ich > 11)
            continue;
        int ilevel = getLevel(ich);
        if (ilevel > 3)
            continue;
        nlevel[ilevel] += 1;
        if (ilevel == 0)
            singletAve[ilevel] += singletValByFile[0][ich];
        else if (ilevel == 1)
            singletAve[ilevel] += singletValByFile[0][ich];
        else if (ilevel == 2)
            singletAve[ilevel] += singletValByFile[0][ich];
        else if (ilevel == 3)
            singletAve[ilevel] += singletValByFile[0][ich];

        printf("ich %i ilevel %i nlevel %i singlet %.3E \n", ich, ilevel, nlevel[ilevel], singletValByFile[0][ich]);
    }
    for (unsigned ilevel = 0; ilevel < singletAve.size(); ++ilevel)
    {
        singletAve[ilevel] /= double(nlevel[ilevel]);
    }

    singletRelative.resize(NCHAN, 1.);
    for (unsigned ich = 0; ich < NCHAN; ++ich)
    {
        int ilevel = getLevel(ich);
        if (ich == 0 || ich == 8 || ich > 11)
            continue;
        singletRelative[ich] = 1. + (singletValByFile[0][ich] - singletAve[ilevel]) / singletAve[ilevel];
        printf("ich %i ilevel %i nlevel %i singletAve %.3E  val %.3E rel %.3E diff %.3E\n", ich, ilevel, nlevel[ilevel], singletAve[ilevel], singletValByFile[0][ich], singletRelative[ich], 1. + (singletValByFile[0][ich] - singletAve[ilevel]) / singletAve[ilevel]);
    }

    for (unsigned ich = 0; ich < NCHAN; ++ich)
    {
        int ilevel = getLevel(ich);
        if (ich == 0 || ich == 8 || ich > 11)
            continue;
        printf("singletNorm[%i] = %.3E ; \n", ich, singletRelative[ich]);
    }

    for (int ic = 0; ic < singletVal.size(); ++ic)
    { // loop over channels{
        for (int ifile = 0; ifile < singletVal[ic].size(); ++ifile)
        {
            printf("line93 singletVal ich %i ifile %i val %f \n", ic, ifile, singletVal[ic][ifile]);
        }
    }

    printf("line95 singlet val size %lu\n", singletVal.size());

    lateVal.resize(glate[0]->GetN());           // number of channels
    for (unsigned i = 0; i < glate.size(); ++i) // loop over all files
    {
        if (glate[i] == nullptr)
            continue;
        Double_t x, y;
        for (int ich = 0; ich < glate[i]->GetN(); ++ich) // loop over channels
        {
            if (ich == 0)
                continue;
            glate[i]->GetPoint(ich, x, y);
            lateVal[ich].push_back(y);
            // printf("file %i ich %i PPM %f late %f \n", i, ich, ppmFile[i], lateVal[ich].back());
        }
    }
    printf("line111 late  val size %lu\n", lateVal.size());

    totalVal.resize(gtotal[0]->GetN());          // number of channels
    for (unsigned i = 0; i < gtotal.size(); ++i) // loop over all files
    {
        if (gtotal[i] == nullptr)
            continue;
        Double_t x, y;
        for (int ich = 0; ich < gtotal[i]->GetN(); ++ich) // loop over channels
        {
            if (ich == 0)
                continue;
            gtotal[i]->GetPoint(ich, x, y);
            totalVal[ich].push_back(y);
            printf("line204 file %i ich %i PPM %f total %f \n", i, ich, ppmFile[i], totalVal[ich].back());
        }
    }
    printf("line111 total  val size %lu\n", totalVal.size());

    peakVal.resize(gpeakGraph.size());                           // number of files
    for (unsigned ifile = 0; ifile < gpeakGraph.size(); ++ifile) // loop over all files
    {
        Double_t x, y;
        if (gpeakGraph[ifile] == nullptr)
            continue;
        for (int ic = 0; ic < gpeakGraph[ifile]->GetN(); ++ic) // loop over channels
        {
            if (ic == 0)
                continue;
            gpeakGraph[ifile]->GetPoint(ic, x, y);
            peakVal[ifile].push_back(y);
            printf("line311 file %i ic %i PPM %f peak %f\n", ifile, ic, ppmFile[ifile], peakVal[ifile].back());
        }
    }
    printf("line311 peak val size %lu\n", peakVal.size());

    /* setup multigraph peak*/
    TMultiGraph *mgPeak[4]; // one for each level
    for (int ilevel = 0; ilevel < 4; ++ilevel)
    {
        mgPeak[ilevel] = new TMultiGraph(Form("mgPeakLevel%i", ilevel), Form("Peak values level %i vs PPM", ilevel));
        mgPeak[ilevel]->SetName(Form("mgPeakLevel%i", ilevel));
        printf("level %i %s \n", ilevel, mgPeak[ilevel]->GetName());
    }

    /* peak value vs PPM, one canvas per channel */
    int nPeakChan = peakVal.size() > 0 ? peakVal[0].size() : 0;
    for (int ich = 0; ich < nPeakChan; ++ich)
    {
        vector<double> ppmForChan;
        vector<double> peakForChan;
        for (unsigned ifile = 0; ifile < peakVal.size(); ++ifile)
        {
            if (ich >= (int)peakVal[ifile].size())
                continue;
            ppmForChan.push_back(ppmFile[ifile]);
            peakForChan.push_back(peakVal[ifile][ich]);
        }
        if (ppmForChan.size() == 0)
            continue;
        int ichan = ich + 1; // channel 0 is skipped when filling peakVal
        int ilevel = getLevel(ichan);
        TGraph *gPeakPPM = new TGraph(ppmForChan.size(), &ppmForChan[0], &peakForChan[0]);
        gPeakPPM->SetName(Form("gPeakPPMChan%i", ichan));
        gPeakPPM->SetTitle(Form("gPeakPPMChan%i", ichan));
        gPeakPPM->SetMarkerStyle(21);
        gPeakPPM->SetMarkerColor(ichan + 1);
        gPeakPPM->SetLineColor(ichan + 1);
        fout->Add(gPeakPPM);
        if (ilevel <= 3 && ichan != 0 && ichan != 8)
            mgPeak[ilevel]->Add(gPeakPPM, "p");

        TCanvas *canPeakPPM = new TCanvas(Form("canPeakPPMChan%i", ichan), Form("canPeakPPMChan%i", ichan));
        canPeakPPM->SetGrid();
        canPeakPPM->SetLogx();
        gPeakPPM->GetHistogram()->GetXaxis()->SetTitle("PPM");
        gPeakPPM->GetHistogram()->GetYaxis()->SetTitle("peak value");
        gPeakPPM->Draw("ap");
        canPeakPPM->Print(Form("peakPPMChan%i.pdf", ichan));
        fout->Append(canPeakPPM);
    }

    /* peak level plots: channels for each level superimposed */
    for (int ilevel = 0; ilevel < 4; ++ilevel)
    {
        TCanvas *cmgPeak = new TCanvas(Form("canPeakMultiLevel%i", ilevel), Form("canPeakMultiLevel%i", ilevel));
        cmgPeak->SetGrid();
        cmgPeak->SetLogx();
        mgPeak[ilevel]->GetXaxis()->SetTitle("PPM");
        mgPeak[ilevel]->GetYaxis()->SetTitle(Form("peak value level %i", ilevel));
        mgPeak[ilevel]->Draw("a");
        cmgPeak->BuildLegend();
        cmgPeak->Print(Form("peakMultiLevel%i.pdf", ilevel));
        fout->Add(mgPeak[ilevel]);
        fout->Append(cmgPeak);
    }

    // get singlet integrals from model make graphs for each level
    /* setup multigraph singlet*/
    TGraph *singletPPM;
    TMultiGraph *mgSinglet[4]; // one for each level
    for (int ilevel = 0; ilevel < 4; ++ilevel)
    {
        mgSinglet[ilevel] = new TMultiGraph(Form("mgSingletLevel%i", ilevel), Form("Singlet Integrals level %i vs PPM", ilevel));
        mgSinglet[ilevel]->SetName(Form("mgSingletLevel%i", ilevel));
        printf("level %i %s \n", ilevel, mgSinglet[ilevel]->GetName());
    }
    /* setup multigraph late*/
    TGraph *latePPM;
    TMultiGraph *mgLate[4]; // one for each level
    for (int ilevel = 0; ilevel < 4; ++ilevel)
    {
        mgLate[ilevel] = new TMultiGraph(Form("mgLateLevel%i", ilevel), Form("Late Integrals level %i vs PPM", ilevel));
        mgLate[ilevel]->SetName(Form("mgLateLevel%i", ilevel));
        printf("level %i %s \n", ilevel, mgLate[ilevel]->GetName());
    }

    TGraph *totalPPM;
    TMultiGraph *mgTotal[4]; // one for each level
    for (int ilevel = 0; ilevel < 4; ++ilevel)
    {
        mgTotal[ilevel] = new TMultiGraph(Form("mgTotalLevel%i", ilevel), Form("Total Integrals level %i vs PPM", ilevel));
        mgTotal[ilevel]->SetName(Form("mgTotalLevel%i", ilevel));
        printf("level %i %s \n", ilevel, mgTotal[ilevel]->GetName());
    }

    /*
        make model graphs
    */

    /* scale singlet model graph */
    for (int ich = 0; ich < modelSingletVal.size(); ++ich)
    {
        if (ich == 0)
            continue;
        // printf("line 149 modelSingletVal ich %i size %lu \n", ich, modelSingletVal[ich].size());
        if (modelSingletVal[ich].size() == 0)
            continue;
        /* normalize */
        double ratio = singletVal[ich][0] / modelSingletVal[ich][0];
        for (int ifile = 0; ifile < modelSingletVal[ich].size(); ++ifile)
        { // loop over files}
            modelSingletVal[ich][ifile] = modelSingletVal[ich][ifile] * ratio;
            printf("line 215 modelSingletVal ifile %i ppm %.3f modelSingletVal %f  \n", ifile, ppmFile[ifile], modelSingletVal[ich][ifile]);
        }
    }

    /* scale late integral has to be compared do model no (1-A) */
    for (int ich = 0; ich < modelLateVal.size(); ++ich)
    {
        if (ich == 0)
            continue;
        printf("line 149 modelLateVal ich %i size %lu \n", ich, modelLateVal[ich].size());
        if (modelLateVal[ich].size() == 0)
            continue;
        /* normalize */
        double ratio = lateVal[ich][0] / modelLateVal[ich][0];
        for (int ifile = 0; ifile < modelLateVal[ich].size(); ++ifile)
        { // loop over files}
            modelLateVal[ich][ifile] = modelLateVal[ich][ifile] * ratio;
            // printf("line 151 modelLateVal ifile %i ppm %.3f modelLateVal %f  \n", ifile, ppmFile[ifile], modelLateVal[ich][ifile]);
        }
    }

    /* scale total integral has to be compared do model no (1-A) */
    for (int ich = 0; ich < modelTotalVal.size(); ++ich)
    {
        if (ich == 0)
            continue;
        printf("line 149 modelTotalVal ich %i size %lu \n", ich, modelTotalVal[ich].size());
        if (modelTotalVal[ich].size() == 0)
            continue;
        /* normalize */
        double ratio = totalVal[ich][0] / modelTotalVal[ich][0];
        for (int ifile = 0; ifile < modelTotalVal[ich].size(); ++ifile)
        { // loop over files}
            modelTotalVal[ich][ifile] = modelTotalVal[ich][ifile] * ratio;
            // printf("line 151 modelTotalVal ifile %i ppm %.3f modelTotalVal %f  \n", ifile, ppmFile[ifile], modelTotal    Val[ich][ifile]);
        }
    }

    /* singlet model graph */
    TGraph *modelSinglet = nullptr;
    int oldLevel = -1;
    for (int ich = 0; ich < gsinglet[0]->GetN(); ++ich)
    {
        int ilevel = getLevel(ich);
        if (ilevel > 3)
            continue;
        if (ilevel != oldLevel)
        {
            oldLevel = ilevel;
        }
        // else
        //     continue;
        printf("line 134 %i %lu %lu \n", ich, singletVal[ich].size(), modelSingletVal[ich].size());
        modelSinglet = new TGraph(modelSingletVal[ich].size(), &ppmFile[0], &modelSingletVal[ich][0]);
        modelSinglet->SetName(Form("modelSingletLevel%ichan%i", ilevel, ich));
        modelSinglet->SetTitle(Form("modelSingletLevel%ichan%i", ilevel, ich));
        modelSinglet->SetLineColor(kBlue);
        modelSinglet->SetLineWidth(2);
        if (ich == 0 || ich == 8)
            continue;
        mgSinglet[ilevel]->Add(modelSinglet, "line");
        printf("level %i %s \n", ilevel, modelSinglet->GetName());
        fout->Add(modelSinglet);
    }

    /* scale late model graph */
    oldLevel = -1;
    TGraph *modelLate = nullptr;
    for (int ich = 0; ich < gsinglet[0]->GetN(); ++ich)
    {
        int ilevel = getLevel(ich);
        if (ilevel > 3)
            continue;
        if (ilevel != oldLevel)
        {
            oldLevel = ilevel;
        }
        // else
        //     continue;
        //  printf("line 134 %i %lu %lu \n", ich, singletVal[ich].size(), modelLateVal[ich].size());
        modelLate = new TGraph(modelLateVal[ich].size(), &ppmFile[0], &modelLateVal[ich][0]);
        modelLate->SetName(Form("modelLateLevel%ichan%i", ilevel, ich));
        modelLate->SetTitle(Form("modelLateLevel%ichan%i", ilevel, ich));
        modelLate->SetLineColor(kBlue);
        modelLate->SetLineWidth(2);
        if (ich == 0 || ich == 8)
            continue;
        mgLate[ilevel]->Add(modelLate, "line");
        printf("level %i %s \n", ilevel, modelLate->GetName());
        fout->Add(modelLate);
    }

    oldLevel = -1;
    TGraph *modelTotal = nullptr;
    for (int ich = 0; ich < gsinglet[0]->GetN(); ++ich)
    {
        int ilevel = getLevel(ich);
        if (ilevel > 3)
            continue;
        if (ilevel != oldLevel)
        {
            oldLevel = ilevel;
        }
        // else
        //     continue;
        //  printf("line 134 %i %lu %lu \n", ich, singletVal[ich].size(), modelLateVal[ich].size());
        modelTotal = new TGraph(modelTotalVal[ich].size(), &ppmFile[0], &modelTotalVal[ich][0]);
        modelTotal->SetName(Form("modelTotalLevel%ichan%i", ilevel, ich));
        modelTotal->SetTitle(Form("modelTotalLevel%ichan%i", ilevel, ich));
        modelTotal->SetLineColor(kBlue);
        modelTotal->SetLineWidth(2);
        if (ich == 0 || ich == 8)
            continue;
        mgTotal[ilevel]->Add(modelTotal, "line");
        printf("level %i %s \n", ilevel, modelTotal->GetName());
        fout->Add(modelTotal);
    }

    for (int ich = 0; ich < singletVal.size(); ++ich)
    {
        if (ich == 0) // channel 0 has no data: skipped when filling singletVal
            continue;
        int ilevel = getLevel(ich);
        if (ilevel > 3)
            continue;
        if (singletVal[ich].size() == 0)
            continue;
        printf("level %i %s \n", ilevel, mgSinglet[ilevel]->GetName());
        singletPPM = new TGraph(singletVal[ich].size(), &ppmFile[0], &singletVal[ich][0]);
        singletPPM->SetName(Form("singletPPMChan%i", ich));
        singletPPM->SetTitle(Form("singletPPMChan%i", ich));
        singletPPM->SetMarkerStyle(21);
        singletPPM->SetMarkerColor(ich + 1);
        singletPPM->SetLineColor(ich + 1);
        if (ich == 9)
            singletPPM->SetMarkerColor(kGreen);
        fout->Add(singletPPM);
        if (ich != 0 && ich != 8)
            mgSinglet[ilevel]->Add(singletPPM, "p");
        TCanvas *csinglet = new TCanvas(Form("canSingletPPMChan%i", ich), "canSingletPPM");
        csinglet->SetName(Form("canSingletPPMChan%i", ich));
        csinglet->SetTitle(Form("canSingletPPMChan%i", ich));
        singletPPM->SetMarkerStyle(21);
        singletPPM->GetHistogram()->GetXaxis()->SetTitle("PPM");
        singletPPM->GetHistogram()->GetYaxis()->SetTitle("singlet integral");
        csinglet->SetGrid();
        csinglet->SetLogx();
        singletPPM->Draw("ap");
        csinglet->Print(".pdf");
        fout->Append(csinglet);
    }
    /* singlet level plots */
    for (int ilevel = 0; ilevel < 4; ++ilevel)
    {
        TCanvas *cmg = new TCanvas(Form("canSingletMultiLevel%i", ilevel), Form("canSingletMultiLevel%i", ilevel));
        cmg->SetGrid();
        cmg->SetLogx();
        mgSinglet[ilevel]->GetXaxis()->SetTitle("PPM");
        mgSinglet[ilevel]->GetYaxis()->SetTitle(Form("singlet integral level %i", ilevel));
        mgSinglet[ilevel]->Draw("a");
        cmg->BuildLegend();
        cmg->Print(Form("singletMultiLevel%i.pdf", ilevel));
        fout->Add(mgSinglet[ilevel]);
        fout->Append(cmg);
    }

    // put everything together
    for (int ich = 0; ich < gsinglet[0]->GetN(); ++ich)
    {
        if (ich == 0)
            continue;
        int ilevel = getLevel(ich);
        if (ilevel > 3)
            continue;
        latePPM = new TGraph(lateVal[ich].size(), &ppmFile[0], &lateVal[ich][0]);
        latePPM->SetName(Form("canLatePPMChan%i", ich));
        latePPM->SetTitle(Form("canLatePPMChan%i", ich));
        latePPM->GetHistogram()->GetXaxis()->SetTitle("PPM");
        latePPM->GetHistogram()->GetYaxis()->SetTitle("late integral");
        latePPM->SetMarkerStyle(21);
        latePPM->SetMarkerStyle(21);
        latePPM->SetMarkerColor(ich + 1);
        latePPM->SetLineColor(ich + 1);
        if (ich == 9)
            latePPM->SetMarkerColor(kGreen);
        // printf("line 213 add ich %i level %i %s\n", ich, ilevel, latePPM->GetName());
        fout->Add(latePPM);

        if (ich != 0 && ich != 8)
        {
            mgLate[ilevel]->Add(latePPM, "p");
            // printf("line 220 add ich %i level %i %s\n", ich, ilevel, latePPM->GetName());
        }
        TCanvas *clate = new TCanvas(Form("canlatePPMChan%i", ich), "canlatePPM");
        clate->SetTitle(Form("canLatePPMChan%i", ich));
        latePPM->SetMarkerStyle(21);
        latePPM->GetHistogram()->GetXaxis()->SetTitle("PPM");
        latePPM->GetHistogram()->GetYaxis()->SetTitle("late integral");
        clate->SetGrid();
        clate->SetLogx();
        latePPM->Draw("ap");
        clate->Print(".pdf");

        fout->Append(clate);
    }
    /*totals */
    for (int ich = 0; ich < gsinglet[0]->GetN(); ++ich)
    {
        if (ich == 0)
            continue;
        int ilevel = getLevel(ich);
        if (ilevel > 3)
            continue;
        totalPPM = new TGraph(totalVal[ich].size(), &ppmFile[0], &totalVal[ich][0]);
        if (ich != 0 && ich != 8)
        {
            mgTotal[ilevel]->Add(totalPPM, "p");
            // printf("line 220 add ich %i level %i %s\n", ich, ilevel, totalPPM->GetName());
        }
        totalPPM->SetName(Form("canTotalPPMChan%i", ich));
        totalPPM->SetTitle(Form("canTotalPPMChan%i", ich));
        totalPPM->GetHistogram()->GetXaxis()->SetTitle("PPM");
        totalPPM->GetHistogram()->GetYaxis()->SetTitle("total integral");
        totalPPM->SetMarkerStyle(21);
        totalPPM->SetMarkerStyle(21);
        totalPPM->SetMarkerColor(ich + 1);
        totalPPM->SetLineColor(ich + 1);
        if (ich == 9)
            totalPPM->SetMarkerColor(kGreen);
        // printf("line 213 add ich %i level %i %s\n", ich, ilevel, totalPPM->GetName());
        TCanvas *cTotal = new TCanvas(Form("canTotalPPMChan%i", ich), "canTotalPPM");
        cTotal->SetTitle(Form("canTotalPPMChan%i", ich));
        totalPPM->SetMarkerStyle(21);
        totalPPM->GetHistogram()->GetXaxis()->SetTitle("PPM");
        totalPPM->GetHistogram()->GetYaxis()->SetTitle("total integral");

        totalPPM->Draw("ap");
        cTotal->SetGrid();
        cTotal->SetLogx();
        cTotal->Print(".pdf");
        fout->Add(totalPPM);

        // printf("line 220 add ich %i level %i %s\n", ich, ilevel, totalPPM->GetName());
    }

    /* late level canvas */
    for (int ilevel = 0; ilevel < 4; ++ilevel)
    {
        TCanvas *cmg = new TCanvas(Form("canLateMultiLevel%i", ilevel), Form("canLateMultiLevel%i", ilevel));
        cmg->SetGrid();
        cmg->SetLogx();
        mgLate[ilevel]->GetXaxis()->SetTitle("PPM");
        mgLate[ilevel]->GetYaxis()->SetTitle(Form("late integral level %i", ilevel));
        mgLate[ilevel]->Draw("a");
        cmg->BuildLegend();
        cmg->Print(Form("lateMultiLevel%i.pdf", ilevel));
        fout->Add(mgLate[ilevel]);
        fout->Append(cmg);
    }

    /* total level canvas */
    for (int ilevel = 0; ilevel < 4; ++ilevel)
    {
        TCanvas *cmg = new TCanvas(Form("canTotalMultiLevel%i", ilevel), Form("canTotalMultiLevel%i", ilevel));
        cmg->SetGrid();
        cmg->SetLogx();
        mgTotal[ilevel]->GetXaxis()->SetTitle("PPM");
        mgTotal[ilevel]->GetYaxis()->SetTitle(Form("total integral level %i", ilevel));
        mgTotal[ilevel]->Draw("a");
        cmg->BuildLegend();
        cmg->Print(Form("totalMultiLevel%i.pdf", ilevel));
        fout->Add(mgTotal[ilevel]);
        fout->Append(cmg);
    }

    /* effNorm integral vs channel number, one curve per file */
    int effNormColors[] = {kRed, kGreen + 2, kBlue, kMagenta, kCyan + 2, kOrange + 7,
                           kSpring, kTeal + 2, kAzure + 2, kViolet, kPink + 9, kBlack};
    TMultiGraph *mgEffNorm = new TMultiGraph("mgEffNormByChannel", "effNorm integral vs channel number");

    for (int ifile = 0; ifile < ppmFile.size(); ++ifile)
    {
        printf("line 569  EffNormByFile for file %i \n", ifile);
        vector<double> chanNum;
        vector<double> chanNumError;
        vector<double> effInt;
        vector<double> effIntError;
        for (int ichan = 0; ichan < hEffNormByFile[ifile].size(); ++ichan)
        {
            if (!hEffNormByFile[ifile][ichan])
                continue;
            chanNum.push_back(ichan);
            chanNumError.push_back(0);
            double error = 0;
            effInt.push_back(hEffNormByFile[ifile][ichan]->IntegralAndError(0, 7500, error));
            effIntError.push_back(error);
            printf("\t \t line 580 ichan %i effNorm integral %f ± %f \n", ichan, effInt.back(), effIntError.back());
        }
        if (chanNum.size() == 0)
            continue;
        TGraphErrors *gEffNorm = new TGraphErrors(chanNum.size(), &chanNum[0], &effInt[0], &chanNumError[0], &effIntError[0]);

        gEffNorm->SetName(Form("gEffNormByChanFile%i", ifile));
        gEffNorm->SetTitle(Form("%.3f PPM", ppmFile[ifile]));
        gEffNorm->SetMarkerStyle(21);
        gEffNorm->SetMarkerColor(effNormColors[ifile % 12]);
        gEffNorm->SetLineColor(effNormColors[ifile % 12]);
        mgEffNorm->Add(gEffNorm, "lp");
        fout->Add(gEffNorm);
    }
    TCanvas *canEffNormByChan = new TCanvas("canEffNormByChan", "effNorm vs channel");
    canEffNormByChan->SetGrid();
    mgEffNorm->GetXaxis()->SetTitle("channel number");
    mgEffNorm->GetYaxis()->SetTitle("effNorm integral");
    mgEffNorm->Draw("a");
    canEffNormByChan->BuildLegend();
    canEffNormByChan->Print("effNormByChannel.pdf");
    fout->Add(mgEffNorm);
    fout->Append(canEffNormByChan);

    /* singlet integral vs channel number, one curve per file */
    TMultiGraph *mgSingletByFile = new TMultiGraph("mgSingletValByFile", "singlet integral vs channel number");
    for (int ifile = 0; ifile < ppmFile.size(); ++ifile)
    {
        vector<double> chanNum;
        vector<double> singletInt;
        for (int ichan = 1; ichan < singletValByFile[ifile].size(); ++ichan) // channel 0 has no data
        {
            chanNum.push_back(ichan);
            singletInt.push_back(singletValByFile[ifile][ichan]);
        }
        if (chanNum.size() == 0)
            continue;
        TGraph *gSingletByFile = new TGraph(chanNum.size(), &chanNum[0], &singletInt[0]);
        gSingletByFile->SetName(Form("gSingletValByChanFile%i", ifile));
        gSingletByFile->SetTitle(Form("%.3f PPM", ppmFile[ifile]));
        gSingletByFile->SetMarkerStyle(21);
        gSingletByFile->SetMarkerColor(effNormColors[ifile % 12]);
        gSingletByFile->SetLineColor(effNormColors[ifile % 12]);
        mgSingletByFile->Add(gSingletByFile, "lp");
        fout->Add(gSingletByFile);
    }
    TCanvas *canSingletValByFile = new TCanvas("canSingletValByFile", "singlet integral vs channel");
    canSingletValByFile->SetGrid();
    mgSingletByFile->GetXaxis()->SetTitle("channel number");
    mgSingletByFile->GetYaxis()->SetTitle("singlet integral");
    mgSingletByFile->Draw("a");
    canSingletValByFile->BuildLegend();
    canSingletValByFile->Print("singletValByChannel.pdf");
    fout->Add(mgSingletByFile);
    fout->Append(canSingletValByFile);

    /* hEffNormByFile waveform for each channel, overlaying every file */
    for (int ichan = 0; ichan < NCHAN; ++ichan)
    {
        bool anyData = false;
        for (int ifile = 0; ifile < ppmFile.size(); ++ifile)
            if (hEffNormByFile[ifile][ichan])
            {
                anyData = true;
                break;
            }
        if (!anyData)
            continue;

        TCanvas *canEffNormWave = new TCanvas(Form("canEffNormWaveChan%i", ichan), Form("effNorm waveform chan %i", ichan));
        canEffNormWave->SetGrid();
        canEffNormWave->SetLogy();
        TLegend *legEffNormWave = new TLegend(0.7, 0.55, 0.9, 0.9);
        bool firstDrawn = false;
        for (int ifile = 0; ifile < ppmFile.size(); ++ifile)
        {
            if (!hEffNormByFile[ifile][ichan])
                continue;
            TH1D *h = hEffNormByFile[ifile][ichan];
            h->SetLineColor(effNormColors[ifile % 12]);
            h->SetTitle(Form("effNorm chan %i", ichan));
            h->GetXaxis()->SetRangeUser(1000, 7500);
            h->GetXaxis()->SetTitle("time [ns]");
            h->GetYaxis()->SetTitle("effNorm yield");
            legEffNormWave->AddEntry(h, Form("%.3f PPM", ppmFile[ifile]), "l");
            h->Rebin(10);
            h->Draw(firstDrawn ? "hist same" : "hist");
            firstDrawn = true;
        }
        legEffNormWave->Draw();
        canEffNormWave->Print(Form("effNormWaveChan%i.pdf", ichan));
        fout->Append(canEffNormWave);
    }

    /* hEffNormByFile waveform for each file overlaying every channel */
    for (int ifile = 0; ifile < ppmFile.size(); ++ifile)
    {
        TCanvas *canEffNormByFile = new TCanvas(Form("canEffNormByFile%i", ifile), Form("effNorm waveform file %i", ifile));
        canEffNormByFile->SetGrid();
        canEffNormByFile->SetLogy();
        TLegend *legEffNormWave = new TLegend(0.7, 0.55, 0.9, 0.9);
        bool firstDrawn = false;
        for (int ichan = 0; ichan < NCHAN; ++ichan)
        {
            if (!hEffNormByFile[ifile][ichan])
                continue;
            if (ichan == 0 || ichan > 8)
                continue; // channel 0 has no data
            TH1D *h = hEffNormByFile[ifile][ichan];
            h->SetLineColor(effNormColors[ichan % 12]);
            h->SetMarkerColor(effNormColors[ichan % 12]);
            h->SetTitle(Form("effNorm chan %i", ichan));
            h->GetXaxis()->SetRangeUser(1000, 15000);
            h->GetXaxis()->SetTitle("time [ns]");
            h->GetYaxis()->SetTitle("effNorm yield");
            legEffNormWave->AddEntry(h, Form("%.3f PPM", ppmFile[ifile]), "l");
            h->Draw(firstDrawn ? "hist same" : "hist");
            firstDrawn = true;
        }
        // legEffNormWave->Draw();
        canEffNormByFile->BuildLegend();
        canEffNormByFile->Print(Form("effNormByFile%i.pdf", ifile));
        fout->Append(canEffNormByFile);
    }

    /* calculate relative singletVal for file 0 */

    printf("MESSAGE singletVal by channel:\n");
    for (int ich = 0; ich < singletVal.size(); ++ich)
    {
        for (int ifile = 0; ifile < singletVal[ich].size(); ++ifile)
            printf("  chan %i file %i PPM %f singletVal %f\n", ich, ifile, ppmFile[ifile], singletVal[ich][ifile]);
    }

    fout->Write();
}