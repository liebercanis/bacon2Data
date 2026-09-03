
/*
    pickup graphs from postMacroFile.C
*/
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "modelAllFit.hh"
#include "distanceLevels.hh"
TFile *fout;
TString fileName;
TGraph *gr;
vector<TGraph *> gsinglet;
vector<TGraph *> glate;
vector<double> ppmFile;
vector<vector<double>> absorb;
vector<vector<double>> singletVal;
vector<vector<double>> modelSingletVal;
vector<vector<double>> modelLateVal;
vector<vector<double>> lateVal;

void makeIntegralGraphs()
{
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
    for (int ifile = 0; ifile < ppmFile.size(); ++ifile)
    {
        fileName.Form("postMacroAllFile%i.root", ifile);
        TFile *fin = new TFile(fileName, "readonly");
        TString graphName;
        graphName.Form("gSingletIntegralsFile%i", ifile);
        fin->GetObject(graphName, gr);
        gsinglet.push_back(gr);
        fout->Add(gr);
        cout << ifile << " got graph " << gr->GetName() << endl;
        graphName.Form("gLateIntegralsFile%i", ifile);
        fin->GetObject(graphName, gr);
        // cout << " got graph " << gr->GetName() << endl;
        fout->Add(gr);
        glate.push_back(gr);
        fin->Close();
    }
    printf("line 81 ppm file %li gsinglet %li glate %li \n", ppmFile.size(), gsinglet.size(), glate.size());
    // fout->ls();
    modelSingletVal.resize(gsinglet[0]->GetN()); // number of channels
    modelLateVal.resize(glate[0]->GetN());       // number of channels
    for (int ifile = 0; ifile < 11; ++ifile)
    {
        fileName.Form("tbModelPPM%0.3f.root", ppmFile[ifile]);
        // printf("open %s \n", fileName.Data());
        TFile *fin = new TFile(fileName, "readonly");
        for (int ichan = 0; ichan < 12; ++ichan)
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
                continue;
            }
            // TH1D *hcopy = (TH1D *)hist->Clone(Form("fitWaveFitChan%iFile%i", ichan, ifile));
            // hcopy->SetTitle(Form("fitWaveChan%iFile%i %.3fPPM", ichan, ifile, ppmFile[ifile]));
            //  integral is over bins
            int earlyBin = 1350 / 2;
            int lateBin = 1460 / 2;
            // printf("get %s \n", hist->GetName());
            double singletIntegral = hist->Integral(earlyBin, lateBin);
            double lateIntegral = hist->Integral(lateBin, 3000 / 2);
            modelSingletVal[ichan].push_back(singletIntegral);
            modelLateVal[ichan].push_back(lateIntegral);

            // printf("line 106 chan %i file %i PPM %.3f singlet integral %.3E late integral %.3E model singlet %.3E size %lu \n", ichan, ifile, ppmFile[ifile], singletIntegral, lateIntegral, modelSingletVal[ichan].back(), modelSingletVal[ichan].size());
            // int jfile = modelSingletVal[ichan].size() - 1;
            // printf("\t line 109 modelSingletVal ich %i ifile %i val %f \n", ichan, jfile, modelSingletVal[ichan][jfile]);
            //  hist->SetLineColor(ifile + 1);
            //  fout->Add(hist);
        }
        fin->Close();
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
        Double_t x, y;
        for (int ic = 0; ic < gsinglet[ifile]->GetN(); ++ic) // loop over channels
        {
            if (ic == 0)
                continue;
            gsinglet[ifile]->GetPoint(ic, x, y);
            singletVal[ic].push_back(y);
            printf("line150 file %i ic %i PPM %f point %lu singlet %f\n", ifile, ic, ppmFile[ifile], singletVal[ic].size() - 1, singletVal[ic].back());
        }
    }

    for (int ic = 0; ic < singletVal.size(); ++ic)
    { // loop over channels{
        for (int ifile = 0; ifile < singletVal[ic].size(); ++ifile)
        {
            printf("line93 singletVal ich %i ifile %i val %f\n", ic, ifile, singletVal[ic][ifile]);
        }
    }

    printf("line95 singlet val size %lu\n", singletVal.size());

    lateVal.resize(glate[0]->GetN());           // number of channels
    for (unsigned i = 0; i < glate.size(); ++i) // loop over all files
    {
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

    /*
        make model graphs
    */

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

    /* late integral has to be compared do model no (1-A) */
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

    /* singlet model graph */
    TGraph *modelSinglet = nullptr;
    int oldLevel = -1;
    for (int ich = 0; ich < gsinglet.size(); ++ich)
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
        //  printf("line 134 %i %lu %lu \n", ich, singletVal[ich].size(), modelSingletVal[ich].size());
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

    /* late model graph */
    oldLevel = -1;
    TGraph *modelLate = nullptr;
    for (int ich = 0; ich < glate.size(); ++ich)
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

    for (int ich = 0; ich < gsinglet.size(); ++ich)
    {
        if (ich == 0) // channel 0 has no data: skipped when filling singletVal
            continue;
        int ilevel = getLevel(ich);
        if (ilevel > 3)
            continue;
        printf("level %i %s \n", ilevel, mgSinglet[ilevel]->GetName());
        singletPPM = new TGraph(gsinglet.size(), &ppmFile[0], &singletVal[ich][0]);
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

    for (int ich = 0; ich < glate.size(); ++ich)
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

    fout->Write();
}