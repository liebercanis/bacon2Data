
/*
    pickup graphs from postMacroFile.C
*/
#include "TMultiGraph.h"
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
                absorb[ifile].push_back(Absorbtion(ppmFile[ifile], distanceLevel[ilevel]));
            }
        }
    }

    for (unsigned ifile = 0; ifile < absorb.size(); ++ifile)
    {
        for (unsigned il = 0; il < absorb[ifile].size(); ++il)
            printf("file %u level %u absorb %E\n", ifile, il, absorb[ifile][il]);
    }

    // for (unsigned ifile = 0; ifile < ppmFile.size(); ++ifile)
    //   tbDraw(-1, ppmFile[ifile]);

    fout = new TFile("integralGraphs.root", "recreate");
    for (int ifile = 0; ifile < 11; ++ifile)
    {
        fileName.Form("postMacroAllFile%i.root", ifile);
        TFile *fin = new TFile(fileName, "readonly");
        TString graphName;
        graphName.Form("gSingletIntegralsFile%i", ifile);
        fin->GetObject(graphName, gr);
        gsinglet.push_back(gr);
        fout->Add(gr);
        cout << " got graph " << gr->GetName() << endl;
        graphName.Form("gLateIntegralsFile%i", ifile);
        fin->GetObject(graphName, gr);
        cout << " got graph " << gr->GetName() << endl;
        fout->Add(gr);
        glate.push_back(gr);
        fin->Close();
    }
    // fout->ls();

    // each one of these corresponds to a file
    for (unsigned i = 0; i < gsinglet.size(); ++i)
        printf("%f %s %s \n", ppmFile[i], gsinglet[i]->GetName(), glate[i]->GetName());

    // get points from graphs

    singletVal.resize(gsinglet[0]->GetN());        // number of channels
    for (unsigned i = 0; i < gsinglet.size(); ++i) // loop over all files
    {
        Double_t x, y;
        for (int ich = 0; ich < gsinglet[i]->GetN(); ++ich) // loop over channels
        {
            gsinglet[i]->GetPoint(ich, x, y);
            singletVal[ich].push_back(y);
            // printf("line89 file %i ich %i PPM %f point %lu singlet %f lateb%f\n", i, ich, ppmFile[i], singletVal[ich].size() - 1, singletVal[ich].back(), singletVal[ich].back());
        }
        // for() printf("line 92 file %i point %lu  singletVal  %f\n", i, singletVal.size());
    }
    printf("line95 singlet val size %lu\n", singletVal.size());

    lateVal.resize(glate[0]->GetN());           // number of channels
    for (unsigned i = 0; i < glate.size(); ++i) // loop over all files
    {
        Double_t x, y;
        for (int ich = 0; ich < glate[i]->GetN(); ++ich) // loop over channels
        {
            glate[i]->GetPoint(ich, x, y);
            lateVal[ich].push_back(y);
            printf("file %i ich %i PPM %f late %f lateb%f\n", i, ich, ppmFile[i], lateVal[ich].back(), lateVal[ich].back());
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
    std::vector<vector<double>> modelSingletVal;
    modelSingletVal.resize(singletVal.size()); // first index channel
    printf("modelSingletVal size %lu singletVal %lu\n", modelSingletVal.size(), singletVal.size());
    for (int ich = 0; ich < singletVal.size(); ++ich) // loop over channels
    {
        // printf("line 139 ch %i %lu \n", ich, singletVal[ich].size());
        int ilevel = getLevel(ich);
        if (ilevel >= NLEVELS)
            break;
        for (int ifile = 0; ifile < singletVal[ich].size(); ++ifile)
        { // loop over files
            // printf("line 143 ch %i file %i size %lu single %f abs %f \n", ich, ifile, modelSingletVal[ich].size(), singletVal[ich][0], absorb[ifile][ilevel]);
            modelSingletVal[ich].push_back(singletVal[ich][4] * (1 - absorb[ifile][ilevel]));
            if (ifile == 0)
                printf("line 123 modelSingletVal ch %i ilevel %i ppm %.3f singletVal %E size %lu  \n", ich, ilevel, ppmFile[ifile], modelSingletVal[ich].back(), modelSingletVal[ich].size());
        }
    }
    std::vector<vector<double>> modelLateVal;
    modelLateVal.resize(lateVal.size()); // first index channel
    /*
        late integral has to be compared do model no (1-A)
    */
    for (int ich = 0; ich < lateVal.size(); ++ich) // loop over channels
    {
        int ilevel = getLevel(ich);
        if (ilevel >= NLEVELS)
            break;
        for (int ifile = 0; ifile < lateVal[ich].size(); ++ifile)
        { // loop over files
            modelLateVal[ich].push_back(lateVal[ich][4] * (1 - absorb[ifile][ilevel]));
            if (ifile == 0)
                printf("line 144 modelTripletVal ch %i ilevel %i ppm %.3f late val %E size %lu  \n", ich, ilevel, ppmFile[ifile], modelLateVal[ich].back(), modelLateVal[ich].size());
        }
    }

    /* singlet model graph */
    TGraph *modelSinglet = nullptr;
    int oldLevel = -1;
    for (int ich = 0; ich < singletVal.size(); ++ich)
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
        modelSinglet = new TGraph(singletVal[ich].size(), &ppmFile[0], &modelSingletVal[ich][0]);
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
    for (int ich = 0; ich < singletVal.size(); ++ich)
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
        modelLate = new TGraph(lateVal[ich].size(), &ppmFile[0], &modelLateVal[ich][0]);
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

    for (int ich = 0; ich < singletVal.size(); ++ich)
    {
        int ilevel = getLevel(ich);
        if (ilevel > 3)
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

    for (int ich = 0; ich < lateVal.size(); ++ich)
    {
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