/* pickup graphs from postMacroFile.C*/
TFile *fout;
TString fileName;
TGraph *gr;
vector<TGraph *> gsinglet;
vector<TGraph *> glate;
vector<double> ppmFile;
void makeIntegralGraphs()
{
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
    vector<vector<double>> singletVal;
    vector<vector<double>> lateVal;
    singletVal.resize(gsinglet[0]->GetN()); // number of channels
    lateVal.resize(glate[0]->GetN());       // number of channels

    for (unsigned i = 0; i < gsinglet.size(); ++i) // loop over all files
    {
        Double_t x, y;
        for (int ich = 0; ich < gsinglet[i]->GetN(); ++ich) // loop over channels
        {
            gsinglet[i]->GetPoint(ich, x, y);
            singletVal[ich].push_back(y);
            glate[i]->GetPoint(ich, x, y);
            lateVal[ich].push_back(y);
            printf("file %i ich %i PPM %f singlet %f lateb%f\n", i, ich, ppmFile[i], singletVal[ich].back(), lateVal[ich].back());
        }
    }

    // make graphs for each channel
    TGraph *singletPPM;
    for (int ich = 0; ich < singletVal.size(); ++ich)
    {
        singletPPM = new TGraph(singletVal.size(), &ppmFile[0], &singletVal[ich][0]);
        singletPPM->SetName(Form("singletPPMChan%i", ich));
        singletPPM->SetTitle(Form("singletPPMChan%i", ich));
        singletPPM->SetMarkerStyle(21);
        fout->Add(singletPPM);
        TCanvas *csinglet = new TCanvas(Form("canSingletPPMChan%i", ich), "canSingletPPM");
        csinglet->SetName(Form("canSingletPPMChan%i", ich));
        csinglet->SetTitle(Form("canSingletPPMChan%i", ich));
        singletPPM->SetMarkerStyle(21);
        singletPPM->GetHistogram()->GetXaxis()->SetTitle("PPM");
        singletPPM->GetHistogram()->GetYaxis()->SetTitle("singlet integral");
        csinglet->SetGrid();
        singletPPM->Draw("ap");
        csinglet->Print(".pdf");
        fout->Append(csinglet);
    }

    // make graphs for each channel
    TGraph *latePPM;
    for (int ich = 0; ich < lateVal.size(); ++ich)
    {
        latePPM = new TGraph(lateVal.size(), &ppmFile[0], &lateVal[ich][0]);
        latePPM->SetName(Form("canLatePPMChan%i", ich));
        latePPM->SetTitle(Form("canLatePPMChan%i", ich));
        latePPM->GetHistogram()->GetXaxis()->SetTitle("PPM");
        latePPM->GetHistogram()->GetYaxis()->SetTitle("late integral");
        latePPM->SetMarkerStyle(21);
        fout->Add(latePPM);

        TCanvas *clate = new TCanvas(Form("canlatePPMChan%i", ich), "canlatePPM");
        clate->SetTitle(Form("canLatePPMChan%i", ich));
        latePPM->SetMarkerStyle(21);
        latePPM->GetHistogram()->GetXaxis()->SetTitle("PPM");
        latePPM->GetHistogram()->GetYaxis()->SetTitle("late integral");
        clate->SetGrid();
        latePPM->Draw("ap");
        clate->Print(".pdf");
        fout->Append(clate);
    }
    fout->Write();
}