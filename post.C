TBRun *tbrun;

void post(char *fname = "anaRun-No_Doping_1_Digitizer-20000.root")
{
    gStyle->SetOptStat(1001101);
    TString fileName(fname);
    TString fullName = TString("compiled/") + fileName;
    TFile *fin = new TFile(fullName, "readonly");
    if (!fin)
    {
        cout << " didnt find " << fileName << endl;
        return;
    }
    cout << " opened  " << fileName << endl;
    /* get tbrun */
    tbrun = NULL;
    fin->GetObject("TBrun", tbrun);
    if (!tbrun)
    {
        printf("TBRun not found \n");
        return;
    }

    Long64_t ntriggers = tbrun->nevents;
    printf(" total triggers analyzed %llu \n", ntriggers);
    /* get histos from file */
    vector<TH1D *> hsum;
    TIter next(fin->GetListOfKeys());
    TKey *key;
    while (TKey *key = (TKey *)next())
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1D"))
            continue;
        TH1D *h = (TH1D *)key->ReadObj();
        TString hname(h->GetName());
        if (!hname.Contains("sumW"))
            continue;
        hsum.push_back(h);
    }

    /*normalize */
    vector<TH1D *> hscale;
    hscale.resize(hsum.size());
    for (unsigned i = 0; i < hsum.size(); ++i)
    {
        hscale[i] = (TH1D *)hsum[i]->Clone(Form("Scale%s", hsum[i]->GetName()));
        hscale[i]->SetTitle(Form("Scaled %s", hsum[i]->GetTitle()));
        hscale[i]->Sumw2(kFALSE);
        hscale[i]->Scale(1. / double(ntriggers)); //"width"
        Long64_t entries = hsum[i]->GetEntries();
        double inte = hsum[i]->Integral();
        printf(" det %u entries %llu integral %.0f ", i, entries, inte);
        entries = hscale[i]->GetEntries();
        inte = hscale[i]->Integral();
        printf(" det %u entries %llu integral %.0f \n", i, entries, inte);
    }

    vector<int> hcolor;
    hcolor.resize(hsum.size());

    /** find max val */
    double ymax;
    for (unsigned i = 4; i < 8; ++i)
    {
        ymax = TMath::Max(ymax, hscale[i]->GetBinContent(hscale[i]->GetMaximumBin()));
    }
    ymax *= 1.1;

    TCanvas *can0 = new TCanvas("summed0", "summed0");
    can0->Divide(2, 2);
    can0->SetLogy();

    for (unsigned i = 0; i < 8; ++i)
    {
        hcolor[i] = i + 1;
        hscale[i]->SetMarkerColor(hcolor[i]);
        hscale[i]->SetLineColor(hcolor[i]);
        hscale[i]->SetMarkerStyle(7);
        hscale[i]->SetMarkerSize(2);
        hscale[i]->GetXaxis()->SetRangeUser(0, 5000.);
        hscale[i]->GetYaxis()->SetRangeUser(0.1, ymax);
    }
    hscale[3]->SetMarkerStyle(5);
    hscale[3]->SetMarkerSize(.2);
    hscale[4]->SetMarkerColor(kGreen + 3);

    for (unsigned i = 0; i < 4; ++i)
    {
        can0->cd(i + 1);
        gPad->SetLogy();
        hscale[i]->Draw();
    }

    TCanvas *can1 = new TCanvas("summed1", "summed1");
    can1->Divide(2, 2);
    can1->SetLogy();
    for (unsigned i = 4; i < 8; ++i)
    {
        can1->cd(i + 1 - 4);
        gPad->SetLogy();
        hscale[i]->Draw();
    }

    TCanvas *canall = new TCanvas("summed-all", "summed-all");
    canall->SetLogy();
    for (unsigned i = 0; i < hsum.size(); ++i)
    {
        if (i == 0)
            hscale[i]->Draw("");
        else
            hscale[i]->Draw("sames");
    }
}