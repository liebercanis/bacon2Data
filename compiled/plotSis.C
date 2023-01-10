void plotSis(TString file = "sis3316_test_data_2000s.root")
{
    TFile *fin = new TFile(file, "readonly");
    fin->ls();

    TTree *tree;
    fin->GetObject("raw", tree);
    printf("tree has %lld entries \n", tree->GetEntries());

    //
    TIter next(fin->GetListOfKeys());
    TKey *key;
    TCanvas *can;
    while (TKey *key = (TKey *)next())
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1"))
            continue;
        TH1I *h = (TH1I *)key->ReadObj();
        TString hname(h->GetName());
        cout << hname << endl;
        TString canName;
        canName.Form("plot-%s", hname.Data());
        can=new TCanvas(canName, canName);
        can->SetLogy();
        h->Draw();
        can->Print(".pdf");
    }

}