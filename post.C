TTree *btree;
vector<TH1D *> hsum;
vector<TH1D *> hFFT;

void post(TString fileName = TString("anaRun-run-01_12_2023-nev0-100.root"))
{
    gStyle->SetOptStat(1001101);
    TString fullName = TString("myData/") + fileName;
    TFile *fin = new TFile(fullName, "readonly");
    if (!fin)
    {
        cout << " didnt find " << fileName << endl;
        return;
    }
    cout << " opened  " << fileName << endl;
    /* get tbrun */
    btree = NULL;
    fin->GetObject("RunTree", btree);
    if (!btree)
    {
        printf("TBRun not found \n");
        return;
    }

    Long64_t ntriggers = btree->GetEntriesFast();
    printf(" total triggers analyzed %llu \n", ntriggers);
    /* get sum histos from file */
    TDirectory *sumDir;
    fin->GetObject("sumDir",sumDir);
    TIter next(sumDir->GetListOfKeys());
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

    /* get FFT histos from file */
    TDirectory *fftDir;
    fin->GetObject("fftDir",fftDir);
    TIter nextf(fftDir->GetListOfKeys());
    while (TKey *key = (TKey *)nextf())
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1D"))
            continue;
        TH1D *h = (TH1D *)key->ReadObj();
        TString hname(h->GetName());
        cout << (h->GetName()) << endl;
        hFFT.push_back(h);
    }


    /*normalize 
    vector<TH1D *> hscale;
    hscale.resize(hsum.size());
    for (unsigned i = 0; i < hsum.size(); ++i)
    {
        hscale[i] = (TH1D *)hsum[i]->Clone(Form("Scale%s", hsum[i]->GetName()));
        hscale[i]->SetTitle(Form("Scaled %s", hsum[i]->GetTitle()));
        hscale[i]->Sumw2(kFALSE);
        hscale[i]->Scale(1. / double(ntriggers)); //"width"
        hscale[i]->Rebin(20);
        Long64_t entries = hsum[i]->GetEntries();
        double inte = hsum[i]->Integral();
        printf(" det %u entries %llu integral %.0f ", i, entries, inte);
        entries = hscale[i]->GetEntries();
        inte = hscale[i]->Integral();
        printf(" det %u entries %llu integral %.0f \n", i, entries, inte);
    }

    vector<int> hcolor;
    hcolor.resize(hsum.size());
    */
}

TCanvas *plot1(int iev = 0, int ichan = 12)
{
  TString tchan; 
  TString tev;
  TString thit;
  TString tder;
  TString tcross;
  tchan.Form("chan%i",ichan);
  tev.Form("EvWave%i_",iev);
  thit.Form("EvHitWave%i_",iev);
  tder.Form("EvDerWave%i_",iev);
  tcross.Form("EvCross%i_",iev);


  TH1D* hwave=NULL;
  TH1D* hhit=NULL;
  TH1D* hder=NULL;
  TH1D *hcross = NULL;
  for (unsigned i = 0; i < hFFT.size(); ++i)
  {
        TString tname = TString(hFFT[i]->GetName());
        if (tname.Contains(tchan) && tname.Contains(tev))
            hwave = hFFT[i];
        if (tname.Contains(tchan) && tname.Contains(thit))
            hhit = hFFT[i];
        if (tname.Contains(tchan) && tname.Contains(tder))
            hder = hFFT[i];
        if (tname.Contains(tchan) && tname.Contains(tcross))
            hcross = hFFT[i];
  }

  if(hwave) cout << "got " << hwave->GetName() << endl ;
  else cout << tev << " not found " << endl;

  if(hhit) cout << "got "  << hhit->GetName() << endl ;
  else cout << thit << " not found " << endl;

  if(hder) cout << "got "  << hder->GetName() << endl ;
  else cout << tder << " not found " << endl;

  if (hcross)
        cout << "got " << hcross->GetName() << endl;
  else
        cout << tcross << " not found " << endl;

  if(!(hwave&&hhit&&hder)) return NULL;

  TString canName;
  canName.Form("Chan-%i-Event-%i",ichan,iev);
  TCanvas *can1 = new TCanvas(canName,canName);
  can1->Divide(1,2);
  can1->cd(1);
  hwave->Draw();
  hhit->SetFillColor(kRed);
  hhit->Draw("same");
  can1->cd(2);
  hcross->SetFillColor(kRed);
  hder->Draw();
  hcross->Draw("same");
  return can1;
}

void summed()
{
    TCanvas *canAll = new TCanvas("summed-all", "summed-all");
    canAll->SetLogy();
    canAll->Divide(4,3);
    for (unsigned i = 0; i < hsum.size(); ++i)
    {
      canAll->cd(i+1);
      hsum[i]->Draw("");
    }
    //canall->BuildLegend();

    TCanvas *canFFT = new TCanvas("FFT-all", "FFT-all");
    canFFT->Divide(4,3);
    for (unsigned i = 0; i < hFFT.size(); ++i)
    {
      canFFT->cd(i+1);
      hFFT[i]->Draw("");
    }
}

void multiPlot(int max=20) {

    for (int iev=0 ; iev < max ; ++iev){
        for (int ichan = 0; ichan<13; ++ichan){
            TCanvas *can = plot1(ichan, iev);
            if(can)
                can->Print(".pdf");
        }
    }
}