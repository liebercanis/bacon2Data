TTree *btree;
vector<TH1D *> hsum;
vector<TH1D *> hHitSum;
vector<TH1D *> hFFT;
vector<double> slope;
vector<double> eslope;
vector<double> chan;
vector<double> echan;
TString tag;



TCanvas *plot1(int iev = 0, int ichan = 6)
{
  TString tchan; 
  TString tev;
  TString thit;
  TString tder;
  TString tcross;
  TString thitpeak;
  tchan.Form("chan%i",ichan);
  tev.Form("EvWave%i_",iev);
  thit.Form("EvHitWave%i_",iev);
  tder.Form("EvDerWave%i_",iev);
  tcross.Form("EvCross%i_",iev);
  thitpeak.Form("EvHitPeakWave%i_chan%i",iev,ichan);
  


  TString twave = tev+tchan;
  TString thwave = thit+tchan;

  TH1D* hwave=NULL;
  TH1D* hhit=NULL;
  TH1D* hder=NULL;
  TH1D *hcross = NULL;
  TH1D *hitpeak = NULL;
  for (unsigned i = 0; i < hFFT.size(); ++i)
  {
        TString tname = TString(hFFT[i]->GetName());
        if (tname==twave)
            hwave = hFFT[i];
        if (tname == thwave)
            hhit = hFFT[i];
        if (tname == thitpeak)
            hitpeak = hFFT[i];

        
        /*
        if (tname.Contains(tchan) && tname.Contains(tder))
            hder = hFFT[i];
        if (tname.Contains(tchan) && tname.Contains(tcross))
            hcross = hFFT[i];
            */
  }

  if(hwave) cout << "got " << hwave->GetName() << endl ;
  else cout << tev << " not found " << endl;
  if(hitpeak) cout << "got " << hwave->GetName() << endl ;
  else cout << tev << " not found " << endl;


  if(hhit) cout << "got "  << hhit->GetName() << endl ;
  //else cout << thit << " not found " << endl;

  if(hder) cout << "got "  << hder->GetName() << endl ;
  //else cout << tder << " not found " << endl;

  //if (hcross)
    //    cout << "got " << hcross->GetName() << endl;
  //else
    //    cout << tcross << " not found " << endl;

  if(!(hwave&&hhit)) return NULL;

  TString canName;
  canName.Form("Chan-%i-Event-%i",ichan,iev);
  TCanvas *can1 = new TCanvas(canName,canName);
  //can1->Divide(1,2);
  //can1->cd(1);
  hwave->Draw();
  hitpeak->SetLineColor(kRed);
  hitpeak->SetFillColor(kRed);
  hitpeak->Draw("same");
  //can1->cd(2);
  //hcross->SetFillColor(kRed);
  //hder->Draw();
  //hcross->Draw("same");
  return can1;
}

void summed()
{
  TCanvas *canSum = new TCanvas("wavesummed-all", "hitsummed-all");

  canSum->Divide(4, 3);
  for (unsigned i = 0; i < hsum.size(); ++i)
  {
        canSum->cd(i + 1);
        gPad->SetLogy();
        hsum[i]->Draw("");
  }
    canSum->Print(".pdf");
   // canall->BuildLegend();

   TCanvas *canAll = new TCanvas("hitsummed-all", "hitsummed-all");

  canAll->Divide(4, 3);
  for (unsigned i = 0; i < hHitSum.size(); ++i)
  {
        canAll->cd(i + 1);
        gPad->SetLogy();
        hHitSum[i]->Fit("expo","","",200,600);
        TF1 *g = (TF1*)hHitSum[i]->GetListOfFunctions()->FindObject("expo");
        printf("%s %E %E \n",hHitSum[i]->GetName(),g->GetParameter(1),g->GetParError(1));
        slope.push_back(g->GetParameter(1));
        eslope.push_back(g->GetParError(1));
        hHitSum[i]->Draw("");
    }
    canAll->Print(".pdf");
    //canall->BuildLegend();

    TCanvas *canFFT = new TCanvas("FFT-all", "FFT-all");
    canFFT->Divide(4,3);
    for (unsigned i = 0; i < hFFT.size(); ++i)
    {
      canFFT->cd(i+1);
      hFFT[i]->Draw("");
    }
    canFFT->Print(".pdf");


    // slope graph
    TGraphErrors *grslope = new TGraphErrors(chan.size()-4, &chan[3],&slope[3], &eslope[3],&echan[3]);
    TString canSlopeName;
    canSlopeName.Form("slope-%s",tag.Data());
    grslope->SetTitle(canSlopeName);
    TCanvas *canSlope = new TCanvas(canSlopeName,canSlopeName);
    grslope->SetMarkerStyle(22);
    grslope->SetMarkerSize(1.0);
    grslope->Draw("ap");
    canSlope->Print(".pdf");



}

void multiPlot(int max=20) {

    for (int iev=0 ; iev < max ; ++iev){
        for (int ichan = 0; ichan<9; ++ichan){
            TCanvas *can = plot1(iev,ichan);
            if(can)
                can->Print(".pdf");
        }
    }
}

void post(TString fileName = TString("anaRun-run-01_12_2023-file0-100.root"))
{
    tag = fileName;
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
      if (hname.Contains("sumWave")) {
        hsum.push_back(h);
      } else if (hname.Contains("sumHit")) {
        hHitSum.push_back(h);
        string sname = string(hname.Data());
        chan.push_back(stof(sname.substr(10,2))); 
        echan.push_back(0);
        //cout << hname << " " << chan[chan.size()-1] << endl;
      }

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
        //cout << (h->GetName()) << endl;
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
