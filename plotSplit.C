TTree *btree;
vector<TH1D *> hwave;
vector<TH1D *> hsplit;
TString tag;


TCanvas *plot1(unsigned iev=0)
{
  TString tchan; 
  if(hwave.size()<iev) {
    printf(" %u max is %lu\n", iev, hwave.size());
    return NULL;
  }

  TString canName;
  canName.Form("Split-Wave-%i",iev);
  TCanvas *can1 = new TCanvas(canName,canName);
  hwave[iev]->Draw();
  hsplit[iev]->SetLineColor(kRed);
  hsplit[iev]->SetFillColor(kRed);
  hsplit[iev]->Draw("same");
  return can1;
}

void plotSplit(TString fileName = TString("anaRun-run-01_12_2023-file0-100.root"))
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
    TDirectory *splitDir;
    fin->GetObject("splitDir",splitDir);
    TIter next(splitDir->GetListOfKeys());
    TKey *key;
    while (TKey *key = (TKey *)next())
    {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1D"))
        continue;
      TH1D *h = (TH1D *)key->ReadObj();
      TString hname(h->GetName());
      if (hname.Contains("EvSplitWave")) {
        hwave.push_back(h);
      } else if (hname.Contains("EvSplitPeak")) {
        hsplit.push_back(h);
      }
    }

    cout << " number split events " << hsplit.size() << endl;

  }
