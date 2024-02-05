////////////////////////////////////////////////////////
TString tag;
TFile *fin;
TDirectory *fftDir;
TCanvas *can;
TH1D*  hwave;
TH1D* hhit;
TH1D* hder;
int currentEvent;

void openFile() {
  // open ouput file and make some histograms
  //TString fileName; fileName.Form("rootData/anaUn-%s.root",tag.Data());
  TString fileName("caenData/anaCRun-run-11_26_2023-file-100.root");
  printf(" ***** looking for file %s  ***** \n",fileName.Data());
  fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }
  fin->GetObject("fftDir",fftDir);
  fftDir->ls();
}

bool getHists(int ievent, int ichan) {

  TString waveName;
  waveName.Form("EvWave%ichan%i",ievent,ichan);
  TString derWaveName;
  derWaveName.Form("EvDerWave%ichan%i",ievent,ichan);

  TString hitName;
  hitName.Form("EvHitPeakWave%ichan%i",ievent,ichan);
  printf(" looking for %s %s \n",waveName.Data(),hitName.Data());

  hwave=NULL;
  hhit=NULL;
  hder=NULL;
  TIter next(fftDir->GetListOfKeys());
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    TString hname(h->GetName());
    if (hname == waveName)
      hwave = h;
    if(hname == hitName)
      hhit = h;
    if(hname == derWaveName)
      hder = h;

  }
  if(hwave)
    cout << hwave->GetName() << endl;
  if (hhit)
    cout << hhit->GetName() << endl;
  if (hhit)
    cout << hder->GetName() << endl;

  return hwave != NULL && hhit != NULL;
}

void plote(int ievent,int ichan){
  if(!getHists(ievent,ichan))
    return;
  TString cname;
  cname.Form("Event%ichan%i", ievent, ichan);
  TCanvas *cane = new TCanvas(cname, cname);
  hwave->Draw();
  hhit->SetLineColor(kRed);
  hhit->SetFillColor(kRed);
  hhit->SetFillStyle(3002);
  hder->SetLineColor(kGreen);
  hhit->Draw("same");
  hder->Draw("same");
}

void next(int ichan=4) {
  plote(++currentEvent,ichan);
}

void plotEv(TString rtag = "Test") {
  currentEvent = -1;
  tag = rtag;
  openFile();
}
