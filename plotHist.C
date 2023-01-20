////////////////////////////////////////////////////////
TString tag;
TFile *fin;
TDirectory *fftDir;
TCanvas *can;
TH1D*  hwave;
TH1D* hhit;
int currentEvent;

void openFile() {
  // open ouput file and make some histograms
  //TString fileName; fileName.Form("rootData/anaUn-%s.root",tag.Data());
  TString fileName("compiled/anaRun-No_Doping_1_Digitizer-10.root");
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
  waveName.Form("EvWave%i_chan%i",ievent,ichan);
  TString hitName;
  hitName.Form("EvHitWave%i_chan%i",ievent,ichan);
  printf(" looking for %s %s \n",waveName.Data(),hitName.Data());

  hwave=NULL;
  hhit=NULL;
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
  }
  if(hwave)
    cout << hwave->GetName() << endl;
  if (hhit)
    cout << hhit->GetName() << endl;
  return hwave != NULL && hhit != NULL;
}

void plote(int ievent,int ichan){
  if(!getHists(ievent,ichan))
    return;
  TString cname;
  cname.Form("Event%i-Chan%i", ievent, ichan);
  TCanvas *cane = new TCanvas(cname, cname);
  hwave->Draw();
  hhit->SetLineColor(kRed);
  hhit->SetFillColor(kRed);
  hhit->SetFillStyle(3002);
  hhit->Draw("same");
}

void next(int ichan=4) {
  plote(++currentEvent,ichan);
}

void plotEv(TString rtag = "Test") {
  currentEvent = -1;
  tag = rtag;
  openFile();
}
