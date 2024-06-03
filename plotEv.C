////////////////////////////////////////////////////////
TString tag;
TFile *fin;
TDirectory *sumWaveDir;
TCanvas *can;
TH1D*  hwave;
TH1D* hhit;
TH1D* hder;
int currentEvent;
vector<TH1D*> vwave;
vector<TH1D*> vhit;


void plot1(TH1D* hwave,TH1D* hhit){
  TString cname;
  cname.Form("PLot%s",hwave->GetName());
  TCanvas *cane = new TCanvas(cname, cname);
  cout << "  canvas " << cname << endl;
  hwave->Draw();
  hhit->SetLineColor(kRed);
  hhit->SetFillColor(kRed);
  hhit->SetFillStyle(3002);
  //hder->SetLineColor(kGreen);
  hhit->Draw("same");
  //hder->Draw("same");
  cane->Print(".pdf");
}


void openFile() {
  // open ouput file and make some histograms
  //TString fileName; fileName.Form("rootData/anaUn-%s.root",tag.Data());
  TString fileName("caenData/anaCRun-run-12_28_2023-file-save.root");
  printf(" ***** looking for file %s  ***** \n",fileName.Data());
  fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }
  fin->GetObject("sumWaveDir",sumWaveDir);
  //sumWaveDir->ls();
}

void getAll()
{

  TIter next(sumWaveDir->GetListOfKeys());
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;

    TH1D *h = (TH1D *)key->ReadObj();
    TString hname(h->GetName());
    int ichan = TString(hname(hname.Last('n') + 1, hname.Length())).Atoi();
    int ievent = TString(hname(hname.Last('e') + 1, hname.Length())).Atoi();

    printf(" ev %i chan %i \n",ievent,ichan);

    TString waveName;
    waveName.Form("EvWave%ichan%i",ievent,ichan);
    TString hitName;
    hitName.Form("EvHitPeakWave%ichan%i",ievent,ichan);
    if (hname == waveName)
      vwave.push_back(h);
    if(hname == hitName)
      vhit.push_back(h);
  }

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
  TIter next(sumWaveDir->GetListOfKeys());
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
  /*
  if(hwave)
    cout << hwave->GetName() << endl;
  if (hhit)
    cout << hhit->GetName() << endl;
  if (hhit)
    cout << hder->GetName() << endl;
  */

 
  return hwave != NULL && hhit != NULL;
}

void plote(int ievent,int ichan){
  if(!getHists(ievent,ichan)){
    cout << " no hist for event " << ievent << "  " << ichan << endl;
    return;
  }
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
  getAll();
  for(unsigned i=0; i< vhit.size(); ++i ) {
    printf("%u %s %s \n",i, vwave[i]->GetName(),vhit[i]->GetName());
    plot1(vwave[i],vhit[i]);
    if(i>10) break;
  }
}
