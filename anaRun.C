vector<TBWave*> waveList;
TFile *fout;
TString tag;
TBRun* tbrun;

void plotWave( int ifile, int theRun, Long64_t jentry) 
{
  TString hname;
  hname.Form("set-%s-run%i-event-%lli",tag.Data(),theRun,jentry);
  Long64_t nsamples =  waveList[ifile]->Samples->GetSize();
  TH1S* hist = new TH1S(hname,hname,nsamples,0,nsamples);
  for(int i=0; i<nsamples; ++i) hist->SetBinContent(i+1,hist->GetBinContent(i+1) + waveList[ifile]->Samples->GetAt(i));
  TCanvas *can = new TCanvas(hname,hname);
  hist->Draw();
  can->Print(".gif");

}

void makeHits(ifile) 
{
   tbrun->detList[ifile]->clear();
   vector<double> wvec;
   for(int i=0; i<nsamples; ++i)  wvec.push_back(waveList[ifile]->Samples->GetAt(i));

}

void anaRun(int theRun=1 , const char* theTag="test")
{
  fout = new TFile("anaRun.root","recreate");
  TString dirName;
  tag = TString(theTag);
  dirName.Form("SiPM_Data_UW_Board/%s/DAQ/run_%i/RAW",tag.Data(),theRun);
  cout << "dirName " << dirName << endl;
  TSystemDirectory dir("RAW",dirName);
  TList *files = dir.GetListOfFiles();
  tbrun = new TBRun(tag);
  cout << " create TBRun  " << tbrun->GetName() << endl;

  TIter next(files);
  TSystemFile *file;
  while( (file = (TSystemFile*) next()) ) {
    string name = string(file->GetName());
    string exten  = name.substr( name.find_last_of(".")+1 );
    if(exten!=string("root")) continue;
    string tag = name.substr( 0, name.find_last_of(".") );
    string chan = tag.substr( tag.find_first_of("CH"), tag.find_first_of("@") - tag.find_first_of("CH"));
    cout << " tag " << tag.find_first_of("CH") <<  "   " <<tag << " can  " << chan << endl;
    string fullName =  string( dirName.Data())  + string("/")+name;
    cout << " open " << fullName << endl;
    TFile *f = new TFile(fullName.c_str(),"readonly");
    if(!f) continue;
    TTree *dtree=NULL;
    f->GetObject("Data_R",dtree);
    if(!dtree) continue;
    TBranch *b_Samples = dtree->FindBranch("Samples");
    if(!b_Samples) continue;
    if(b_Samples->GetEntries()<1) continue;
    cout << f->GetName()  << "  " << b_Samples->GetEntries() << endl;
    TBWave* wave = new TBWave(chan.c_str());
    TTree *tree = NULL;
    f->GetObject("Data_R",tree);
    wave->Init(tree);
    waveList.push_back(wave);
  }
  

  cout << " files = " << waveList.size() << endl;
  for(unsigned ifile = 0; ifile<waveList.size(); ++ ifile) {
    tbrun->detList[ifile]->SetName(Form("SIPM-%s", waveList[ifile]->GetName()));
    cout << " for wave  " << waveList[ifile]->GetName() << " det " <<  tbrun->detList[ifile]->GetName() << endl;
    Long64_t nbytes = 0, nb = 0;
    Long64_t ientry = waveList[ifile]->LoadTree(0);
    nb = waveList[ifile]->fChain->GetEntry(0);  
    Long64_t nsamples =  waveList[ifile]->Samples->GetSize();
    cout << ientry << "  " << nsamples << endl;

    fout->cd();

    Long64_t nentries =  waveList[ifile]->fChain->GetEntries();
    cout << " entries " << nentries << endl;
    nentries = 3;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      ientry = waveList[ifile]->LoadTree(jentry);
      if (ientry < 0) break;
      nb = waveList[ifile]->GetEntry(jentry);  
      nbytes += nb;
      nsamples =  waveList[ifile]->Samples->GetSize();
      plotWave(ifile,theRun,jentry);
      makeHits(ifile);
    }

  }

  fout->ls();
  return;
}
