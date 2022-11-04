#include "TBRun.hxx"
ClassImp(TBRun)

TBRun::TBRun(TString runName ): TNamed(runName,runName)
{
  btree = new TTree("BRunTree"," bacon data " );
  bevent = new TBEvent(runName);
  btree->Branch("bev",&bevent);
  det0.SetName(Form("wave%i",1));
  det1.SetName(Form("wave%i",4));
  det2.SetName(Form("wave%i",2));
  det3.SetName(Form("wave%i",6));
  detList.push_back(&det0);
  detList.push_back(&det1);
  detList.push_back(&det2);
  detList.push_back(&det3);

  for(unsigned i=0; i<NDET; ++i ) {
    cout << " btree adding branch " << detList[i]->GetName() << endl;
    btree->Branch(detList[i]->GetName(),detList[i]);
  }
  cout << " TBRun tree " << btree->GetName() << endl;
  btree->GetListOfBranches()->ls();
  clear();
}

// get det names from pre-existing btree
TBRun::TBRun(TTree* btreeInit, TString runName ): TNamed(runName,runName)
{
  TString brunName = TString("BRun"); 
  btree = new TTree("BRunTree","bacon data" );
  bevent = new TBEvent(runName);
  btree->Branch("bev",&bevent);
  det0.SetName(Form("wave%i",1));
  det1.SetName(Form("wave%i",4));
  det2.SetName(Form("wave%i",2));
  det3.SetName(Form("wave%i",6));
  detList.push_back(&det0);
  detList.push_back(&det1);
  detList.push_back(&det2);
  detList.push_back(&det3);


  TObjArray *brList = btreeInit->GetListOfBranches();
  brList->ls();
  TIter next(brList);
  TBranch *aBranch=NULL;
  int idet=0;
  while( ( aBranch = (TBranch *) next() ) ) {
    if(!TString(aBranch->GetName()).Contains("wave") ) continue;
    detList[idet]->SetName(aBranch->GetName());
    ++idet;
  }

  for(unsigned i=0; i<NDET; ++i ) {
    cout << " btree adding branch " << detList[i]->GetName() << endl;
    btree->Branch(detList[i]->GetName(),detList[i]);
  }
  cout << " TBRun tree " << btree->GetName() << " detList size " << detList.size() << endl;
  btree->GetListOfBranches()->ls();
  clear();


}


//TBRun::~TBRun(){}

void TBRun::clear()
{
  bevent->clear();
  detListClear();
}

void TBRun::detListClear() 
{
  for(unsigned i=0; i<NDET; ++i ) detList[i]->clear(); 
}

void TBRun::dumpEvent(ofstream& dumpFile) {
  if (!dumpFile.is_open()||!dumpFile.good()) {
    cout  << " bad dumpFile " << endl;
    return;
  }
  // event header lines
  dumpFile  << "C**** dump run " << bevent->GetName()  << " event " <<  bevent->event << endl;
  dumpFile << bevent->event << " " << bevent->trigTime << endl;

  for(unsigned id=0; id<detList.size(); ++id)  {
    dumpFile << id << " " <<  detList[id]->GetName() << " "  << detList[id]->hits.size() << " " << detList[id]->nspe << "  " << detList[id]->qSum  << endl;
    for(unsigned ih=0; ih<detList[id]->hits.size(); ++ih)  {
      // these are all in a row 
      dumpFile << ih <<  " " 
        << detList[id]->hits[ih].firstBin << " " 
        << detList[id]->hits[ih].lastBin << " " 
        << detList[id]->hits[ih].startTime << " "
        << detList[id]->hits[ih].peakWidth << " "
        << detList[id]->hits[ih].qpeak << " "
        << detList[id]->hits[ih].peakt << " "
        << detList[id]->hits[ih].peakMaxTime << " "
        << detList[id]->hits[ih].peakBin << " "
        << detList[id]->hits[ih].qsum << " "
        << detList[id]->hits[ih].qerr << " "
        << detList[id]->hits[ih].good << " "
        << detList[id]->hits[ih].kind << " "
        << endl;
    }
  }
}



