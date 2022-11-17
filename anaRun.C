#include "TBRawEvent.hxx"
#include "TBRun.hxx"
TFile *fout;
TFile *fin;
TTree *rtree;
TBRawRun *tbrun;
TObjArray *brList;
vector<TBRawEvent *> detList;
// TTree *ftree;

bool openFile(TString tag)
{
  // open ouput file and make some histograms
  TString fileName; fileName.Form("%s.root",tag.Data());
  printf(" looking for file %s\n",fileName.Data());
  fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return false;
  }

  rtree = NULL;
  brList = NULL;
  fin->GetObject("RawTree", rtree);
  if(!rtree) {
    printf("!! RawTree not found \n");
    return false;
  }

  printf(" RawTree entries = %llu\n", rtree->GetEntries() );
  brList = rtree->GetListOfBranches();
  brList->ls();

  TIter next(brList);
  TBranchElement *aBranch=NULL;
  detList.clear();

  while( ( aBranch = (TBranchElement *) next() ) ) {
    TString brname = TString(aBranch->GetName()) ;
    int chan = TString(brname(brname.First('n')+1,brname.Length() - brname.First('n'))).Atoi();
    //aBranch->GetListOfLeaves()->ls();
    //aBranch->SetAddress( &(detList[detList.size()-1]) );
    aBranch->SetAddress(0);
    //cout << " info " << aBranch->GetFullName()  << " class " << aBranch->GetClassName() << " chan " << chan << " det name  " << detList[detList.size() - 1]->GetName() << endl;
  }
  printf("det list total =  %lu \n", detList.size());
  return true;
}


void anaRun(TString tag= TString("Test"))
{
  TString outFileName = TString("anaRun") + tag + TString(".root");
  fout = new TFile(outFileName,"recreate");
  TString rawRunName = TString("raw-") + tag;
  if(!openFile(tag)) return;

  
  //rtree->GetListOfBranches()->ls();
  //printf("det list total =  %lu \n", detList.size());
  //for (unsigned i = 0; i < detList.size(); ++i)
     //printf(" %i name  %s chan %i \n ", i, detList[i]->GetName(), detList[i]->channel);
  printf(" loop over entries \n");
  for (Long64_t entry = 0; entry < 2;  ++entry)
  {
    cout << " .... entry " << entry << endl;

    rtree->GetEntry(entry);

    TIter next(brList);
    TBranchElement *aBranch = NULL;
    detList.clear();

    while ((aBranch = (TBranchElement *) next())) {
      TBRawEvent *det = (TBRawEvent *) aBranch->GetObject();
      detList.push_back(det);
    }
    

    for (int ichan = 0; ichan < detList.size(); ++ichan) 
    {
      cout << detList[ichan]->channel << " " << detList[ichan]->time << " " << detList[ichan]->rdigi.size() << "; ";
      cout << endl;
      TString hname;
      hname.Form("EvRawWaveEv%ich%i", int(entry), ichan);
      int nbins = detList[ichan]->rdigi.size();
      if(nbins<1)
        continue;
      TH1D *hEvRawWave = new TH1D(hname, hname, nbins, 0, nbins);
      for (unsigned j = 0; j < detList[ichan]->rdigi.size(); ++j)
        hEvRawWave->Fill(j + 1, (double) detList[ichan]->rdigi[j]);
    }
  }
    fout->Close();
    return;
  }
