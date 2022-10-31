//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <valarray>
#include <TString.h>
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TString.h"
#include "TObjString.h"
#include "TSystemDirectory.h"
#include "TFile.h"
#include "TBWave.hxx"

TBWave *bwave;
TTree  *btree;
using namespace std;

void readEvents(ifstream *stream)
{
  char line[256];
  std::vector<Short_t> vshort;
  // cout << is << " for "  << rawEv->description  << endl;
  Long64_t nentries = 0;
  unsigned nlines = 0;
  while (stream->good())
  { 
    stream->getline(line, 256);
    TString tline(line);
    TObjArray *tokenArray = tline.Tokenize(' ');

    //if(nlines/100000*100000== nlines)
      //cout << "  " << nlines << "  " << vshort.size()  << " " << btree->GetEntries() << endl;
    if (tokenArray->GetEntries() <1)
      continue;
    TObjString *sv = (TObjString *)tokenArray->At(0);
    if (tokenArray->GetEntries() > 1)
    {
      if (sv->GetString() == TString("DC") && vshort.size() > 0)
      {
        bwave->Timestamp = ++nentries;
        //bwave->Samples->Set(vshort.size());
        for (unsigned is = 0; is < vshort.size(); ++is) 
          bwave->Samples->SetAt(vshort[is],is);
        btree->Fill();
        //cout << " end of event " << bwave->Timestamp << " samples " << vshort.size() << " , " << bwave->Samples->GetSize() << endl;
        bwave->clear();
        vshort.clear();
      }
    }
    else
    {
      vshort.push_back( Short_t(sv->GetString().Atoi()) );
      //if(vshort.size()<10) cout << " filling " << vshort.size() <<  "  " <<  Short_t(sv->GetString().Atoi()) << endl;
    }
    ++nlines;
    if(nentries>10000)
      break;
  }
  printf (" readEvents for %s digits  %lld \n",bwave->GetName(),btree->GetEntries());
}


unsigned openFiles(TString tag) 
{
  TString dirName= TString("data/")+tag;
  cout << "dirName " << dirName << endl;
  TSystemDirectory dir("rawdir",dirName);
  TList *files = dir.GetListOfFiles();

  TIter next(files);
  TSystemFile *file;
  unsigned idet = 0;
  while ((file = (TSystemFile *)next()))
  {
    string name = string(file->GetName());
    string exten  = name.substr( name.find_last_of(".")+1 );
    if(exten!=string("txt")) continue;
    string tag = name.substr( 0, name.find_last_of(".") );
    string fullName =  string( dirName.Data())  + string("/")+name;
    ifstream* in = new ifstream(fullName,std::ios::in);
    if(! in->is_open())
      continue;

    TFile *fout = new TFile(Form("%s/det%u.root", dirName.Data(), idet), "recreate");
    bwave = new TBWave(Form("Data-det%i",idet));
    btree = new TTree("Data_R","Data_R");
    //bwave->makeBranches(btree);
    cout << fullName << " run name " << bwave->GetName() << endl;
    readEvents(in);
    in->close();
    fout->ls();
    fout->Write();
    fout->Close();
    ++idet;
  }
  printf("opened %u files \n",idet);
  return idet;
}

int main(int argc, char* argv[])
{
  cout << "executing " << argv[0] << endl;
  if(argc<1) {
    printf(" usage: convertRaw <tag>  \n ");
    exit(0);
   }
  cout << argv[1] << endl;
  TString tag(argv[1]);

  if (openFiles(tag) < 1)
    return 1;

  
  exit(0);
}
