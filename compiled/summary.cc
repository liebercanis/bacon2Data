#include <TString.h>
#include <TROOT.h>
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TGraphErrors.h"
#include "TList.h"
//
#include "anaRun.cc"

std::vector<TString> fileList;
std::vector<double> filenum;
std::vector<double> efilenum;
std::vector<double> vqsum8;
std::vector<double> veqsum8;


// count subruns and channels
unsigned long  countFiles()
{
  TString dirName = TString("myData/");
  TSystemDirectory dir("mydata", dirName); // TSystemDirectory
  TList* files = dir.GetListOfFiles();            //
  TIter next(files);
  TSystemFile *file;
  while ((file = (TSystemFile *)next()))
  {
    string name = string(file->GetName());
    cout << name << endl;
    string exten = name.substr(name.find_last_of(".") + 1);
    if (exten != string("root"))
      continue;
    if (name.find("anaRun-run") != std::string::npos)
      fileList.push_back(TString(name.c_str()));
  }
  return fileList.size();
}

int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << " make summary plots  " << endl;
  printf(" usage: summary  <max files 0=all>  \n ");
  if(argc<1)
    exit(0);
  TString tag("run");

  TFile *fout = new TFile("summary.root","recreate");

  countFiles();
  printf(" for %s found %lu files \n", tag.Data(), fileList.size());
  Long64_t maxFiles = fileList.size();
  if (argc > 1)
  {
    maxFiles = atoi(argv[2]);
  }

  cout << " starting summary for  " << maxFiles << " files " << endl;


  for (unsigned ifile = 0; ifile < maxFiles; ++ifile){
    cout << " starting anaRunFile " << fileList[ifile] << endl;
    TString fullName = TString("myData/") + fileList[ifile];
    TFile *fin = new TFile(fullName);
    TBFile *bf;
    fin->GetObject("tbfile", bf);
    if(!bf) {
      //printf(" no timestamp in file %s \n",fileList[i].Data());
      continue;
    }
    cout << bf->GetTitle()  << " modified " << bf->modified << endl;
    filenum.push_back(double(ifile));
    efilenum.push_back(0);
    // use bf->modified a std string
    /* migrate fit from post.C:        TF1 *g = (TF1*)hHitSum[i]->GetListOfFunctions()->FindObject("expo");*/

    TH1D *hqsum;
    TH1D *hqprompt;
    fin->GetObject("histQsum", hqsum);
    fin->GetObject("histQprompt", hqprompt);

    if (!hqsum || !hqprompt)
      continue;

    cout << hqsum->GetName() << " " << hqprompt->GetName() << endl;

    for (int i = 0; i < hqsum->GetNbinsX(); ++ i){
      if(i==4)
        continue;
      cout << "chan " << i + 1 << " qsum " << hqsum->GetBinContent(i + 1) << endl;
      if(i+1==8){
        vqsum8.push_back(hqsum->GetBinContent(i + 1));
        veqsum8.push_back(hqsum->GetBinError(i + 1));
      }
    }
    TGraphErrors *gqsum = new TGraphErrors(filenum.size(), &filenum[0], &vqsum8[0], &efilenum[0], &veqsum8[0]);
    TCanvas *can = new TCanvas("Qsummary","Qsummary");
    gqsum->Draw("ap");
    fout->Append(gqsum);
    fout->Append(can);
  }

  cout << "summary finished " << maxFiles << endl;

  fout->ls();
  fout->Write();
  exit(0);
}
