#include <TString.h>
#include <TROOT.h>
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TList.h"
//
#include "anaRun.cc"

std::vector<TString> fileList;

// count subruns and channels
unsigned long  countFiles(TString dateTag)
{
  TString dirName = TString("rootData/");
  TSystemDirectory dir("rootdir", dirName); // TSystemDirectory
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
    if (name.find(dateTag.Data()) != std::string::npos)
      fileList.push_back(TString(name.c_str()));
  }
  return fileList.size();
}

int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << " do set of files by date tag as 01_12_2023 " << endl;
  printf(" usage: ana  <xx_xx_xxxx> data  <max files 0=all>  \n ");
  if(argc<2)
    exit(0);
  TString tag("run");
  if (argc > 1)
  {
    tag = TString(argv[1]);
  }
  
  countFiles(tag);
  printf(" for %s found %lu files \n", tag.Data(), fileList.size());
  Long64_t maxFiles = fileList.size();
  if (argc > 2)
  {
    maxFiles = atoi(argv[2]);
  }

  anaRun *r = new anaRun(tag);

  cout << " starting ana for  " << maxFiles << " files " << endl;

  for (unsigned i = 0; i < maxFiles; ++i){
    cout << " starting anaRunFile " << fileList[i] << endl;
    r->anaRunFile(fileList[i],0);
  }

  cout << "ana finished " << maxFiles << endl;
  exit(0);
}
