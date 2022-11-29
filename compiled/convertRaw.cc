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
#include "TBRawRun.hxx"

using namespace std;

TString dirName;
TList *files;
std::vector<int> vsub;
std::vector<int> vchan;
TBRawRun *rawRun;
unsigned maxRun;
unsigned totalEvents;
/*
Record Length: 4100
BoardID: 31
Channel: 6
Event Number: 0
Pattern: 0x0000
Trigger Time Stamp: 213649
DC offset (DAC): 0x1999
*/
int readHeader(int chan, ifstream *stream)
{
  if (stream->peek() != 'R')
    return 0;
  int nlines = 0;
  int place[7] = {2, 1, 1, 2, 1, 3, 3};
  char line[256];
  TBRawEvent *rev = rawRun->getDet(chan);
  if (!rev)
  {
    printf(" no TBRawEvent for chan %i !!!\n", chan);
    return nlines;
  }
  rev->clear();
  for (int i = 0; i < 7; ++i)
  {
    stream->getline(line, 256);
    if (!stream->good())
      return nlines;
    TString tline(line);
    TObjArray *tokenArray = tline.Tokenize(' ');
    TObjString *sv = (TObjString *)tokenArray->At(place[i]);
    // printf(" header line %i %s \n .... place %i = %s \n", nlines, line , place[i], sv->GetString().Data());
    if (nlines == 0)
      rev->length = sv->GetString().Atoi();
    if (nlines == 1)
      rev->boardID = sv->GetString().Atoi();
    // if(nlines==2) rev->channel = sv->GetString().Atoi();
    if (nlines == 3)
      rev->event = sv->GetString().Atoi();
    if (nlines == 5)
      rev->time = sv->GetString().Atoll();
    if (nlines == 6)
      rev->dcOffset = sv->GetString().Atoi();
    ++nlines;
  }
  //rev->printHeader();
  return nlines;
}

void parseName(TString name, int &irun, int &ichan)
{
  TObjArray *tokenArray = name.Tokenize('_');
  // int narray = tokenArray->GetEntries();
  //  for (int i = 0; i < narray; ++i)
  //    printf("token %i %s \n", i, ((TObjString *)tokenArray->At(i))->GetString().Data());
  TObjString *srun = (TObjString *)tokenArray->At(1);
  irun = srun->GetString().Atoi();
  TObjString *schan = (TObjString *)tokenArray->At(3);
  ichan = schan->GetString().Atoi();
  // printf("%s run %i chan %i \n", name.Data(), irun, ichan);
}

// count subruns and channels
void count(TString tag)
{
  TIter next(files);
  TSystemFile *file;
  while ((file = (TSystemFile *)next()))
  {
    string name = string(file->GetName());
    string exten = name.substr(name.find_last_of(".") + 1);
    if (exten != string("dat"))
      continue;
    string tag = name.substr(0, name.find_last_of("."));
    string fullName = string(dirName.Data()) + string("/") + name;
    int subRun, chan;
    TString tname(name.c_str());
    parseName(tname, subRun, chan);
    // printf("\t\t %s %i %i \n",name.c_str(),subRun,chan);
    if (std::find(std::begin(vsub), std::end(vsub), subRun) == std::end(vsub))
      vsub.push_back(subRun);
    if (std::find(std::begin(vchan), std::end(vchan), chan) == std::end(vchan))
      vchan.push_back(chan);
  }
}


bool readChan(int chan, ifstream *stream)
{
  char line[256];
  int ievent = 0;

  while (1)
  {
    // read header for event
    if (stream->peek() == 'R')
    {
      ++ievent;
      if (ievent > 1) return stream->good();
      readHeader(chan, stream);
    }
    else
    {
      stream->getline(line, 256);
      if (!stream->good())
        break;

      TString tline(line);
      UShort_t val = UShort_t(tline.Atoi());
      rawRun->getDet(chan)->rdigi.push_back(val);
    }
    //printf(" chan %i stream pos %i \n", chan,pos);
  }
  return stream->good();
}

void openSubRun(int isub)
{
  TIter next(files);
  TSystemFile *file;
  std::vector<int> vchan;
  std::vector<ifstream *> fileList;
  std::vector<string> fileName;

  // open all files for this subrun
  while ((file = (TSystemFile *)next()))
  {
    string name = string(file->GetName());
    string exten = name.substr(name.find_last_of(".") + 1);
    if (exten != string("dat"))
      continue;
    string tag = name.substr(0, name.find_last_of("."));
    string fullName = string(dirName.Data()) + string("/") + name;
    int subRun, chan;
    TString tname(name.c_str());
    parseName(tname, subRun, chan);
    // printf(" %s subrun %i chan %i \n",tname.Data(),subRun,chan);
    if (subRun != isub) continue;

    ifstream *in = new ifstream(fullName, std::ios::in);
    if (!in->is_open()) continue;

    in->seekg(0, std::ios::end); // go to the end
    printf(" subrun %i opened file %s length %i \n", isub, tname.Data(), int(in->tellg()));
    in->seekg(0, std::ios::beg);

    // count channels
    if (std::find(std::begin(vchan), std::end(vchan), chan) == std::end(vchan))
    {
      vchan.push_back(chan);
      fileList.push_back(in);
      fileName.push_back(name);
    }
  }

  printf(" subrun %i number of files opened %lu \n", isub, fileList.size());
  //for (unsigned i = 0; i < fileName.size(); ++i)
  //  printf("\t %i %s \n", i, fileName[i].c_str());

  bool good = true;
  int nevents = 0;
  int nbytes = 0;
  // loop over events in subrun
  while (good)
  {
    ++nevents;
    /* read all channels for this event */
    for (unsigned i = 0; i < fileList.size(); ++i)
    {
      good = readChan(vchan[i], fileList[i]);
      if (!good)
        break;
      // end of reading subrun, close all files for this subrun
    }
    nbytes += rawRun->fill();
    rawRun->clear();
    if (nevents / 1000 * 1000 == nevents)
      printf(" events read subrun %i event %i entries %llu bytes %i \n", isub, nevents, rawRun->btree->GetEntries(),nbytes);
  }

  // end of reading subrun, close all files for this subrun
  for (unsigned i = 0; i < fileList.size(); ++i)
    fileList[i]->close();

  totalEvents += nevents;
  printf(" end of reading subrun %i total events  %i entries %llu \n", isub,totalEvents ,rawRun->btree->GetEntries());
}

int main(int argc, char *argv[])
{
  totalEvents = 0;
  cout << "executing " << argv[0] << " argc " << argc << endl;
  if (argc < 2)
  {
    printf(" usage: convertRaw <tag>  <maxruns> \n ");
    exit(0);
  }
  cout << argv[1] << endl;
  TString tag(argv[1]);

  maxRun = 0;
  if (argc > 2)
  {
    cout << "max runs = " << argv[2] << endl;
    maxRun = (unsigned)atoi(argv[2]);
  }

  dirName = TString("data/") + tag;
  TSystemDirectory dir("rawdir", dirName); // TSystemDirectory
  files = dir.GetListOfFiles();            // TList
  cout << "dirName " << dirName << " has " << files->GetEntries() << " files " << endl;

  count(tag);
  if (vsub.size() < 1)
    return 1;
  printf(" *** found %lu subruns %lu chan *** \n", vsub.size(), vchan.size());
  for (unsigned i = 0; i < vchan.size(); ++i)
    printf("\t chan %i \n", vchan[i]);

  TFile *fout = new TFile(Form("data/rootData/%s.root", tag.Data()), "recreate");
  rawRun = new TBRawRun(tag);
  //rawRun->btree->SetAutoSave(-300000000);
  /* autof is the number of consecutive entries after which TTree::Fill will flush all branch buffers to disk.*/
  rawRun->btree->SetAutoFlush(100);
  printf(" new rawRun %s auto %llu \n", rawRun->GetName(), rawRun->btree->GetAutoFlush());

  for (unsigned i = 0; i < vchan.size(); ++i)
    rawRun->addDet(vchan[i]);

  std::sort(vsub.begin(), vsub.end());

  if (maxRun == 0)
    maxRun = vsub.size();

  for (unsigned isub = 0; isub < maxRun; ++isub)
  {
    printf("\t *****  starting subrun %u  number %u \n ****** ", isub, vsub[isub]);
    openSubRun(vsub[isub]);
  }

  fout->ls();
  fout->Write();
  fout->Close();

  exit(0);
}
