/*  read raw data from SIS3316 125MHz MG, Dec 15 2022 */
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
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
enum
{
  NBYTES = 4,
  BYTE = 8,
  WORD = 4 * BYTE
};

std::uint32_t eventLength;
ifstream input;
std::vector<std::uint32_t> btest;
TBRawRun *rawRun;
TBRawEvent *rawEvent;

std::uint32_t
getLeft(uint32_t a)
{
  return a >> 16;
}

std::uint32_t getRight(uint32_t a)
{
  return (a << 16) >> 16;
}
std::uint32_t getRight25(uint32_t a)
{
  return (a << 6) >> 6;
}
uint32_t readChannel()
{
  std::uint32_t line;
  //cout << "before...at byte " << input.tellg() << endl;
  input.read(reinterpret_cast<char *>(&line), sizeof line);
  uint32_t channel = getRight(line);
  channel = channel >> 4;
  input.seekg(std::streampos(-(sizeof line)), ios_base::cur);
  //cout << "after...at byte " << input.tellg() << endl;
  return channel;
}
  void readEventHeader()
  {
    std::uint32_t line;
    input.read(reinterpret_cast<char *>(&line), sizeof line);
    // std::cout << std::hex << std::showbase << line << dec << '\n';
    // std::bitset<WORD> bits(value);
    // std::cout << bits.to_string() << std::endl;

    if (int(input.tellg()) == 36)
    {
      uint32_t bits[4];
      for (int ib = 0; ib < 4; ++ib)
      {
        bits[ib] = (btest[ib] & line) >> ib;
        printf(" bit %i = %i ", ib, bits[ib]);
      }
      cout << endl;
    }
    uint64_t time2 = getLeft(line);
    uint32_t channel = getRight(line);
    channel = channel >> 4;

    // next line
    // cout << " at byte " << input.tellg() << endl;

    input.read(reinterpret_cast<char *>(&line), sizeof line);
    uint64_t time1 = getLeft(line);
    uint64_t time0 = getRight(line);
    uint64_t time = time0 + (time1 << 16) + (time2 << 32);
    
    //  skip the 9 lines
    input.seekg(std::streampos(NBYTES * 7), ios_base::cur);
    // read energy
    input.read(reinterpret_cast<char *>(&line), sizeof line);
    uint32_t startEnergy = line;
    input.read(reinterpret_cast<char *>(&line), sizeof line);
    uint32_t maxEnergy = line;
    //printf("startEnergy %i maxEnergy %i \n", startEnergy, maxEnergy);
    // read MAW
    input.read(reinterpret_cast<char *>(&line), sizeof line);
    // cout << " at byte " << input.tellg() <<" "<<  hex << line << dec <<  endl;
    
    uint32_t samples = getRight25(line);
    uint32_t mawflag = btest[27] & line;
    //printf("\t header channel %u samples %u \n ", channel, samples);
    rawEvent->channel = channel;
    rawEvent->samples = samples;
    rawEvent->startEnergy = startEnergy;
    rawEvent->maxEnergy = maxEnergy;
    rawEvent->time = time;
    return;
  }

  int main(int argc, char *argv[])
  {
    // fill bitmask
    for (unsigned ib = 0; ib < WORD; ++ib)
    {
      btest.push_back(pow(2, ib));
    }
    cout << "executing " << argv[0] << " argc " << argc << endl;
    if (argc < 2)
    {
      printf(" usage: convertRaw <tag>  <first event> <number  of events> \n ");
      exit(0);
    }
    int first=0;
    int maxEvents=0;
    TString tag(argv[1]);
    if(argc>2)
      first = atoi(argv[2]);
    if (argc > 3)
      maxEvents = atoi(argv[3]);
    printf(" starting convertRaw %s first %i maxEvents %i \n ", tag.Data(), first, maxEvents);
    TString inFileName = TString("data/") + tag + TString(".dat");
    cout << "Open " << inFileName << endl;
    input.open(inFileName.Data(), ios::binary);
    if (input.is_open() == false)
    {
      cout << "Failed! Quit" << endl;
      return -1;
    }
    cout << "Check input file size... ";
    input.seekg(0, ios::end);     // move getter to the end of file
    int fileSize = input.tellg(); // get input file size
    input.seekg(0, ios::beg);     // move getter back to the beginning
    cout << fileSize << " bytes ," << fileSize / 1024 << " kbytes  , " << fileSize / 1024 / 1024 << " Mbytes \n"
         << endl;
    if (fileSize < 400)
    {
      cout << "Too small! Quit" << endl;
      return -1;
    }

    TFile *fout = new TFile(Form("data/rootData/%s.root", tag.Data()), "recreate");
    rawRun = new TBRawRun(tag);
    rawRun->btree->SetAutoFlush(100);
    printf(" new rawRun %s auto %llu \n", rawRun->GetName(), rawRun->btree->GetAutoFlush());

    cout << "Get input file header" << endl;
    unsigned eventLength;
    vector<TString> headerNames;
    headerNames.push_back("Header Marker       ");
    headerNames.push_back("SIS identifier      ");
    headerNames.push_back("bank buffer counter ");
    headerNames.push_back("module channel id   ");
    headerNames.push_back("number of events    ");
    headerNames.push_back("event length        ");
    headerNames.push_back("encoded short event ");
    headerNames.push_back("multi event length  ");
    input.seekg(0, ios::beg);
    for (int iline = 0; iline < 8; ++iline)
    {
      std::uint32_t line;
      input.read(reinterpret_cast<char *>(&line), sizeof line);
      cout << "# " << iline << "  " << headerNames[iline] << " : 0x" << hex << line << dec << " dec " << line << endl;
      if (iline == 5)
        eventLength = line;
    }
    cout << " event length = " << eventLength << " =  " << eventLength * 2 << " samples " << endl;
    cout << " filesize " << fileSize << " bytes " << endl;

    /* read event headers and data */
    uint32_t iline = 0;
    uint32_t eventCount = 0;
    uint32_t eventsRead = 0;
    bool ifBreak = false;
    while (input.good() && input.tellg() < fileSize)
    {
      uint32_t channel = readChannel();
      if (channel == 0)
        ++eventCount;
      if (eventCount-1<first) {
        input.seekg(std::streampos(2128), ios_base::cur);
        //cout << "... skip  "<< channel << " event "<< eventCount << " at byte " << input.tellg() << endl;
        continue;
      }
      if (channel == 0 && rawRun->detList.size() > 1)
      {
        rawRun->fill();
        eventsRead = rawRun->btree->GetEntries();
        cout << " fill rawRun " << eventCount << " " << rawRun->btree->GetEntries() << endl;
      }
      if(maxEvents>0 && eventsRead >= maxEvents){
        cout << "... break at event  " << eventCount << " at byte " << input.tellg() << endl;
        ifBreak = true;
        break;
      }
      if(channel ==0) {
        cout << "... start of event  " << eventCount << " at byte " << input.tellg() << endl;
      }
      // start of channel reads.  Fill tree if not the first time at channel 0 
      
      // fill TBRawEvent header
      rawEvent = rawRun->getDet(channel);
      rawEvent->clear();
      readEventHeader();
      //printf("\t channel %i samples %i bytes %i \n", rawEvent->channel, rawEvent->samples, NBYTES * rawEvent->samples);
      // read data buff
      std::uint32_t line;
      for (uint32_t idigi = 0; idigi < rawEvent->samples; ++idigi)
      {
        input.read(reinterpret_cast<char *>(&line), sizeof line);
        ushort digi0 = getRight(line);
        ushort digi1 = getLeft(line);
        rawEvent->rdigi.push_back(digi0);
        rawEvent->rdigi.push_back(digi1);
        //printf(" idigi %i %u %u %lu \n", idigi, digi0, digi1, rawEvent->rdigi.size());
      }
      // for MAW buff
      input.seekg(std::streampos(WORD), ios_base::cur);
      //rawEvent->print();
      //cout << "... end of event  "<< eventCount << " at byte " << input.tellg() << endl;
      //cout << " _____________ " << endl << endl;
    }
    // fill last event
    if(!ifBreak) rawRun->fill();
    cout << " fill rawRun " << eventCount << " " << rawRun->btree->GetEntries() << endl;
    printf(" end of convertRaw at count %i  ",eventCount);
    printf(" nchannels %lu events written  %llu \n", rawRun->detList.size(), rawRun->btree->GetEntries());
    input.close();
    fout->ls();
    fout->Write();
    fout->Close();

    exit(0);
  }
