/***************************************************************************/
/*                                                                         */
/*  Filename: sis3316_offline.cpp                                          */
/*                                                                         */
/*  Funktion:                                                              */
/*                                                                         */
/*  Autor:                TH                                               */
/*  date:                 27.06.2018                                       */
/*  last modification:    11.11.2021                                       */
/*                                                                         */
/* ----------------------------------------------------------------------- */
/*                                                                         */
/*  SIS  Struck Innovative Systeme GmbH                                    */
/*                                                                         */
/*  Harksheider Str. 102A                                                  */
/*  22399 Hamburg                                                          */
/*                                                                         */
/*  Tel. +49 (0)40 60 87 305 0                                             */
/*  Fax  +49 (0)40 60 87 305 20                                            */
/*                                                                         */
/*  http://www.struck.de                                                   */
/*                                                                         */
/*  � 2021                                                                 */
/*                                                                         */
/***************************************************************************/

// #define SIS3316_DATA_FORMAT_2021
// #define CERN_ROOT_PLOT
// outputs one file for each read file
#include <iostream>
#include <sys/stat.h>
#include <time.h>
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
#include <TH1I.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TFile.h"
// in bobj
#include "TBFile.hxx"
#include "TBRawEvent.hxx"
using namespace std;

/*===========================================================================*/
/* Globals					  			     */
/*===========================================================================*/
#define MAX_NUMBER_LWORDS_64MBYTE 0x1000000 /* 64MByte */
std::uint64_t maxFiles;
TString dirName;
TList *files;
std::map<int, std::string> vfile;
TFile *fout;
TBFile *tbfile;
/* globals	                               		  			     */
int NCHAN = 13;
int NSAMPLES = 1024;
TString tag;
vector<TH1D *> hChan;
vector<TBRawEvent *> rawEvent;
vector<TTree *> ftrees;
std::vector<unsigned short> wave;
FILE *filePointer;
std::uint64_t totalEvents;
bool verbose = false;

unsigned int gl_ch_data[MAX_NUMBER_LWORDS_64MBYTE];

/* methods */
std::uint32_t getLeft(uint32_t a)
{
	return a >> 16;
}

std::uint32_t getRight(uint32_t a)
{
	return (a << 16) >> 16;
}

// count subruns and channels
void countFiles()
{
  TIter next(files);
	TSystemFile *file;
	while ((file = (TSystemFile *)next()))
	{
		string name = string(file->GetName());
		string exten = name.substr(name.find_last_of(".") + 1);
    if (exten != string("dat"))
      continue;
    int s = name.find_last_of(".") - name.find_last_of("_")-1;
    string num = name.substr(name.find_last_of("_")+1,s);
    int inum = atoi(num.c_str());
    //cout << " s= " << s << " " << name << " " << num << endl;
    vfile.insert(std::pair<int, std::string>(inum, name));
  }
}



/* Prototypes */
std::uint64_t processFile(TString fullName);
void getWave(unsigned int *ch_rawdata_buffer, std::vector<unsigned short> &data);

int ReadBufferHeaderCounterNofChannelToDataFile(unsigned int *header_marker, unsigned int *sis3316_indentifier, unsigned int *bank_loop_no, unsigned int *channel_no, unsigned int *nof_events, unsigned int *event_length, unsigned int *maw_length, unsigned int *reserved_ch_bankx_buffer_length);
int ReadEventsFromDataFile(unsigned int *memory_data_array, unsigned int nof_write_length_lwords);

/* ***************************************************************************************************************** */
int main(int argc, char *argv[])
{
	cout << "executing " << argv[0] << " argc " << argc << endl;
	char dirTag[512];
	int int_ch;
	if(argc==1) {
		printf("Usage: %s  [-?h] [-D <dir> maxFiles (default -1=all) \n", argv[0]);
		exit(0);
	}

	if (argc > 1) 
	{
		while ((int_ch = getopt(argc, argv, "?hD:")) != -1)
		{
			switch (int_ch)
			{
			printf("ch %c   %s  \n", int_ch, optarg );

			case 'D':
				sscanf(optarg, "%s", dirTag);
				break;

			case '?':
			case 'h':
			default:
				// printf("Usage: %s  [-?h] [-I ip] [-A num]  ", argv[0]);
				printf("Usage: %s  [-?h] [-D <dir> maxFiles (default -1=all) ", argv[0]);
				printf("   \n");
				printf("   \n");
				printf("   -D <dir>  read from directory  data/<dir> \n");
				printf("   \n");
				printf("   -h     ..........  print this message only\n");
				printf("   \n");
				exit(1);
			}
		}
	}
	maxFiles= 0;
	if (argc > 3)
		maxFiles = std::uint64_t(atoi(argv[3]));

	printf(" running %s tag %s maxFiles = %llu \n", argv[0], dirTag, maxFiles);

	tag = TString(dirTag);

	dirName = TString("data/") + tag;
  TSystemDirectory dir("rawdir", dirName); // TSystemDirectory
	files = dir.GetListOfFiles();			 // TList

  countFiles();
	cout << "dirName " << dirName << " has "  << vfile.size() <<" .dat files:" << endl;

  for (std::map<int,string>::iterator it = vfile.begin(); it != vfile.end(); ++it)
    std::cout << it->first << " => " << it->second << '\n';
	

	/* loop over files  */
  totalEvents = 0;
  unsigned filesRead=0;
  unsigned totalFiles = vfile.size();
  if(maxFiles>0)
    totalFiles = maxFiles;
  for (std::map<int, string>::iterator it = vfile.begin(); it != vfile.end(); ++it)
    {
    string fname = it->second;
    TString fullName = dirName + TString("/") + TString(fname.c_str());
    printf("\n\n\t starting file %s \n", fullName.Data());
    // open output file
    fout = new TFile(Form("rootData/run-%s-file%u.root",dirTag,filesRead), "recreate");
    printf("opened output file %s \n", fout->GetName());
    // output trees for each channel
    tbfile = new TBFile(fullName);
    fout->Append(tbfile);
    rawEvent.clear();
    ftrees.clear();
    hChan.clear();
    rawEvent.resize(NCHAN);
    ftrees.resize(NCHAN);
    for (int ichan = 0; ichan < NCHAN; ++ichan)
    {
      rawEvent[ichan] = new TBRawEvent(ichan);
      ftrees[ichan] = new TTree(Form("tchan%i", ichan), Form("channel %i", ichan));
      ftrees[ichan]->Branch(rawEvent[ichan]->GetName(), &rawEvent[ichan]);
      printf(" new tree  %s \n", ftrees[ichan]->GetName());
    }

    for (unsigned ich = 0; ich < NCHAN; ++ich)
    {
      TString hname;
      TString htitle;
      hname.Form("channel%i", ich);
      htitle.Form("channel %i", ich);
      hChan.push_back(new TH1D(hname, htitle, NSAMPLES, 1, NSAMPLES));
    }

 
    totalEvents += processFile(fullName);
    printf("\t FINISHED  file %s totalEvents %lld \n", fullName.Data(), totalEvents);
    printf("Entries by channel \n");
    for (int ichan = 0; ichan < NCHAN; ++ichan)
         printf("%s %llu \n", ftrees[ichan]->GetName(), ftrees[ichan]->GetEntries());
    fout->Write();
    fout->Close();
    ++filesRead;
    if(filesRead>=totalFiles)
      break;
    }
	printf(" FINISHED sis after file %i total events %llu  \n",filesRead,totalEvents);

	return 0;
}

//---------------------------------------------------------------------------

int ReadBufferHeaderCounterNofChannelToDataFile(unsigned int *header_marker, unsigned int *sis3316_indentifier, unsigned int *bank_loop_no, unsigned int *channel_no, unsigned int *nof_events, unsigned int *event_length, unsigned int *maw_length, unsigned int *reserved_ch_bankx_buffer_length)
{
	int nof_read = 0;
	// header
	nof_read = fread(header_marker, 0x4, 0x1, filePointer);
	nof_read = nof_read + fread(sis3316_indentifier, 0x4, 0x1, filePointer);
	nof_read = nof_read + fread(bank_loop_no, 0x4, 0x1, filePointer);
	nof_read = nof_read + fread(channel_no, 0x4, 0x1, filePointer);
	nof_read = nof_read + fread(nof_events, 0x4, 0x1, filePointer);
	nof_read = nof_read + fread(event_length, 0x4, 0x1, filePointer);
	nof_read = nof_read + fread(maw_length, 0x4, 0x1, filePointer);
	nof_read = nof_read + fread(reserved_ch_bankx_buffer_length, 0x4, 0x1, filePointer);

	return nof_read;
}

//---------------------------------------------------------------------------
int ReadEventsFromDataFile(unsigned int *memory_data_array, unsigned int nof_write_length_lwords)
{
	int nof_read;

	nof_read = fread(memory_data_array, 0x4, nof_write_length_lwords, filePointer);
	return nof_read;
}

void getWave(unsigned int *ch_rawdata_buffer, std::vector<unsigned short> &wave)
{
	unsigned short *ushort_adc_buffer_ptr;
	ushort_adc_buffer_ptr = (unsigned short *)ch_rawdata_buffer;
	// module_ch_no = ch_no & 0xf;
	for (unsigned i = 0; i < wave.size(); i++)
	{
		wave[i] = (ushort_adc_buffer_ptr[i]);
	}
}

std::uint64_t processFile(TString fullName)
{

  std::uint64_t fileEvents=0;
  //
  unsigned valid_BankBufferHeader_valid_flag;
  int nof_read;
  unsigned buffer_no;
  unsigned i_event;
  unsigned buffer_length;
  unsigned event_length;
  unsigned header_length;
  unsigned sample_length;
  unsigned channel_no;
  unsigned short_event_and_maw_buffer_length;
  unsigned i_ch;
  unsigned headerformat;
  unsigned header_marker;
  unsigned header_indentifier;
  unsigned header_reserved_ch_bankx_buffer_length;
  unsigned int uint_save_raw_sample_first_event_only_mode;
  unsigned int maw_buffer_length;
  unsigned uint_current_bankx_buffer_counter;
  unsigned nof_events;

  filePointer = fopen(fullName.Data(), "rb");
  if (filePointer == NULL)
  {
    printf("file %s filePointer == NULL \n", fullName.Data());
    return fileEvents;
  }

  if (filePointer != NULL){

    uint_current_bankx_buffer_counter = 0; //

    do // valid_BankBufferHeader_valid_flag = 1;
    {
      valid_BankBufferHeader_valid_flag = 0;
      nof_read = ReadBufferHeaderCounterNofChannelToDataFile(&header_marker, &header_indentifier, &buffer_no, &channel_no, &nof_events, &event_length, &short_event_and_maw_buffer_length, &header_reserved_ch_bankx_buffer_length);

      if (nof_read != 8)
      {
        valid_BankBufferHeader_valid_flag = 0;
        printf("Stop Condition: file %s no Header data (EOF) nof_read = %i \n", fullName.Data(), int(nof_read));
        return fileEvents;
      }
      else
      { // valid header length
        if(verbose) printf(" buffer_no %i header: marker = 0x%08x    identifier = 0x%08x  header length = %d  totalEvents %llu  \n", buffer_no, header_marker, header_indentifier, nof_read, totalEvents);

        if (header_marker != 0xDEADBEEF)
        {
          valid_BankBufferHeader_valid_flag = 0;
          printf("Stop Condition: In file %s no valid header-marker (0xDEADBEEF) after %llu file events \n", fullName.Data(),fileEvents);
        }
        else
        {
          valid_BankBufferHeader_valid_flag = 1;
          if (uint_current_bankx_buffer_counter != buffer_no) uint_current_bankx_buffer_counter = buffer_no;

          // print Buffer Header
          /*
             printf("Chx Bankx buffer header information: %i \n", nof_events);
             printf("\tindentifier       = %d   \tbuffer_no    = %d    \tchannel_id         = %d   \n", header_indentifier, buffer_no, channel_no);
             printf("\tnof_events        = %d   \tevent_length = %d    \tshort_event_length = %d   \traw_data_first_event_only_flag = %ld  \n", nof_events, event_length, (short_event_and_maw_buffer_length & 0x7fff) >> 16, (short_event_and_maw_buffer_length & 0x780000000) >> 31);
             printf("\tmaw_buffer_length = %d    \n", (short_event_and_maw_buffer_length & 0x7fff));
             printf("\treserved_ch_bankx_buffer_length = %d (0x%08x) \n", header_reserved_ch_bankx_buffer_length, header_reserved_ch_bankx_buffer_length);
             */

          maw_buffer_length = (short_event_and_maw_buffer_length & 0x7fff);

          // read event by event from file (has to check "sample length" of each event)
          nof_read = 0;
          for (i_event = 0; i_event < nof_events; i_event++)
          {
            nof_read = nof_read + ReadEventsFromDataFile(&gl_ch_data[0], 1); // read 4 bytes = 32 bits
            i_ch = (gl_ch_data[0] & 0xfff0) >> 4;
            headerformat = (gl_ch_data[0] & 0xf);

            header_length = 3; // if headerformat == 0
            if ((headerformat & 0x1) == 1)
            {
              header_length = header_length + 7;
            }
            if ((headerformat & 0x2) == 2)
            {
              header_length = header_length + 2;
            }
            if ((headerformat & 0x4) == 4)
            {
              header_length = header_length + 3;
            }
            if ((headerformat & 0x8) == 8)
            {
              header_length = header_length + 2;
            }

            // read rest of Event-Header
            nof_read = nof_read + ReadEventsFromDataFile(&gl_ch_data[1], header_length - 1);
            /* make the time  */
            uint64_t time2 = getLeft(gl_ch_data[0]);
            uint64_t time1 = getLeft(gl_ch_data[1]);
            uint64_t time0 = getRight(gl_ch_data[1]);
            uint64_t time = time0 + (time1 << 16) + (time2 << 32);

            sample_length = 2 * (gl_ch_data[header_length - 1] & 0x3ffffff); // bits 0 - 25
            if (sample_length != 0)
            {
              nof_read = nof_read + ReadEventsFromDataFile(&gl_ch_data[header_length], sample_length / 2); // read Raw Data
              wave.clear();
              wave.resize(sample_length);
              getWave(&gl_ch_data[header_length], wave);
              rawEvent[i_ch]->channel = i_ch;
              rawEvent[i_ch]->buffer = buffer_no;
              rawEvent[i_ch]->trigger = i_event;
              rawEvent[i_ch]->time = time;
              rawEvent[i_ch]->rdigi = wave;
              ftrees[i_ch]->Fill();

              // simple baseline
              double base = 0;
              unsigned nsum = 50;
              for (unsigned ib = 0; ib < nsum; ++ib)
                base += double(wave[ib]);
              base /= double(nsum);
              // hChan[i_ch]->Reset("ICESM");
              for (unsigned ib = 0; ib < wave.size(); ++ib)
                hChan[i_ch]->SetBinContent(ib + 1, double(wave[ib]) - base + hChan[i_ch]->GetBinContent(ib + 1));

              // if (i_event == nof_events - 1)
              //	printf("\t ...chan %i i_event %i length %lu  baseline %f (%f) \n", i_ch, i_event + 1, wave.size(), base, double(wave[0]));
            }
            if (maw_buffer_length != 0)
            {
              nof_read = nof_read + ReadEventsFromDataFile(&gl_ch_data[header_length + (sample_length / 2)], maw_buffer_length); // read MAW Data
            }
          }
          fileEvents += nof_events;
          if(verbose) printf("...finished fileEvents %llu nof_events = %i  buffer %u  \tChannel ID = %d   \tFormat bits = 0x%02X \n", fileEvents, nof_events, buffer_no, i_ch, headerformat);
        }
      }

    } while (valid_BankBufferHeader_valid_flag == 1);

    fclose(filePointer);
  }
  printf(" return from processEvents file %s fileEvents %lld\n",fullName.Data(),fileEvents );
  return fileEvents;
} 
