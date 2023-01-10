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

#include <iostream>
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
using namespace std;

#ifdef CERN_ROOT_PLOT
sis_root_graph *gl_graph_raw;
sis_root_graph_maw *gl_graph_maw;
#endif

/*===========================================================================*/
/* Globals					  			     */
/*===========================================================================*/
#define MAX_NUMBER_LWORDS_64MBYTE 0x1000000 /* 64MByte */

unsigned int gl_ch_data[MAX_NUMBER_LWORDS_64MBYTE];

FILE *gl_FILE_DataEvenFilePointer;
TFile *fout;
TTree *ftree;

/*===========================================================================*/
/* Prototypes			                               		  			     */
/*===========================================================================*/
int NCHAN = 13;
int NSAMPLES = 1024;
vector<TH1I *> hChan;
void getWave(unsigned int *ch_rawdata_buffer, std::vector<unsigned short>& data);

int ReadBufferHeaderCounterNofChannelToDataFile(unsigned int *header_marker, unsigned int *sis3316_indentifier, unsigned int *bank_loop_no, unsigned int *channel_no, unsigned int *nof_events, unsigned int *event_length, unsigned int *maw_length, unsigned int *reserved_ch_bankx_buffer_length);
int ReadEventsFromDataFile(unsigned int *memory_data_array, unsigned int nof_write_length_lwords);

/* ***************************************************************************************************************** */
int main(int argc, char *argv[])
{
	cout << "sis3316_offline" << endl; // prints sis3316_offline

	unsigned valid_BankBufferHeader_valid_flag;
	int nof_read;
	unsigned buffer_no;
	unsigned i_event;
	unsigned nof_events;
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

	// unsigned int i ;
	// unsigned int adc_offset ;
	// unsigned int adc_length ;
	// unsigned int maw_offset ;
	// unsigned int maw_length ;
	unsigned int maw_buffer_length;

	unsigned uint_plot_axis_flag;
	unsigned uint_current_bankx_buffer_counter;

	unsigned bank_buffer_counter;

	char filetag[512];
	// unsigned int i_file;
	int int_ch;

	if (argc > 1)
	{

		while ((int_ch = getopt(argc, argv, "?hF:")) != -1)
		{
			switch (int_ch)
			{
				// printf("ch %c    \n", int_ch );

			case 'F':
				sscanf(optarg, "%s", filetag);
				break;

			case '?':
			case 'h':
			default:
				// printf("Usage: %s  [-?h] [-I ip] [-A num]  ", argv[0]);
				printf("Usage: %s  [-?h] [-F data_filename.dat] ", argv[0]);
				printf("   \n");
				printf("   \n");
				printf("   -F file tag.  read from file  \n");
				printf("   \n");
				printf("   -h     ..........  print this message only\n");
				printf("   \n");
				printf("   \n");
				printf("   date:               11. November 2021 \n");
				printf("   \n");
				printf("   \n");
				exit(1);
			}
		}
	} // if (argc > 2)

	printf("\n");


	bank_buffer_counter = 0;

	// sprintf(filename,"../data_files/sample_test_gui/sis3316_test_data_%d.dat",i_file ) ;
	char filename[120];
	sprintf(filename, "data/%s.dat", filetag);
	printf("tag %s filename %s\n", filetag, filename);
	gl_FILE_DataEvenFilePointer = fopen(filename, "rb");
	if (gl_FILE_DataEvenFilePointer == NULL)
	{
		printf("gl_FILE_DataEvenFilePointer == NULL \n");
		return -1;
	}

	// open output file
	fout  = new TFile(Form("%s.root",filetag),"recreate");

	for (unsigned ich = 0; ich < NCHAN;++ ich){
		TString hname;
		TString htitle;
		hname.Form("channel%i",ich);
		htitle.Form("channel %i",ich);
		hChan.push_back(new TH1I(hname, htitle, NSAMPLES, 1, NSAMPLES));
	}

	ftree = new TTree("raw","raw data");
	int bch;
	ftree->Branch("channel", &bch);
	int bnbuf;
	ftree->Branch("buff", &bnbuf);
	int bevent;
	ftree->Branch("event", &bevent);
	std::vector<unsigned short> wave;
	auto branch = ftree->Branch("data", &wave, NSAMPLES);


	if (gl_FILE_DataEvenFilePointer != NULL)
	{
		uint_plot_axis_flag = 1;			   // set first time:  clear graph and draw x-y axis flag
		uint_current_bankx_buffer_counter = 0; //

		do
		{
			valid_BankBufferHeader_valid_flag = 0;
			nof_read = ReadBufferHeaderCounterNofChannelToDataFile(&header_marker, &header_indentifier, &buffer_no, &channel_no, &nof_events, &event_length, &short_event_and_maw_buffer_length, &header_reserved_ch_bankx_buffer_length);
			if (buffer_no > 0)
				goto CLOSE;
			// printf("\n");
			printf(" bank_buffer_counter %i header: marker = 0x%08x    identifier = 0x%08x      length = %d    \n", bank_buffer_counter, header_marker, header_indentifier, nof_read);
			++bank_buffer_counter;

			if (nof_read != 8)
			{
				valid_BankBufferHeader_valid_flag = 0;
				printf("Stop Condition: no Header data (EOF) nof_read = %i \n", int(nof_read));
				goto CLOSE;
			}
			else
			{ // valid header length
				if (header_marker != 0xDEADBEEF)
				{
					valid_BankBufferHeader_valid_flag = 0;
					printf("Stop Condition: no valid header-marker (0xDEADBEEF) \n");
				}
				else
				{
					valid_BankBufferHeader_valid_flag = 1;
					if (uint_current_bankx_buffer_counter != buffer_no)
					{							 // new bankx buffer no/counter
						uint_plot_axis_flag = 1; // set clear graph and draw x-y axis flag
						uint_current_bankx_buffer_counter = buffer_no;
						// printf("Set uint_plot_axis_flag  \n");
					}

					// print Buffer Header
					printf("Chx Bankx buffer header information: %i \n", nof_events);
					printf("\tindentifier       = %d   \tbuffer_no    = %d    \tchannel_id         = %d   \n", header_indentifier, buffer_no, channel_no);
					printf("\tnof_events        = %d   \tevent_length = %d    \tshort_event_length = %d   \traw_data_first_event_only_flag = %ld  \n", nof_events, event_length, (short_event_and_maw_buffer_length & 0x7fff) >> 16, (short_event_and_maw_buffer_length & 0x780000000) >> 31);
					printf("\tmaw_buffer_length = %d    \n", (short_event_and_maw_buffer_length & 0x7fff));
					printf("\treserved_ch_bankx_buffer_length = %d (0x%08x) \n", header_reserved_ch_bankx_buffer_length, header_reserved_ch_bankx_buffer_length);

					maw_buffer_length = (short_event_and_maw_buffer_length & 0x7fff);

					// read event by event from file (has to check "sample length" of each event)
					nof_read = 0;
					for (i_event = 0; i_event < nof_events; i_event++)
					{
						nof_read = nof_read + ReadEventsFromDataFile(&gl_ch_data[0], 1);
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
						nof_read = nof_read + ReadEventsFromDataFile(&gl_ch_data[1], header_length - 1); // read rest of Event-Header

						sample_length = 2 * (gl_ch_data[header_length - 1] & 0x3ffffff); //
						if (sample_length != 0)
						{
							nof_read = nof_read + ReadEventsFromDataFile(&gl_ch_data[header_length], sample_length / 2); // read Raw Data
							wave.clear();
							wave.resize(sample_length);
							getWave(&gl_ch_data[header_length],wave);
							bch = i_ch;
							bnbuf = buffer_no;
							bevent = i_event;
							ftree->Fill();

							for (unsigned ib = 0; ib < wave.size(); ++ib)
							{
								hChan[i_ch]->SetBinContent(ib + 1, wave[ib]  + hChan[i_ch]->GetBinContent(ib + 1));
								}
							if (i_event == nof_events-1)
								printf("\t ...chan %i i_event %i length %lu  \n",i_ch, i_event+1, wave.size());
						}
						if (maw_buffer_length != 0)
						{
							nof_read = nof_read + ReadEventsFromDataFile(&gl_ch_data[header_length + (sample_length / 2)], maw_buffer_length); // read MAW Data
						}
					}
					printf("...finished nof_events = %d   nof_read = %d  \tChannel ID = %d   \tFormat bits = 0x%02X \n\n", nof_events, nof_read, (gl_ch_data[0] & 0xf0) >> 4, headerformat);
				}
			}

			//gSystem->ProcessEvents(); // handle GUI events
				
		} while (valid_BankBufferHeader_valid_flag == 1);
//fclose(gl_FILE_DataEvenFilePointer);
//printf("file closed and finished   \n");
//bank_buffer_counter++;
} // good fil

//}

CLOSE :

	printf(" close with  bank_buffer_counter %i buffer_no %i nof_events  = %i ftree entries %lld \n", bank_buffer_counter, buffer_no , nof_events,ftree->GetEntries());
	fclose(gl_FILE_DataEvenFilePointer);
	fout->ls();
	fout->Write();
	fout->Close();
	return 0;
}

//---------------------------------------------------------------------------

int ReadBufferHeaderCounterNofChannelToDataFile(unsigned int *header_marker, unsigned int *sis3316_indentifier, unsigned int *bank_loop_no, unsigned int *channel_no, unsigned int *nof_events, unsigned int *event_length, unsigned int *maw_length, unsigned int *reserved_ch_bankx_buffer_length)
{
	int nof_read = 0;
	// header
	nof_read = fread(header_marker, 0x4, 0x1, gl_FILE_DataEvenFilePointer);
	nof_read = nof_read + fread(sis3316_indentifier, 0x4, 0x1, gl_FILE_DataEvenFilePointer);
	nof_read = nof_read + fread(bank_loop_no, 0x4, 0x1, gl_FILE_DataEvenFilePointer);
	nof_read = nof_read + fread(channel_no, 0x4, 0x1, gl_FILE_DataEvenFilePointer);
	nof_read = nof_read + fread(nof_events, 0x4, 0x1, gl_FILE_DataEvenFilePointer);
	nof_read = nof_read + fread(event_length, 0x4, 0x1, gl_FILE_DataEvenFilePointer);
	nof_read = nof_read + fread(maw_length, 0x4, 0x1, gl_FILE_DataEvenFilePointer);
	nof_read = nof_read + fread(reserved_ch_bankx_buffer_length, 0x4, 0x1, gl_FILE_DataEvenFilePointer);

	return nof_read;
}

//---------------------------------------------------------------------------
int ReadEventsFromDataFile(unsigned int *memory_data_array, unsigned int nof_write_length_lwords)
{
	int nof_read;

	nof_read = fread(memory_data_array, 0x4, nof_write_length_lwords, gl_FILE_DataEvenFilePointer);
	return nof_read;
}

void getWave(unsigned int *ch_rawdata_buffer, std::vector<unsigned short> &wave )
{
	unsigned short *ushort_adc_buffer_ptr;
	ushort_adc_buffer_ptr = (unsigned short *)ch_rawdata_buffer;
	//module_ch_no = ch_no & 0xf;
	for (unsigned i = 0; i < wave.size(); i++)
	{
		wave[i]=(ushort_adc_buffer_ptr[i]);
	}
}
