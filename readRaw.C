// Parse data file  from a Struck ADC
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdint>
#include <bitset>
#include <algorithm>
using namespace std;
enum
{
	NBYTES = 4 ,
	BYTE = 8 ,
	WORD = 4 * BYTE
};
std::uint32_t eventLength;
std::uint32_t multiEventLength;

ifstream input;
int format;
int timestamp;
//std::vector<bitset<WORD>> bitMask;
std::vector<std::uint32_t> btest;
std::string
toBinary(std::uint32_t n)
{
	string result;
	do
		result.push_back('0' + (n & 1));
	while (n >>= 1);
	reverse(result.begin(), result.end());
	return result;
}

std::uint32_t getLeft(uint32_t a)
{
	return a >> 16;
}

std::uint32_t getRight(uint32_t a)
{
	return (a << 16)>>16 ;
}
std::uint32_t getRight25(uint32_t a)
{
	return (a << 6) >> 6;
}

uint32_t readEventHeader(int iline, std::uint32_t &channel)
{
	std::uint32_t nlines=0;
	cout << ">>>> header for "<< iline << " at byte " << input.tellg() << endl;
	std::uint32_t line;
	input.read(reinterpret_cast<char *>(&line), sizeof line);
	++nlines;
	// std::cout << std::hex << std::showbase << line << dec << '\n';
	// std::bitset<WORD> bits(value);
	// std::cout << bits.to_string() << std::endl;

	if (int(input.tellg())==36)
	{
		uint32_t bits[4];
		for (int ib = 0; ib < 4; ++ib)
		{
			bits[ib] = (btest[ib] & line)>>ib;
			printf(" bit %i = %i ", ib, bits[ib]);
		}
		cout << endl;
	}
	uint64_t time2 = getLeft(line);
	channel = getRight(line);
	channel = channel>>4;
	// next line
	//cout << " at byte " << input.tellg() << endl;

	input.read(reinterpret_cast<char *>(&line), sizeof line);
	++nlines;
	uint64_t time1 = getLeft(line);
	uint64_t time0 = getRight(line);
	uint64_t time = time0+ (time1<<16)  + (time2 <<32);
	printf("channel %u time2 %llx time1 %llx  time0 %llx time %llx %llu \n ", channel, time2, time1, time0, time,time);

	//cout << " at byte " << input.tellg() << " "<< hex  << line << dec  << endl;
	// skip the 9 lines
	input.seekg(std::streampos(NBYTES*7), ios_base::cur);
	nlines += 7;
	// read energy
	input.read(reinterpret_cast<char *>(&line), sizeof line);
	++nlines;
	uint32_t startEnergy = line;
	input.read(reinterpret_cast<char *>(&line), sizeof line);
	++nlines;
	uint32_t maxEnergy = line;
	printf("startEnergy %i maxEnergy %i \n", startEnergy, maxEnergy);
	// read MAW
	input.read(reinterpret_cast<char *>(&line), sizeof line);
	++nlines; 
	// cout << " at byte " << input.tellg() <<" "<<  hex << line << dec <<  endl;
	uint32_t samples = getRight25(line);
	uint32_t mawflag = btest[27] & line;
	cout << " nlines " << nlines << " number raw samples " << samples << " MAW test flag " << mawflag << endl;
	return samples;
}

void readRaw(const char *inputFileName  = "sis3316_test_data.dat",
			 const char *index_folder = "./")
{
	// fill bitmask
	for (unsigned ib = 0; ib < WORD; ++ib)
	{
		//bitMask.push_back(pow(2, ib));
		btest.push_back(pow(2, ib));
		//cout<< hex  << bitMask[ib].to_ullong() << "   "  << btest[ib] << dec << endl;
	}

  string input_file = string("compiled/data/")+string(inputFileName);
	cout << "Open " << input_file << endl;
	input.open(input_file, ios::binary);
	if (input.is_open() == false)
	{
		cout << "Failed! Quit" << endl;
		return;
	}

	cout << "Check input file size... ";
	input.seekg(0, ios::end);	  // move getter to the end of file
	int fileSize = input.tellg(); // get input file size
	input.seekg(0, ios::beg);	  // move getter back to the beginning
	cout << fileSize << " bytes ," << fileSize / 1024 << " kbytes  , " << fileSize / 1024 / 1024 << " Mbytes \n"
		 << endl;
	if (fileSize < 400)
	{
		cout << "Too small! Quit" << endl;
		return;
	}

	
	cout << "Get input file header" << endl;
	int word[100];			   // read not more than 100 words at a time
	char *byte = (char *)word; // index in unit of byte instead of word
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
		input.read(byte, NBYTES); // get 1st line of file header
		cout << "# " << iline << "  " << headerNames[iline] << " : 0x" << hex << word[0] << dec << " dec " << word[0] << endl;
		if (iline == 5)
			eventLength = word[0];
		if (iline == 7)
			multiEeventLength = word[0];
	}
	cout << " event length = " << eventLength << " =  " << eventLength * 2 << " samples " << endl;
	cout << " filesize "<< fileSize << " bytes " << " multiEventLength " << multiEventLength << << endl;
	uint32_t iline = 0;
	uint32_t nchannels = 0;
	uint32_t nevents = 0;
	while (input.good() && input.tellg() < fileSize)
	{
		uint32_t channel;
		uint32_t samples = readEventHeader(iline++,channel);
		// count channels and events
		if (channel == 0)
			++nevents;
		if(nevents==1)
			++nchannels;

		// cout << "\t at  byte " << input.tellg() << endl;
		printf("\t channel %i samples %i bytes %i \n", channel, samples, NBYTES * samples);
		// data buff 
		input.seekg(std::streampos(NBYTES * samples), ios_base::cur);
		// for MAW buff
		//cout << "before...at byte " << input.tellg() << endl;
		input.seekg(std::streampos(WORD), ios_base::cur);
		//cout << "after...at byte " << input.tellg() << endl;
		input.seekg(std::streampos(-WORD), ios_base::cur);
		//cout << "backward...at byte " << input.tellg() << endl;
		input.seekg(std::streampos(WORD), ios_base::cur);
		//cout << "foreward...at byte " << input.tellg() << endl;
	}

	printf(" nchannels %i nevents %i \n",nchannels,nevents );

	input.close();
}
