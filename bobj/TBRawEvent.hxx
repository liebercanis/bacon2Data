/**
** MG, July 2020 
**/
#ifndef TBRAWEVENT_DEFINED
#define TBRAWEVENT_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>

using namespace std;

// class to store info for the event

class TBRawEvent: public TNamed {
	public:
    TBRawEvent(int ichannel = 0);
	//		~TBRawEvent();
	// data elements
	unsigned   channel;
	unsigned   samples;
	unsigned   startEnergy;
	unsigned   maxEnergy;
	Long64_t time;
	std::vector<unsigned> rdigi;

	// methods
	void clear()
	{
		channel = 0;
		samples = 0;
		startEnergy = 0;
		maxEnergy = 0;
		time = 0;
		rdigi.clear();
	}

	void print() {
		printf(" TBRawRunEvent chan %u samples %u time %llu startEnergy  %u time maxEnergy %u rdigi size %lu \n", 
		channel, samples, time, startEnergy, maxEnergy,rdigi.size()); 
	}

	ClassDef(TBRawEvent,1)
};
#endif

