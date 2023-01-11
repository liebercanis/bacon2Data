/**
** MG, July 2020 
** Modified for SIS data Jan 11, 2023
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
    TBRawEvent(unsigned  ichannel = 0);
	//		~TBRawEvent();
	// data elements
	unsigned   channel;
	unsigned   buffer;
	unsigned   trigger;
	Long64_t time;
	std::vector<unsigned short> rdigi;

	// methods
	void clear()
	{
		channel = 0;
		buffer = 0;
		trigger = 0;
		time = 0;
		rdigi.clear();
	}

	void print() {
		printf(" TBRawRunEvent chan %u buff %u trigger %u time %lld rdigi size %lu \n", 
		channel, buffer, trigger, time, rdigi.size()); 
	}

	ClassDef(TBRawEvent,2)
};
#endif

