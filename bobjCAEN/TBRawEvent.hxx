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
	// methods
	void clear();
	void print();
	// data elements
	int length;
	int boardID;
	int channel;
	int event;
	Long64_t time;
	int dcOffset;

	void printHeader() {
		printf(" TBRawRunEvent %i chan %i length %i ID %i off %i time %llu \n", 
		event, channel, length, boardID,dcOffset, time); 
	}

	std::vector<UShort_t> rdigi;
	ClassDef(TBRawEvent,1)
};
#endif

