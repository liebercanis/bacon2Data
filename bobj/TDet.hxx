/**
** MG, March 2020 
**/
#ifndef TDET_DEFINED
#define TDET_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>
#include "TDetHit.hxx"

using namespace std;

// class to store info for the det  

class TDet: public TNamed{
	public:
		TDet(TString name=TString("UNKNOWN"), const TString title=TString(""));
    //TDet::~TDet(){}

		// data elements
    Long64_t event;
    Int_t    flags;
    Double_t ave;
    Double_t sigma;
    Int_t    nspe;
    Double_t qPrompt;
    Double_t qSum;
    Double_t hitPrompt;
    Double_t hitSum;

    std::vector<TDetHit> hits; 

    unsigned nhits() { return hits.size(); }

    void clear()
    {
      event=0;
      flags=0;
      ave=0;
      sigma=0;
      nspe=0;
      qPrompt=0;
      qSum=0;
      hitPrompt=0;
      hitSum=0;
      hits.clear();
    }

    ClassDef(TDet,3)
};
#endif

