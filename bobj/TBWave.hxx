/**
** MG, May 2022
**/
#ifndef TBWAVE_DEFINED
#define TBWAVE_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <TArrayS.h>
#include <TTree.h>
#include <vector>

using namespace std;

// class to store info for the event

class TBWave: public TNamed {
	public:
   TBWave(TString detName = "SIPM0");
	 //~TBWave();
	 // methods
   Long64_t LoadTree(Long64_t entry);
   void Init(TTree *tree);
   Int_t  GetEntry(Long64_t entry);


	void clear();
	void print()
	{ 
    cout << this->GetName() << "  " << description << endl; 
  }

  Bool_t Notify()
  {
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
  }


  // data elements
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t fCurrent;

  TString description;

  UShort_t        Channel;
  ULong64_t       Timestamp;
  UShort_t        Board;
  UShort_t        Energy;
  UShort_t        EnergyShort;
  UInt_t          Flags;
  Int_t           Probe;
  TArrayS         *Samples;

  // List of branches
  TBranch        *b_Channel;   //!
  TBranch        *b_Timestamp;   //!
  TBranch        *b_Board;   //!
  TBranch        *b_Energy;   //!
  TBranch        *b_EnergyShort;   //!
  TBranch        *b_Flags;   //!
  TBranch        *b_Probe;   //!
  TBranch        *b_Samples;   //!

  ClassDef(TBWave,1)
};
#endif

