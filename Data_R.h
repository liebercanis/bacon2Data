//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May  9 14:11:11 2022 by ROOT version 6.12/06
// from TTree Data_R/CoMPASS RAW events TTree
// found on file: DataR_CH0@DT5730SB_2198_run_1.root
//////////////////////////////////////////////////////////

#ifndef Data_R_h
#define Data_R_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TArrayS.h"

class Data_R {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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

   Data_R(TTree *tree=0);
   virtual ~Data_R();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Long64_t nloop=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Data_R_cxx
Data_R::Data_R(TTree *tree) : fChain(0) 
{
  return;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DataR_CH0@DT5730SB_2198_run_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("DataR_CH0@DT5730SB_2198_run_1.root");
      }
      f->GetObject("Data_R",tree);

   }
   Init(tree);
}

Data_R::~Data_R()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Data_R::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Data_R::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Data_R::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Samples = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Channel", &Channel, &b_Channel);
   fChain->SetBranchAddress("Timestamp", &Timestamp, &b_Timestamp);
   fChain->SetBranchAddress("Board", &Board, &b_Board);
   fChain->SetBranchAddress("Energy", &Energy, &b_Energy);
   fChain->SetBranchAddress("EnergyShort", &EnergyShort, &b_EnergyShort);
   fChain->SetBranchAddress("Flags", &Flags, &b_Flags);
   fChain->SetBranchAddress("Probe", &Probe, &b_Probe);
   fChain->SetBranchAddress("Samples", &Samples, &b_Samples);
   Notify();
}

Bool_t Data_R::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Data_R::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Data_R::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Data_R_cxx
