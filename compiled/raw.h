//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 10 15:08:07 2023 by ROOT version 6.24/06
// from TTree raw/raw data
// found on file: sis3316_test_data_2000s.root
//////////////////////////////////////////////////////////

#ifndef raw_h
#define raw_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "c++/v1/vector"

class raw {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           channel;
   Int_t           buff;
   Int_t           event;
   vector<unsigned short> *data;

   // List of branches
   TBranch        *b_channel;   //!
   TBranch        *b_buff;   //!
   TBranch        *b_event;   //!
   TBranch        *b_data;   //!

   raw(TTree *tree=0);
   virtual ~raw();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   TH1I*  histEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef raw_cxx
raw::raw(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sis3316_test_data_2000s.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sis3316_test_data_2000s.root");
      }
      f->GetObject("raw",tree);

   }
   Init(tree);
}

raw::~raw()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}
TH1I* raw::histEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain)
      return NULL;
   
   fChain->GetEntry(entry);
   TString hname;
   hname.Form("Ch%iBuff%iEvent%i", channel, buff, event);
   TH1I *h = new TH1I(hname, hname, data->size(), 0, data->size());
   for (unsigned  i = 0; i < data->size(); ++i ){
      h->SetBinContent(i + 1, int(data->at(i)));
   }
   return h;
}

Int_t raw::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t raw::LoadTree(Long64_t entry)
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

void raw::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   data = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("buff", &buff, &b_buff);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("data", &data, &b_data);
   Notify();
}

Bool_t raw::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void raw::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t raw::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef raw_cxx
