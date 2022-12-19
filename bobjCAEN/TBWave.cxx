#include "TBWave.hxx"
ClassImp(TBWave)

TBWave::TBWave(TString detName) : TNamed(detName, detName)
{
  description = TString("none");
  cout << " new TBWave " << detName << "  " << this->GetName() << " " << description << endl;
  //clear();
}

// TBWave::~TBWave(){}

void TBWave::clear()
{
  fCurrent= -1;
  Channel= 0;
  Timestamp= 0;
  Board= 0;
  Energy= 0;
  EnergyShort= 0;
  Flags= 0;
  Probe= 0;

  Samples->Reset();
}
Long64_t TBWave::LoadTree(Long64_t entry)
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


Int_t TBWave::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

void TBWave::Init(TTree *tree)
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

