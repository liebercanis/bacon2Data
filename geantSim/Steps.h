//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul  4 13:11:17 2025 by ROOT version 6.34.08
// from TTree Processed_ArgonSteps_Summary/
// found on file: 20250627_BACONSim_ProcessedSimulation_250kevents.root
//////////////////////////////////////////////////////////

#ifndef Steps_h
#define Steps_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Steps
{
public:
   TTree *fChain;  //! pointer to the analyzed TTree or TChain
   Int_t fCurrent; //! current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t index;
   Double_t stepenergies_keV;
   Long64_t eventnumbers;
   Char_t processinfomation[12];
   Char_t stepvolume[12];
   Long64_t parentID;
   Bool_t TriggerSiPM1_LOS;
   Double_t TriggerSiPM1_distance_mm;
   Double_t TriggerSiPM1_fluxfraction;
   Double_t TriggerSiPM1_angletonormal_rad;
   Bool_t TriggerSiPM2_LOS;
   Double_t TriggerSiPM2_distance_mm;
   Double_t TriggerSiPM2_fluxfraction;
   Double_t TriggerSiPM2_angletonormal_rad;
   Bool_t TriggerSiPM3_LOS;
   Double_t TriggerSiPM3_distance_mm;
   Double_t TriggerSiPM3_fluxfraction;
   Double_t TriggerSiPM3_angletonormal_rad;
   Bool_t Rod1_Sensor6_LOS;
   Double_t Rod1_Sensor6_distance_mm;
   Double_t Rod1_Sensor6_fluxfraction;
   Double_t Rod1_Sensor6_angletonormal_rad;
   Bool_t Rod1_Sensor3_LOS;
   Double_t Rod1_Sensor3_distance_mm;
   Double_t Rod1_Sensor3_fluxfraction;
   Double_t Rod1_Sensor3_angletonormal_rad;
   Bool_t Rod1_Sensor0_LOS;
   Double_t Rod1_Sensor0_distance_mm;
   Double_t Rod1_Sensor0_fluxfraction;
   Double_t Rod1_Sensor0_angletonormal_rad;
   Bool_t Rod2_Sensor2_LOS;
   Double_t Rod2_Sensor2_distance_mm;
   Double_t Rod2_Sensor2_fluxfraction;
   Double_t Rod2_Sensor2_angletonormal_rad;
   Bool_t Rod2_Sensor4_LOS;
   Double_t Rod2_Sensor4_distance_mm;
   Double_t Rod2_Sensor4_fluxfraction;
   Double_t Rod2_Sensor4_angletonormal_rad;
   Bool_t Rod2_Sensor1_LOS;
   Double_t Rod2_Sensor1_distance_mm;
   Double_t Rod2_Sensor1_fluxfraction;
   Double_t Rod2_Sensor1_angletonormal_rad;
   Bool_t Rod3_Sensor3_LOS;
   Double_t Rod3_Sensor3_distance_mm;
   Double_t Rod3_Sensor3_fluxfraction;
   Double_t Rod3_Sensor3_angletonormal_rad;
   Bool_t Rod3_Sensor5_LOS;
   Double_t Rod3_Sensor5_distance_mm;
   Double_t Rod3_Sensor5_fluxfraction;
   Double_t Rod3_Sensor5_angletonormal_rad;
   Bool_t Rod3_Sensor2_LOS;
   Double_t Rod3_Sensor2_distance_mm;
   Double_t Rod3_Sensor2_fluxfraction;
   Double_t Rod3_Sensor2_angletonormal_rad;
   Bool_t PMT_Sensor12_LOS;
   Double_t PMT_Sensor12_distance_mm;
   Double_t PMT_Sensor12_fluxfraction;
   Double_t PMT_Sensor12_angletonormal_rad;
   Bool_t OR_LOSdata_allSiPMs;
   Double_t Fluxfraction_allSiPMs;
   Bool_t AND_LOSdata_allSiPMs;
   Double_t position_x_mm;
   Double_t position_y_mm;
   Double_t position_z_mm;
   Double_t position_r_mm;
   Double_t position_theta_rad;
   Double_t position_phi_rad;

   // List of branches
   TBranch *b_index;                          //!
   TBranch *b_stepenergies_keV;               //!
   TBranch *b_eventnumbers;                   //!
   TBranch *b_processinfomation;              //!
   TBranch *b_stepvolume;                     //!
   TBranch *b_parentID;                       //!
   TBranch *b_TriggerSiPM1_LOS;               //!
   TBranch *b_TriggerSiPM1_distance_mm;       //!
   TBranch *b_TriggerSiPM1_fluxfraction;      //!
   TBranch *b_TriggerSiPM1_angletonormal_rad; //!
   TBranch *b_TriggerSiPM2_LOS;               //!
   TBranch *b_TriggerSiPM2_distance_mm;       //!
   TBranch *b_TriggerSiPM2_fluxfraction;      //!
   TBranch *b_TriggerSiPM2_angletonormal_rad; //!
   TBranch *b_TriggerSiPM3_LOS;               //!
   TBranch *b_TriggerSiPM3_distance_mm;       //!
   TBranch *b_TriggerSiPM3_fluxfraction;      //!
   TBranch *b_TriggerSiPM3_angletonormal_rad; //!
   TBranch *b_Rod1_Sensor6_LOS;               //!
   TBranch *b_Rod1_Sensor6_distance_mm;       //!
   TBranch *b_Rod1_Sensor6_fluxfraction;      //!
   TBranch *b_Rod1_Sensor6_angletonormal_rad; //!
   TBranch *b_Rod1_Sensor3_LOS;               //!
   TBranch *b_Rod1_Sensor3_distance_mm;       //!
   TBranch *b_Rod1_Sensor3_fluxfraction;      //!
   TBranch *b_Rod1_Sensor3_angletonormal_rad; //!
   TBranch *b_Rod1_Sensor0_LOS;               //!
   TBranch *b_Rod1_Sensor0_distance_mm;       //!
   TBranch *b_Rod1_Sensor0_fluxfraction;      //!
   TBranch *b_Rod1_Sensor0_angletonormal_rad; //!
   TBranch *b_Rod2_Sensor2_LOS;               //!
   TBranch *b_Rod2_Sensor2_distance_mm;       //!
   TBranch *b_Rod2_Sensor2_fluxfraction;      //!
   TBranch *b_Rod2_Sensor2_angletonormal_rad; //!
   TBranch *b_Rod2_Sensor4_LOS;               //!
   TBranch *b_Rod2_Sensor4_distance_mm;       //!
   TBranch *b_Rod2_Sensor4_fluxfraction;      //!
   TBranch *b_Rod2_Sensor4_angletonormal_rad; //!
   TBranch *b_Rod2_Sensor1_LOS;               //!
   TBranch *b_Rod2_Sensor1_distance_mm;       //!
   TBranch *b_Rod2_Sensor1_fluxfraction;      //!
   TBranch *b_Rod2_Sensor1_angletonormal_rad; //!
   TBranch *b_Rod3_Sensor3_LOS;               //!
   TBranch *b_Rod3_Sensor3_distance_mm;       //!
   TBranch *b_Rod3_Sensor3_fluxfraction;      //!
   TBranch *b_Rod3_Sensor3_angletonormal_rad; //!
   TBranch *b_Rod3_Sensor5_LOS;               //!
   TBranch *b_Rod3_Sensor5_distance_mm;       //!
   TBranch *b_Rod3_Sensor5_fluxfraction;      //!
   TBranch *b_Rod3_Sensor5_angletonormal_rad; //!
   TBranch *b_Rod3_Sensor2_LOS;               //!
   TBranch *b_Rod3_Sensor2_distance_mm;       //!
   TBranch *b_Rod3_Sensor2_fluxfraction;      //!
   TBranch *b_Rod3_Sensor2_angletonormal_rad; //!
   TBranch *b_PMT_Sensor12_LOS;               //!
   TBranch *b_PMT_Sensor12_distance_mm;       //!
   TBranch *b_PMT_Sensor12_fluxfraction;      //!
   TBranch *b_PMT_Sensor12_angletonormal_rad; //!
   TBranch *b_OR_LOSdata_allSiPMs;            //!
   TBranch *b_Fluxfraction_allSiPMs;          //!
   TBranch *b_AND_LOSdata_allSiPMs;           //!
   TBranch *b_position_x_mm;                  //!
   TBranch *b_position_y_mm;                  //!
   TBranch *b_position_z_mm;                  //!
   TBranch *b_position_r_mm;                  //!
   TBranch *b_position_theta_rad;             //!
   TBranch *b_position_phi_rad;               //!

   Steps(TTree *tree = 0);
   virtual ~Steps();
   virtual Int_t Cut(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void Init(TTree *tree);
   virtual void Loop(Long64_t maxEntry);
   virtual bool Notify();
   virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef Steps_cxx
Steps::Steps(TTree *tree) : fChain(0)
{
   // if parameter tree is not specified (or zero), connect the file
   // used to generate this class and read the Tree.
   if (tree == 0)
   {
      TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject("20250627_BACONSim_ProcessedSimulation_250kevents.root");
      if (!f || !f->IsOpen())
      {
         f = new TFile("20250627_BACONSim_ProcessedSimulation_250kevents.root");
      }
      f->GetObject("Processed_ArgonSteps_Summary", tree);
   }
   Init(tree);
}

Steps::~Steps()
{
   if (!fChain)
      return;
   delete fChain->GetCurrentFile();
}

Int_t Steps::GetEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain)
      return 0;
   return fChain->GetEntry(entry);
}
Long64_t Steps::LoadTree(Long64_t entry)
{
   // Set the environment to read one entry
   if (!fChain)
      return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0)
      return centry;
   if (fChain->GetTreeNumber() != fCurrent)
   {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Steps::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree)
      return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("index", &index, &b_index);
   fChain->SetBranchAddress("stepenergies_keV", &stepenergies_keV, &b_stepenergies_keV);
   fChain->SetBranchAddress("eventnumbers", &eventnumbers, &b_eventnumbers);
   fChain->SetBranchAddress("processinfomation", processinfomation, &b_processinfomation);
   fChain->SetBranchAddress("stepvolume", stepvolume, &b_stepvolume);
   fChain->SetBranchAddress("parentID", &parentID, &b_parentID);
   fChain->SetBranchAddress("TriggerSiPM1_LOS", &TriggerSiPM1_LOS, &b_TriggerSiPM1_LOS);
   fChain->SetBranchAddress("TriggerSiPM1_distance_mm", &TriggerSiPM1_distance_mm, &b_TriggerSiPM1_distance_mm);
   fChain->SetBranchAddress("TriggerSiPM1_fluxfraction", &TriggerSiPM1_fluxfraction, &b_TriggerSiPM1_fluxfraction);
   fChain->SetBranchAddress("TriggerSiPM1_angletonormal_rad", &TriggerSiPM1_angletonormal_rad, &b_TriggerSiPM1_angletonormal_rad);
   fChain->SetBranchAddress("TriggerSiPM2_LOS", &TriggerSiPM2_LOS, &b_TriggerSiPM2_LOS);
   fChain->SetBranchAddress("TriggerSiPM2_distance_mm", &TriggerSiPM2_distance_mm, &b_TriggerSiPM2_distance_mm);
   fChain->SetBranchAddress("TriggerSiPM2_fluxfraction", &TriggerSiPM2_fluxfraction, &b_TriggerSiPM2_fluxfraction);
   fChain->SetBranchAddress("TriggerSiPM2_angletonormal_rad", &TriggerSiPM2_angletonormal_rad, &b_TriggerSiPM2_angletonormal_rad);
   fChain->SetBranchAddress("TriggerSiPM3_LOS", &TriggerSiPM3_LOS, &b_TriggerSiPM3_LOS);
   fChain->SetBranchAddress("TriggerSiPM3_distance_mm", &TriggerSiPM3_distance_mm, &b_TriggerSiPM3_distance_mm);
   fChain->SetBranchAddress("TriggerSiPM3_fluxfraction", &TriggerSiPM3_fluxfraction, &b_TriggerSiPM3_fluxfraction);
   fChain->SetBranchAddress("TriggerSiPM3_angletonormal_rad", &TriggerSiPM3_angletonormal_rad, &b_TriggerSiPM3_angletonormal_rad);
   fChain->SetBranchAddress("Rod1_Sensor6_LOS", &Rod1_Sensor6_LOS, &b_Rod1_Sensor6_LOS);
   fChain->SetBranchAddress("Rod1_Sensor6_distance_mm", &Rod1_Sensor6_distance_mm, &b_Rod1_Sensor6_distance_mm);
   fChain->SetBranchAddress("Rod1_Sensor6_fluxfraction", &Rod1_Sensor6_fluxfraction, &b_Rod1_Sensor6_fluxfraction);
   fChain->SetBranchAddress("Rod1_Sensor6_angletonormal_rad", &Rod1_Sensor6_angletonormal_rad, &b_Rod1_Sensor6_angletonormal_rad);
   fChain->SetBranchAddress("Rod1_Sensor3_LOS", &Rod1_Sensor3_LOS, &b_Rod1_Sensor3_LOS);
   fChain->SetBranchAddress("Rod1_Sensor3_distance_mm", &Rod1_Sensor3_distance_mm, &b_Rod1_Sensor3_distance_mm);
   fChain->SetBranchAddress("Rod1_Sensor3_fluxfraction", &Rod1_Sensor3_fluxfraction, &b_Rod1_Sensor3_fluxfraction);
   fChain->SetBranchAddress("Rod1_Sensor3_angletonormal_rad", &Rod1_Sensor3_angletonormal_rad, &b_Rod1_Sensor3_angletonormal_rad);
   fChain->SetBranchAddress("Rod1_Sensor0_LOS", &Rod1_Sensor0_LOS, &b_Rod1_Sensor0_LOS);
   fChain->SetBranchAddress("Rod1_Sensor0_distance_mm", &Rod1_Sensor0_distance_mm, &b_Rod1_Sensor0_distance_mm);
   fChain->SetBranchAddress("Rod1_Sensor0_fluxfraction", &Rod1_Sensor0_fluxfraction, &b_Rod1_Sensor0_fluxfraction);
   fChain->SetBranchAddress("Rod1_Sensor0_angletonormal_rad", &Rod1_Sensor0_angletonormal_rad, &b_Rod1_Sensor0_angletonormal_rad);
   fChain->SetBranchAddress("Rod2_Sensor2_LOS", &Rod2_Sensor2_LOS, &b_Rod2_Sensor2_LOS);
   fChain->SetBranchAddress("Rod2_Sensor2_distance_mm", &Rod2_Sensor2_distance_mm, &b_Rod2_Sensor2_distance_mm);
   fChain->SetBranchAddress("Rod2_Sensor2_fluxfraction", &Rod2_Sensor2_fluxfraction, &b_Rod2_Sensor2_fluxfraction);
   fChain->SetBranchAddress("Rod2_Sensor2_angletonormal_rad", &Rod2_Sensor2_angletonormal_rad, &b_Rod2_Sensor2_angletonormal_rad);
   fChain->SetBranchAddress("Rod2_Sensor4_LOS", &Rod2_Sensor4_LOS, &b_Rod2_Sensor4_LOS);
   fChain->SetBranchAddress("Rod2_Sensor4_distance_mm", &Rod2_Sensor4_distance_mm, &b_Rod2_Sensor4_distance_mm);
   fChain->SetBranchAddress("Rod2_Sensor4_fluxfraction", &Rod2_Sensor4_fluxfraction, &b_Rod2_Sensor4_fluxfraction);
   fChain->SetBranchAddress("Rod2_Sensor4_angletonormal_rad", &Rod2_Sensor4_angletonormal_rad, &b_Rod2_Sensor4_angletonormal_rad);
   fChain->SetBranchAddress("Rod2_Sensor1_LOS", &Rod2_Sensor1_LOS, &b_Rod2_Sensor1_LOS);
   fChain->SetBranchAddress("Rod2_Sensor1_distance_mm", &Rod2_Sensor1_distance_mm, &b_Rod2_Sensor1_distance_mm);
   fChain->SetBranchAddress("Rod2_Sensor1_fluxfraction", &Rod2_Sensor1_fluxfraction, &b_Rod2_Sensor1_fluxfraction);
   fChain->SetBranchAddress("Rod2_Sensor1_angletonormal_rad", &Rod2_Sensor1_angletonormal_rad, &b_Rod2_Sensor1_angletonormal_rad);
   fChain->SetBranchAddress("Rod3_Sensor3_LOS", &Rod3_Sensor3_LOS, &b_Rod3_Sensor3_LOS);
   fChain->SetBranchAddress("Rod3_Sensor3_distance_mm", &Rod3_Sensor3_distance_mm, &b_Rod3_Sensor3_distance_mm);
   fChain->SetBranchAddress("Rod3_Sensor3_fluxfraction", &Rod3_Sensor3_fluxfraction, &b_Rod3_Sensor3_fluxfraction);
   fChain->SetBranchAddress("Rod3_Sensor3_angletonormal_rad", &Rod3_Sensor3_angletonormal_rad, &b_Rod3_Sensor3_angletonormal_rad);
   fChain->SetBranchAddress("Rod3_Sensor5_LOS", &Rod3_Sensor5_LOS, &b_Rod3_Sensor5_LOS);
   fChain->SetBranchAddress("Rod3_Sensor5_distance_mm", &Rod3_Sensor5_distance_mm, &b_Rod3_Sensor5_distance_mm);
   fChain->SetBranchAddress("Rod3_Sensor5_fluxfraction", &Rod3_Sensor5_fluxfraction, &b_Rod3_Sensor5_fluxfraction);
   fChain->SetBranchAddress("Rod3_Sensor5_angletonormal_rad", &Rod3_Sensor5_angletonormal_rad, &b_Rod3_Sensor5_angletonormal_rad);
   fChain->SetBranchAddress("Rod3_Sensor2_LOS", &Rod3_Sensor2_LOS, &b_Rod3_Sensor2_LOS);
   fChain->SetBranchAddress("Rod3_Sensor2_distance_mm", &Rod3_Sensor2_distance_mm, &b_Rod3_Sensor2_distance_mm);
   fChain->SetBranchAddress("Rod3_Sensor2_fluxfraction", &Rod3_Sensor2_fluxfraction, &b_Rod3_Sensor2_fluxfraction);
   fChain->SetBranchAddress("Rod3_Sensor2_angletonormal_rad", &Rod3_Sensor2_angletonormal_rad, &b_Rod3_Sensor2_angletonormal_rad);
   fChain->SetBranchAddress("PMT_Sensor12_LOS", &PMT_Sensor12_LOS, &b_PMT_Sensor12_LOS);
   fChain->SetBranchAddress("PMT_Sensor12_distance_mm", &PMT_Sensor12_distance_mm, &b_PMT_Sensor12_distance_mm);
   fChain->SetBranchAddress("PMT_Sensor12_fluxfraction", &PMT_Sensor12_fluxfraction, &b_PMT_Sensor12_fluxfraction);
   fChain->SetBranchAddress("PMT_Sensor12_angletonormal_rad", &PMT_Sensor12_angletonormal_rad, &b_PMT_Sensor12_angletonormal_rad);
   fChain->SetBranchAddress("OR_LOSdata_allSiPMs", &OR_LOSdata_allSiPMs, &b_OR_LOSdata_allSiPMs);
   fChain->SetBranchAddress("Fluxfraction_allSiPMs", &Fluxfraction_allSiPMs, &b_Fluxfraction_allSiPMs);
   fChain->SetBranchAddress("AND_LOSdata_allSiPMs", &AND_LOSdata_allSiPMs, &b_AND_LOSdata_allSiPMs);
   fChain->SetBranchAddress("position_x_mm", &position_x_mm, &b_position_x_mm);
   fChain->SetBranchAddress("position_y_mm", &position_y_mm, &b_position_y_mm);
   fChain->SetBranchAddress("position_z_mm", &position_z_mm, &b_position_z_mm);
   fChain->SetBranchAddress("position_r_mm", &position_r_mm, &b_position_r_mm);
   fChain->SetBranchAddress("position_theta_rad", &position_theta_rad, &b_position_theta_rad);
   fChain->SetBranchAddress("position_phi_rad", &position_phi_rad, &b_position_phi_rad);
   Notify();
}

bool Steps::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void Steps::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain)
      return;
   fChain->Show(entry);
}
Int_t Steps::Cut(Long64_t entry)
{
   // This function may be called from Loop.
   // returns  1 if entry is accepted.
   // returns -1 otherwise.
   return 1;
}

void Steps::Loop(Long64_t maxEntry = 0)
{
   //   In a ROOT session, you can do:
   //      root> .L Steps.C
   //      root> Steps t
   //      root> t.GetEntry(12); // Fill t data members with entry number 12
   //      root> t.Show();       // Show values of entry 12
   //      root> t.Show(16);     // Read and show values of entry 16
   //      root> t.Loop();       // Loop on all entries
   //

   //     This is the loop skeleton where:
   //    jentry is the global entry number in the chain
   //    ientry is the entry number in the current Tree
   //  Note that the argument to GetEntry must be:
   //    jentry for TChain::GetEntry
   //    ientry for TTree::GetEntry and TBranch::GetEntry
   //
   //       To read only selected branches, Insert statements like:
   // METHOD1:
   //    fChain->SetBranchStatus("*",0);  // disable all branches
   //    fChain->SetBranchStatus("branchname",1);  // activate branchname
   // METHOD2: replace line
   //    fChain->GetEntry(jentry);       //read all branches
   // by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0)
      return;

   Long64_t nentries = fChain->GetEntriesFast();

   if (maxEntry > 0)
      nentries = maxEntry;

   printf("Steps::Loop over %lld \n", nentries);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      printf("%lld event %lld \n", jentry, eventnumbers);
   }
}

#endif // #ifdef Steps_cxx
