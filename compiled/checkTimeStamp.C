void checkTimeStamp(TString fname = TString("run-01_12_2023-file0.root"))
{
   TString fullName = TString("data/rootData/")+fname;
   TFile* fin = new TFile(fullName,"READONLY");
   //fin->ls();
   TBFile *bf;
   fin->GetObject("tbfile",bf);
   bf->print();
}
