// 
{
  printf("\n this is rootlogon for bacon \n");
  TString arch=gSystem->GetBuildArch();
  cout << " arch is " << arch << endl; 
  gSystem->AddIncludePath(" -I. -I./bobj/");
  gSystem->AddDynamicPath("./bobj/");
  printf(" include path %s \n\n",gSystem->GetIncludePath());
  cout << "DYNAMIC PATH "  << gSystem->GetDynamicPath() << endl;
  cout << "LINKED LIBS "  << gSystem->GetLinkedLibs() << endl;
  gROOT->LoadMacro("util.C");
  //gROOT->LoadMacro("Data_R.C");
  printf("\n this is ROOT \n");

  gStyle->SetOptDate();
  int iload = gSystem->Load("bobj/libBaconAna.so");
  printf(" loaded libBaconAna = %i zero is success! \n",iload);
}
