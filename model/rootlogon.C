// 
{
  printf("\n this is rootlogon for bacon \n");
  TString arch=gSystem->GetBuildArch();
  cout << " arch is " << arch << endl; 
  gSystem->AddIncludePath(" -I.");
  printf(" include path %s \n\n",gSystem->GetIncludePath());
  cout << "DYNAMIC PATH "  << gSystem->GetDynamicPath() << endl;
  cout << "LINKED LIBS "  << gSystem->GetLinkedLibs() << endl;
  gROOT->LoadMacro("util.C");
  gStyle->SetOptDate();
}
