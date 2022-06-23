#include "TDet.hxx"
ClassImp(TDet)

TDet::TDet(TString name, TString title): TNamed(name,title)
{
  clear();
}

