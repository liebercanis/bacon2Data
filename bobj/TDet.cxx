#include "TDet.hxx"
ClassImp(TDet)

TDet::TDet(int ichan, int ilevel)
{
  this->SetName(Form("chan%i",ichan));
  this->SetTitle(Form("channel %i level %i", ichan, ilevel));
  channel = ichan;
  level = ilevel;
  clear();
}

