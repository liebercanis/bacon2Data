#include "TDet.hxx"
ClassImp(TDet)
TDet::TDet()
{
  int ichan = -1;
  int ilevel = -1;
  this->SetName(Form("chan%i", ichan));
  this->SetTitle(Form("channel %i level %i", ichan, ilevel));
  channel = ichan;
  level = ilevel;
  clear();
}
TDet::TDet(int ichan, int ilevel)
{
  this->SetName(Form("chan%i",ichan));
  this->SetTitle(Form("channel %i level %i", ichan, ilevel));
  //printf(" TDET INSTANCE %i %i %s \n", ichan, ilevel, this->GetName());
  channel = ichan;
  level = ilevel;
  clear();
}

