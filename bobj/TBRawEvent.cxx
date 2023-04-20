#include "TBRawEvent.hxx"
ClassImp(TBRawEvent)

TBRawEvent::TBRawEvent(unsigned ichannel): channel(ichannel)
{

  char dname[100];
  sprintf(dname, "chan%i", ichannel);
  this->SetName(dname);
  //cout << " new TBRawEvent " << this->GetName() << endl;
  clear();
}

//TBRawEvent::~TBRawEvent(){}
