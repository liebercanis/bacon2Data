#include "TBRawEvent.hxx"
ClassImp(TBRawEvent)

TBRawEvent::TBRawEvent(int ichannel): channel(ichannel)
{

  char dname[100];
  sprintf(dname, "chan%i", ichannel);
  this->SetName(dname);
  cout << " new TBRawEvent " << this->GetName() << endl;
  clear();
}

//TBRawEvent::~TBRawEvent(){}

void TBRawEvent::clear()
{
  length=0;
  boardID=0;
  event=0;
  time=0;
  dcOffset=0;
  rdigi.clear();	 
}

