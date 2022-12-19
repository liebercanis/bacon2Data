#include "TBEvent.hxx"
ClassImp(TBEvent)

TBEvent::TBEvent(TString runName ): TNamed(runName,runName)
{
  clear();
  cout << " TBEvent instance " << runName << endl;
}
//TBEvent::~TBEvent(){}

