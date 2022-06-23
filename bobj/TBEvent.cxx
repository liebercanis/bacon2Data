#include "TBEvent.hxx"
ClassImp(TBEvent)

TBEvent::TBEvent(TString runName ): TNamed(runName,runName)
{
  clear();
}


//TBEvent::~TBEvent(){}

void TBEvent::clear()
{
  event=0;
  trigTime=0; 
  energy=0;
}


