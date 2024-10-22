#include "TBRawRun.hxx"
ClassImp(TBRawRun)

    TBRawRun::TBRawRun(TString runName) : TNamed(runName, runName)
{
  btree = new TTree("RawTree", " bacon raw data ");
  eventData = new TBEventData();
  eventSummary = new TBRawSummary();
  btree->Branch("eventData", eventData);
}

void TBRawRun::clear()
{
  detListClear();
  eventSummary->clear();
}

void TBRawRun::detListClear()
{
  for (unsigned i = 0; i < detList.size(); ++i)
    detList[i]->clear();
}
