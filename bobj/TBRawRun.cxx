#include "TBRawRun.hxx"
ClassImp(TBRawRun)

TBRawRun::TBRawRun(TString runName ): TNamed(runName,runName){
  btree = new TTree("RawTree"," bacon raw data " );
}

void TBRawRun::clear()
{
  detListClear();
}

void TBRawRun::detListClear() 
{
  for(unsigned i=0; i<detList.size(); ++i ) detList[i]->clear(); 
}


