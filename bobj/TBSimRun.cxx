#include "TBSimRun.hxx"
ClassImp(TBSimRun)

    TBSimRun::TBSimRun(TString runName) : TNamed(runName, runName)
// TBSimRun::TBSimRun()
{
  btree = new TTree("SimTree", " bacon sim data ");
}
