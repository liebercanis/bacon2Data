#include "TBRun.hxx"
ClassImp(TBRun)

TBRun::TBRun(TString tag): TNamed(TString("TBrun"),tag)
{
  btree = new TTree("RunTree", " bacon data ");
  cout << "Instance of TBRun tree " << btree->GetName() << endl;
  btree->GetListOfBranches()->ls();
  clear();
}
// TBRun::~TBRun(){}

void TBRun::clear()
{
  for (unsigned i = 0; i < detList.size(); ++i)
    detList[i]->clear();
}