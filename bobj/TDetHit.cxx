#include "TDetHit.hxx"
ClassImp(TDetHit)

TDetHit::TDetHit(): TNamed("TDetHit","TDetHit")
{
  clear();
}

//TDetHit::~TDetHit(){}

void TDetHit::clear()
{
  firstBin=0;
  lastBin=0;
  startTime=0;
  peakWidth=0;
  qpeak=0;
  peakt=0;
  peakMaxTime=0;
  peakBin=0;
  qsum=0;
  qerr=0;
  good=0;
  kind=-1; // unassigned
  digi.clear();
}

TH1D* TDetHit::plot(){
  if(digi.size()==0)
    return NULL;
  TString name;
  TString title;
  name.Form("hitPlot-start%i-qsum%f",firstBin,qsum);
  title.Form(" hit plot start %i qsum %f",firstBin,qsum);
  TH1D *h = new TH1D(name, title, digi.size(), 0, digi.size());
  for(int i=0; i<digi.size(); ++i)
    h->SetBinContent(i + 1, digi[i]);
  return h;
}
