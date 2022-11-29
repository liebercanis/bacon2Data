#include <numeric>
#include "TBRawEvent.hxx"
#include "TBRun.hxx"
TFile *fout;
TFile *fin;
TTree *rtree;
TBRawRun *tbrun;
TObjArray *brList;
vector<TBRawEvent *> detList;
vector<int> vchan;
vector<double> digi;
vector<TH1D *> noiseHist;
vector<TH1D *> skewHist;
vector<TH1D *> valHist;
vector<TH1D *> sumWave;
// TTree *ftree;

bool openFile(TString tag)
{
  // open ouput file and make some histograms
  TString fileName; fileName.Form("%s.root",tag.Data());
  printf(" looking for file %s\n",fileName.Data());
  fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return false;
  }

  rtree = NULL;
  brList = NULL;
  fin->GetObject("RawTree", rtree);
  if(!rtree) {
    printf("!! RawTree not found \n");
    return false;
  }

  printf(" RawTree entries = %llu\n", rtree->GetEntries() );
  brList = rtree->GetListOfBranches();
  brList->ls();

  TIter next(brList);
  TBranchElement *aBranch=NULL;
  detList.clear();

  while( ( aBranch = (TBranchElement *) next() ) ) {
    TString brname = TString(aBranch->GetName()) ;
    int chan = TString(brname(brname.First('n')+1,brname.Length() - brname.First('n'))).Atoi();
    vchan.push_back(chan);
    aBranch->SetAddress(0);
  }
  printf("total channels  =  %lu \n", vchan.size());
  return true;
}


void anaRun(TString tag= TString("No_Doping_1_Digitizer"))
{
  // hist to fit to baseline offset
  TH1D *hEvGaus = new TH1D("baseline","baseline",200,-100,100);

   if(!openFile(tag)) return;
  TString outFileName = TString("anaRun") + tag + TString(".root");
  fout = new TFile(outFileName,"recreate");

  for (unsigned ichan = 0; ichan < vchan.size();  ++ichan ){
    noiseHist.push_back(new TH1D(Form("NoiseChan%i", ichan),Form("NoiseChan%i", ichan), 100,0,10));
    skewHist.push_back(new TH1D(Form("SkewChan%i", ichan),Form("SkewChan%i", ichan), 100,-20,20));
    valHist.push_back(new TH1D(Form("ValChan%i", ichan),Form("ValChan%i", ichan), 2200,-20,200));
    sumWave.push_back(new TH1D(Form("sumWave%i", ichan),Form("sunWave%i", ichan), 4100,0,4100));

  }

  Long64_t nentries = rtree->GetEntries();
  printf(" loop over %llu entries \n",nentries);
  for (Long64_t entry = 0; entry < nentries; ++entry)
  {
    if(entry/100*100==entry) cout << " .... entry " << entry << endl;
    rtree->GetEntry(entry);

    // find branches and cast as TBRawEvent. each is a channel
    TIter next(brList);
    TBranchElement *aBranch = NULL;
    detList.clear();
    while ((aBranch = (TBranchElement *) next())) {
      TBRawEvent *det = (TBRawEvent *) aBranch->GetObject();
      detList.push_back(det);
    }
    
    // loop over channels
    for (int ichan = 0; ichan < detList.size(); ++ichan) 
    {
      int nbins = detList[ichan]->rdigi.size();
      if (nbins < 1)
        continue;

      /**
      cout << detList[ichan]->channel << " " << detList[ichan]->time << " " << detList[ichan]->rdigi.size() << "; ";
      cout << endl;
      **/

      TString hname;
      hname.Form("EvRawWaveEv%ich%i", int(entry), ichan);
      TH1D *hEvRawWave = NULL;
      
      // simple baseline
      double base = std::accumulate(detList[ichan]->rdigi.begin(), detList[ichan]->rdigi.end(),0)/ detList[ichan]->rdigi.size();

      // baseline correction from fitted Gaussian
      hEvGaus->Reset();
      for (unsigned j = 0; j < detList[ichan]->rdigi.size(); ++j)
        hEvGaus->Fill((double)detList[ichan]->rdigi[j] - base);
  
      TFitResultPtr fitptr = hEvGaus->Fit("gaus","QOS");
      TF1 *gfit = (TF1 *)hEvGaus->GetListOfFunctions()->FindObject("gaus");
      double ave = gfit->GetParameter(1);
      double sigma  = gfit->GetParameter(2);
      double skew   = hEvGaus->GetSkewness();
      // fitptr->Print();
      //cout << hname << " " <<   " ave " << ave << " sigma " << sigma << " skew " << skew <<  endl;
      noiseHist[ichan]->Fill(sigma);
      skewHist[ ichan]->Fill(skew);

      if((abs(skew)>1)) { 
        hEvRawWave = new TH1D(hname, hname, nbins, 0, nbins);
        for (unsigned j = 0; j < detList[ichan]->rdigi.size(); ++j) {
          double val = (double)detList[ichan]->rdigi[j] - base - ave;
          hEvRawWave->SetBinContent(j + 1, val);
          sumWave[ichan]->SetBinContent(j + 1,  sumWave[ichan]->GetBinContent(j+ 1) + val );
          valHist[ichan]->Fill(val);
        }
      }
    }
  }
  fout->ls();
  fout->Write();
  return;
}
