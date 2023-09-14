enum
{
  NFILES = 4
};

TFile *fout;

static bool goodChannel(int ic)
{
  bool val = true;
  if (ic == 5 || ic == 6 || ic == 8 || ic == 3 || ic == 9 || ic == 10 || ic == 11)
    val = false;
  return val;
}

void channelPlot1(int ic = 12)
{
    cout << "channelPlot1 for "<<ic << endl;
   

    int myColor[13] = {kRed, kBlue - 9, kGreen, kOrange, kBlue, kAzure + 3, kCyan - 3, kGreen - 4, kGreen, kSpring, kYellow, kBlue + 9, kTeal - 4};

    TString summaryFile[NFILES];
    TFile *fin[NFILES];
    std::vector<TH1D *> vhist;
    double dopant[NFILES];
    TString hname;
    TH1D *hWave;

    summaryFile[0] = TString("summary-type-1-dir-caenDataZeroPPM-2023-09-14-10-36.root"); // 05_22_2023 05_30_2023
    summaryFile[1] = TString("summary-type-1-dir-caenDataPointOnePPM-2023-09-14-13-28.root");
    summaryFile[2] = TString("summary-type-1-dir-caenDataPointTwoPPM-2023-09-14-13-28.root"); // 07_06_2023  07_07_2023
    summaryFile[3] = TString("summary-type-1-dir-caenDataPointFivePPM-2023-09-14-13-29.root");
    dopant[0] = 0.05;
    dopant[1] = 0.1;
    dopant[2] = 0.2;
    dopant[3] = 0.50;
    // dopant[4] = 1.0;
    // dopant[5] = 2.0;
    vhist.clear();
    for (int i = 0; i < NFILES; ++i)
    {
    fin[i] = new TFile(summaryFile[i], "readonly");
    TDirectory *runSumDir=NULL;
    fin[i]->GetObject("runSumDir", runSumDir);
    if(runSumDir==NULL) {
      printf(" no run sum in %s \n",summaryFile[i].Data());
      continue;
    }
    hname.Form("RunHitWaveChan%i", ic);
    printf("\t look for  chan %i %s \n", ic, hname.Data());
    TH1D *hWave=NULL;
    runSumDir->GetObject(hname, hWave);
    TString newName = Form("Ch%iDop%.2fPPM", ic, dopant[i]);
    vhist.push_back((TH1D *)hWave->Clone(newName));
    int jh = vhist.size() - 1;
    vhist[jh]->SetTitle(newName);
    vhist[jh]->SetLineColor(myColor[i]);
    vhist[jh]->SetMarkerColor(myColor[i]);
    vhist[jh]->SetMarkerStyle(20 + i);
    cout << " add " << vhist[jh]->GetName() << endl;
    fout->Add(vhist[jh]);
    //fin[i]->Close();
  }
  cout << "got " << vhist.size() << endl;
  

  TString canName;
  gStyle->SetOptStat(1);
  gStyle->SetOptLogy();
  canName.Form("channel%i", ic);
  TCanvas *can = new TCanvas(canName, canName);
  for (int i = vhist.size() - 1; i>-1 ; --i)
  {
    cout << "draw " << vhist[i]->GetName() << endl;
    if (i == vhist.size() - 1)
      vhist[i]->Draw();
    else
      vhist[i]->Draw("sames");
  }
  can->BuildLegend();
  can->Print(".pdf");
}

void channelPlot()
{
  fout = new TFile("channelPlotAll.root", "recreate");

  for (int ic = 0; ic < 13; ++ic){
  if (!goodChannel(ic))
      continue;
  channelPlot1(ic);
  }
  fout->ls();
  fout->Write();
}
