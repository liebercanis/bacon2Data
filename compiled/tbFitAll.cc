#include <iostream>
#include <fstream>
#include "modelAFit.hh"
// time is in microseconds
using namespace TMath;
TFile *fin;
TDirectory *runSumDir;
std::vector<TH1D *> vhist;
std::vector<TH1D *> hffit;
std::vector<TH1D *> hmodel;
std::vector<TH1D *> hffitPmt;
std::vector<TH1D *> hffitChan;
TString histSet;
double dopant[10];
double ylow = 1300;
double yhigh = 4000;

void fillHistogram(TF1* func, TH1D* hist) {
  if(!func) {
    printf("!!! fillHistogram called with null func !!!!! \n");
    return;
  }
  if(hist==NULL)
    return;
  for (int ib = 1; ib < hist->GetNbinsX(); ++ib)
  {
    double xbin = hist->GetBinCenter(ib);
    double fbin = func->Eval(xbin);
    hist->SetBinContent(ib,fbin);
    hist->SetBinError(ib,0);
    hist->GetYaxis()->SetTitle("yield");
    hist->GetXaxis()->SetTitle("time [ns]");
  }
}

double fitBack(TH1D *hist)
{
  double lowCut = 8000;
  double highCut = 15000;
  TF1 *gfit = NULL;
  double ave = 0;
  int lowBins = hist->FindBin(lowCut);
  int highBins = hist->FindBin(highCut);
  auto fitBack = new TF1("fitBack", "pol0", lowCut, highCut);
  fitBack->SetParameter(0, 1E-5);
  fitBack->SetParLimits(0, 1.E-12, 1.E2);

  hist->Fit("fitBack", "LF", " ", lowCut, highCut);
  gfit = (TF1 *)hist->GetListOfFunctions()->FindObject("fitBack");
  if (gfit)
  {
    ave = gfit->GetParameter(0);
  }
  else
    printf("P1 Fit to hist fails \n");
  double aveb = hist->Integral(lowCut, highCut) / double(highBins - lowBins);
  printf(" BBBB background fit ave %E aveb %E \n\n", ave, aveb);
  return ave;
}

int openFile(int fileNum=0)
{
  TString summaryFile[10];
  summaryFile[0] = TString("summary-type-1-dir-caenDataZeroPPM-2023-09-14-10-36.root"); // 05_22_2023 05_30_2023
  summaryFile[1] = TString("summary-type-1-dir-caenDataPointOnePPM-2023-09-14-13-28.root");
  summaryFile[2] = TString("summary-type-1-dir-caenDataPointTwoPPM-2023-09-14-13-28.root"); // 07_06_2023  07_07_2023
  summaryFile[3] = TString("summary-type-1-dir-caenDataPointFivePPM-2023-09-14-13-29.root");
  summaryFile[4] = TString("summary-type-1-dir-caenDataOnePPM-2023-09-15-11-11.root");
  summaryFile[5] = TString("summary-type-1-dir-caenDataTwoPPM-2023-09-15-13-14.root");
  dopant[0] = 0.05;
  dopant[1] = 0.1;
  dopant[2] = 1.0;
  dopant[3] = 0.50;
  dopant[4] = 1.0;
  dopant[5] = 2.0;
  TString fileName = summaryFile[fileNum];

  // open input filex
  vhist.resize(13);
  hffit.resize(13);
  hmodel.resize(13);
  hffitPmt.resize(13);
  hffitChan.resize(13);
  for (unsigned k = 0; k < vhist.size(); ++k)
    vhist[k] = NULL;
  printf(" looking for summary file %s\n", fileName.Data());

  bool exists = false;
  FILE *aFile;
  aFile = fopen(fileName.Data(), "r");
  if (aFile)
  {
    fclose(aFile);
    exists = true;
  }
  if (!exists)
  {
    printf(" couldnt open file %s\n", fileName.Data());
    return false;
  }
  // file exists, so open with TFile
  fin = new TFile(fileName, "readonly");
  printf(" opened file %s\n", fileName.Data());
  runSumDir = NULL;
  // get subdirectory pointer
  fin->GetObject("runSumDir", runSumDir);
  if (!runSumDir)
  {
    printf(" no runSumDir in file %s\n", fileName.Data());
    return false;
  }
  runSumDir->ls();
  cout << " making histogram list with " << histSet << "  " << runSumDir->GetName()  << endl;
  int iGot=0;
  for (int ic = 0; ic < 13; ++ic)
  {
    if (!goodChannel(ic))
      continue;
    TString hname;
    hname.Form("RunHitWaveChan%i", ic);
    TH1D *hWave = NULL;
    printf("\t look for  chan %i %s \n",ic,hname.Data());
    // get histogram by name
    runSumDir->GetObject(hname, hWave);
    if (hWave != NULL)
    {
      vhist[ic] = hWave;
      printf("\t got chan %i %s \n",ic,hWave->GetName());
      background[ic] = fitBack(hWave);
      ++iGot;
    }
  }
  return iGot;
}

void tbFitAll(int fileNum=0)
{

  TFile *fout = new TFile("tbFitAll.root", "recreate");
  // get data histograms from file
  int iGot = openFile(fileNum);
  printf(" got %i data histograms \n",iGot);
  if (iGot<1)
    return;

  for (unsigned ic = 0; ic < vhist.size(); ++ic)
    printf(" background fit ch %i value %E\n",ic,background[ic]);
  // fill buff

  // shift PMT by 46 ns = 23 bins
  int binShift = 23;
  for (int ibin = 0; ibin < vhist[12]->GetNbinsX() - binShift; ++ibin)
  {
    vhist[12]->SetBinContent(ibin,vhist[12]->GetBinContent(ibin+binShift)) ;
    vhist[12]->SetBinError(ibin,vhist[12]->GetBinError(ibin+binShift)) ;
  }


  /* fill buffer */
  for (unsigned ih = 0; ih < vhist.size(); ++ih)
    {
      if (vhist[ih] == NULL)
        continue;
      cout << ".... for channel %i " << ih << "  " << vhist[ih]->GetName() << " maximum bin " << vhist[ih]->GetMaximumBin() << endl;
      fout->Add(vhist[ih]);
      // fill data buffer
      for (int ib = 1; ib < vhist[ih]->GetNbinsX(); ++ib)
        buff[ih][ib] = vhist[ih]->GetBinContent(ib);
    }

  // clones to store fit histos
    for (unsigned ih = 0; ih < vhist.size(); ++ih){
      if(vhist[ih]==NULL)
        continue;
      vhist[ih]->GetListOfFunctions()->Clear();
      hffit[ih] = (TH1D *)vhist[ih]->Clone(Form("hFitCh%i", ih));
      fout->Add(hffit[ih]);
      hmodel[ih] = (TH1D *)vhist[ih]->Clone(Form("hModelCh%i", ih));
      fout->Add(hmodel[ih]);
    }

    // clones to store fit components histos 
    for (unsigned ih = 0; ih < NCOMP; ++ih)
    {
      TString histString;
      histString.Form("ModelPmt%s", compName[ih].Data());
      hffitPmt[ih] = (TH1D *)vhist[0]->Clone(histString);
      hffitPmt[ih]->SetTitle(histString);
      fout->Add(hffitPmt[ih]);
      hffitChan[ih] = (TH1D *)vhist[0]->Clone(Form("ModelCh7%s",compName[ih].Data() ));
      hffitChan[ih]->SetTitle(Form("ModelCh7%s",compName[ih].Data() ));
      fout->Add(hffitChan[ih]);
    }

      // Set starting values and step sizes for parameters
      static Double_t vstart[NPARS];
      static Double_t step[NPARS];
      /* added mix lifetime parameter 8
      >>> modelFit start parameters fit ppm 0.050000
        param 0 norm 2.0632E+03
        param 1 ppm 5.0000E-02
        param 2 tau3 1.1777E+03
        param 3 kp 1.0000E+00
        param 4 sfrac 3.0000E-01
        param 5 rfrac 1.0000E-03
        param 6 kxprime 1.0000E-01
        param 7 tmix 2.6647E+03
        param 8 taumix 4.7000E+03
        param 9  taurecomb    
      */

      // fit starting values
      vstart[0] = 2.06322e+03;
      vstart[1] = 0.05; //dopant[fileNum];
      vstart[2] = 1.1777E+03;
      vstart[3] = kqZero; // defined in  modelAFit.hh
      vstart[4] = 0.9;   // sfrac 0.886;
      vstart[5] = 1.E-3; // recomb  1.E-3;
      vstart[6] = 2.0*kqZero;
      vstart[7] = 4.7E+03; // was 2.66E-3
      vstart[8] = 4.7E3;
      vstart[9] = 1.44940e+01/4;

      TMinuit *gMinuit = new TMinuit(NPARS); // initialize TMinuit with a maximum of 5 params
      gMinuit->SetFCN(fcn);

      Double_t arglist[10];
      Int_t ierflg = 0;

      arglist[0] = 0.5; // for likelihood
      gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

      setParNames();
      for (unsigned j = 0; j < NPARS; ++j)
      {
        step[j] = 1.E-6 * vstart[j];
        gMinuit->mnparm(j, lparNames[j].Data(), vstart[j], step[j], 0, 0, ierflg);
        lpar[j] = vstart[j];
      }

      double fval = 0;
      double gin[NPARS];
      int npar = NPARS;
      fcn(npar, gin, fval, vstart, 0);
      printf(" >>>>   fval %E \n", fval);

      double currentValue;
      double currentError;

      
      /******************************/
      // fix most parameters
      /******************************/
      arglist[0] = 2; // par PPM
      //gMinuit->mnexcm("FIX", arglist, 1, ierflg);
      arglist[0] = 3; // par tau3
      gMinuit->mnexcm("FIX", arglist, 1, ierflg);
      // fix most parameters
      arglist[0] = 4; // kp
      //gMinuit->mnexcm("FIX", arglist, 1, ierflg);
      // fix most parameters
      arglist[0] = 5; // sfrac
      gMinuit->mnexcm("FIX", arglist, 1, ierflg);
      arglist[0] = 6; // recomb
      gMinuit->mnexcm("FIX", arglist, 1, ierflg);
      arglist[0] = 7; // kxprime
      //gMinuit->mnexcm("FIX", arglist, 1, ierflg);
      // fix most parameters
      arglist[0] = 8; // tmix
      gMinuit->mnexcm("FIX", arglist, 1, ierflg);
      arglist[0] = 9; // taumix
      gMinuit->mnexcm("FIX", arglist, 1, ierflg);
      arglist[0] = 10; //   taurecomb    
      //gMinuit->mnexcm("FIX", arglist, 1, ierflg);

      /*
          set limits ... here par starts with 1 so add 1
       */
      arglist[0] = 1;      // par norm
      arglist[1] = 1.0E2;  // low
      arglist[2] = 1.0E10; // high
      gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

      // set limits ... here par starts with 1 so add 1
      

      // set limits ... here par starts with 1 so add 1
      arglist[0] = 2;   // par PPM
      arglist[1] = 0.;  // low
      arglist[2] = 1.; // high
      gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

      // set limits ... here par starts with 1 so add 1
      arglist[0] = 4;  // par kp
      arglist[1] = 0.; // low
      arglist[2] = 3.; // high
      gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

      // fix sfrac
      arglist[0] = 5; // par sfrac
      // gMinuit->mnexcm("FIX", arglist, 1, ierflg);
      arglist[1] = 0.1; // low
      arglist[2] = .9;  // high
      gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

      // rfrac
      gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);
      arglist[0] = 6;     // par rfrac
      arglist[1] = 0.;    // low
      arglist[2] = 5.E-3; // high
      gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

      // arglist[0] = 6; // par rfrac
      // gMinuit->mnexcm("FIX", arglist, 1, ierflg);

      arglist[0] = 7;      // kxprime
      arglist[1] = 0.;     // low
      arglist[2] = 30.E-1; // high
      gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

      arglist[0] = 9;     // taumix
      arglist[1] = 1.0E3; // low
      arglist[2] = 10.E3; // high
      gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

      arglist[0] = 10;     // taumix
      arglist[1] = 1.; // low
      arglist[2] = 30; // high
      gMinuit->mnexcm("SET LIM", arglist, 3, ierflg);

      printf(" \n\n >>> modelFit start parameters fit ppm %f \n", vstart[1]);
      for (int ii = 0; ii < NPARS; ++ii)
      {
        gMinuit->GetParameter(ii, currentValue, currentError);
        printf("\t  param %i %s %.4E  \n", ii, lparNames[ii].Data(), currentValue);
      }

      double amin;
      gMinuit->mnprin(1, amin);

      // gMinuit->FixParameter(6);
      // create functions to plot
      int myColor[13] = {kRed, kBlue - 9, kGreen, kOrange, kBlue, kAzure + 3, kCyan - 3, kGreen - 4, kGreen, kSpring, kYellow, kBlue + 9, kTeal - 4};
      int myStyle[13] = {21, 22, 23, 24, 25, 26, 21, 22, 23, 31, 32, 33, 34};
      TString cname;

      for (int ifit = 0; ifit < 13; ++ifit)
      {
        if (vhist[ifit] == NULL)
          continue;
        vhist[ifit]->GetListOfFunctions()->Clear();
        vhist[ifit]->SetLineColor(myColor[ifit]);
        vhist[ifit]->SetMarkerColor(myColor[ifit]);
        vhist[ifit]->SetMarkerStyle(myStyle[ifit]);
        // if (histSet == TString("HitWave"))
        // vhist[ifit]->GetYaxis()->SetRangeUser(1.E-7, 1.);
        // else
        // vhist[ifit]->GetYaxis()->SetRangeUser(1.E-3, 1.E3);
      }

      gStyle->SetOptFit(1111);
      gStyle->SetOptStat(11);

      cname.Form("DataDopant-%.f", lpar[1]);
      TCanvas *cand = new TCanvas(cname, cname);
      cand->SetLogy();
      for (int k = 0; k < 13; ++k)
      {
        if (!goodChannel(k))
          continue;
        if (k == 0)
          vhist[k]->Draw();
        else
          vhist[k]->Draw("sames");
      }
      cand->BuildLegend();
      cand->Print(".png");
      createFunctions();
      // show();

      for (int k = 0; k < 13; ++k)
      {
        fillHistogram(ffit[k], hffit[k]);
      }

      cname.Form("FitChan-Dopant-%.3f", lpar[1]);
      TCanvas *canAllChan = new TCanvas(cname, cname);
      canAllChan->SetLogy();
      for (int k = 0; k < 13; ++k)
      {
        if (!goodChannel(k))
          continue;
        hffit[k]->SetLineColor(myColor[k]);
        hffit[k]->GetYaxis()->SetRangeUser(1.E-10, 1.);
        if (k == 0)
          hffit[k]->Draw();
        else
          hffit[k]->Draw("sames");
      }
      canAllChan->BuildLegend();
      canAllChan->Print(".png");

      for (int k = 0; k < NCOMP-1; ++k)
      {
        fillHistogram(ffitPmt[k], hffitPmt[k]);
      }
      // fill last hist comp with background

      for (int ibin = 0; ibin < hffitPmt[NCOMP - 1]->GetNbinsX(); ++ibin){
        hffitPmt[NCOMP - 1]->SetBinContent(ibin, background[12]);
        hffitPmt[NCOMP - 1]->SetBinError(ibin,0);
      }

      // PMT components
      cname.Form("ModelPmt-Dopant-%.3f-recombination-%.2E", lpar[1], vstart[5]);
      TCanvas *canPmt = new TCanvas(cname, cname);
      canPmt->SetLogy();
      for (int k = 0; k < NCOMP; ++k)
      {

        TString histString;
        histString.Form("ModelPmt%s", compName[k].Data());
        hffitPmt[k]->SetName(histString);
        hffitPmt[k]->SetTitle(histString);
        hffitPmt[k]->SetLineColor(myColor[k]);
        hffitPmt[k]->GetYaxis()->SetRangeUser(1.E-7, 1.);
        fout->Add(hffitPmt[k]);
        cout << hffitPmt[k]->GetName() << " " << hffitPmt[k]->GetTitle() << " name " << compName[k] << endl;
        if (k == 0)
          hffitPmt[k]->Draw("C");
        else
          hffitPmt[k]->Draw("Csames");
      }
      //vhist[12]->Draw("same");
      // hffit[12]->Draw("same");
      canPmt->BuildLegend();
      canPmt->Print(".png");

      //return;

      for (int k = 0; k < 13; ++k)
      {
        fillHistogram(ffitChan[k], hffitChan[k]);
      }

      // channel7  components
      cname.Form("ModelChan7Dopant-%.3f", lpar[1]);
      TCanvas *canChan = new TCanvas(cname, cname);
      canChan->SetLogy();
      for (int k = 0; k < NCOMP-1; ++k)
      {
        // ffitChan[k]->GetYaxis()->SetRangeUser(1.E-10, 3E-1);
        TString histString;
        histString.Form("hModelChan7%s", compName[k].Data());
        hffitChan[k]->SetName(histString);
        hffitChan[k]->SetTitle(histString);
        hffitChan[k]->Print();
        hffitChan[k]->SetLineColor(myColor[k]);
        hffitChan[k]->GetYaxis()->SetRangeUser(1.E-5, 1.);
        fout->Add(hffitChan[k]);
        printf(" comp %i\n", int(ffitChan[k]->GetParameter(1)));
        if (k == 0)
          hffitChan[k]->Draw("C");
        else
          hffitChan[k]->Draw("Csames");
      }
      canChan->BuildLegend();
      canChan->Print(".png");

      /* if just plotting model return here*/
      // return;

      // minimize with MIGRAD
      // Now ready for minimization step
      arglist[0] = 100000; // maxcalls
      arglist[1] = 1.E-2;  // tolerance

      /* MIGrad[maxcalls][tolerance]*/
      gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

      // Print results
      Double_t edm, errdef;
      Int_t nvpar, nparx, icstat;
      gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
      /*
      Prints the values of the parameters at the time of the call.
      According to the value of IKODE, the printout is: IKODE=INKODE=

    0 only info about function value
    1 parameter values, errors, limits
    2 values, errors, step sizes, internal values
    3 values, errors, step sizes, first derivs.
    4 values, parabolic errors, MINOS errors
    when INKODE=5, MNPRIN chooses IKODE=1,2, or 3, according to fISW[1]
      */
      gMinuit->mnprin(1, amin);

      // put parameters into model par array
      for (int k = 0; k < NPARS; ++k)
      {
        gMinuit->GetParameter(k, currentValue, currentError);
        lpar[k] = currentValue;
        printf("\t copy %s new value %f \n", lparNames[k].Data(), lpar[k]);
      }

      for (int k = 0; k < 13; ++k){
        fillHistogram(fmodel[k], hmodel[k]);
      }

      cname.Form("ModelAllChannelsDopant-%.f", lpar[1]);
      TCanvas *canf = new TCanvas(cname, cname);
      canf->SetLogy();
      for (int k = 0; k < 13; ++k)
      {
        if (!goodChannel(k))
          continue;
        hmodel[k]->SetLineColor(myColor[k]);
        hmodel[k]->GetYaxis()->SetRangeUser(1.E-4, 1.);
        if (k == 0)
          hmodel[k]->Draw("C");
        else
          hmodel[k]->Draw("csames");
      }
      canf->BuildLegend();
      canf->Print(".png");

      /* fill hffit histograms */
      for (int k = 0; k < 13; ++k)
        fillHistogram(ffit[k], hffit[k]);

      for (int k = 0; k < NCOMP-1; ++k)
        fillHistogram(ffitPmt[k], hffitPmt[k]);

      for (int k = 0; k < 13; ++k)
      {
        if (!goodChannel(k))
          continue;
        cname.Form("FitToChan%i-Dopant-%.3f-recombination-%.2E", k, lpar[1], vstart[5]);
        TCanvas *can = new TCanvas(cname, cname);
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(11);
        can->SetLogy();
        hffit[k]->SetLineColor(kBlack);
        cname.Form("DataChannel%i", k);
        vhist[k]->SetTitle(cname);
        hffit[k]->SetLineWidth(1);
        vhist[k]->SetMarkerSize(0.4);
        vhist[k]->GetYaxis()->SetTitle("yield");
        vhist[k]->GetXaxis()->SetTitle("time [ns]");
        vhist[k]->GetXaxis()->SetRangeUser(xlow, xhigh);
        vhist[k]->GetYaxis()->SetRangeUser(1E-5, 3E-1);
        hffit[k]->GetXaxis()->SetRangeUser(xlow, xhigh);
        hffit[k]->GetYaxis()->SetRangeUser(1E-5,3E-1);
        hffit[k]->Draw("");
        if(k==12) {
          for (int icomp =1 ; icomp < NCOMP; ++icomp)
            hffitPmt[icomp]->Draw("same");
        }
        vhist[k]->Draw("same");
        fout->Add(ffit[k]);
        can->BuildLegend();
        can->Print(".png");
      }
      fout->Add(ntScan);
      fout->Write();

      for (int k = 0; k < 13; ++k)
        printf("channel %i background fit %E \n",k,background[k]);
    /*
      // PMT components
      cname.Form("FitPmt-Dopant-%.3f-recombination-%.2E", lpar[1], vstart[5]);
      TCanvas *canPmtFit = new TCanvas(cname, cname);
      canPmt->SetLogy();
      for (int k = 0; k < NCOMP; ++k)
      {

        TString histString;
        histString.Form("hModelPmtComp%s", compName[k].Data());
        hffitPmt[k]->SetName(histString);
        hffitPmt[k]->SetTitle(histString);
        hffitPmt[k]->SetLineColor(myColor[k]);
        hffitPmt[k]->GetYaxis()->SetRangeUser(1.E-7, 1.);
        fout->Add(hffitPmt[k]);
        cout << hffitPmt[k]->GetName() << " " << hffitPmt[k]->GetTitle() << " name " << compName[k] << endl;
        if (k == 0)
          hffitPmt[k]->Draw("C");
        else
          hffitPmt[k]->Draw("Csames");
      }
      vhist[12]->Draw("same");
      // hffit[12]->Draw("same");
      canPmtFit->BuildLegend();
      canPmtFit->Print(".png");
      */
}
