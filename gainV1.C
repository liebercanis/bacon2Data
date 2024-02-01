#include <ctime>
std::vector<double> fPositionX;
std::vector<double> fPositionY;
std::vector<double> fFitADC;
std::vector<double> fFitADCY;
std::vector<double> fFitADCError;
std::vector<double> fSpeNumber;
std::vector<double> fSpeNumberError;
std::vector<double> sipmGain;
std::vector<double> sipmGainError;
std::vector<double> sipmNumber;
std::vector<double> sipmNumberError;

std::vector<TH1D *> hlist;
enum
{
  nbins = 4000
};

std::string currentDate()
{
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  char output[30];
  strftime(output, 30, "%Y-%m-%d-%H-%M", timeinfo);
  return std::string(output);
}
// Define a linear fit
double fline(double *x, double *par)
{
  return par[0] + x[0] * par[1];
}

void gain()
{
  TH1D *hFit = new TH1D("GausFit", "GausFit", 4000, -20.E3, 200.E3);
  TFile *fin = new TFile("post-12_28_2023-10127754.root", "readonly");
  std::string sdate = currentDate();
  TFile *fout = new TFile(Form("gains-%s.root", sdate.c_str()), "recreate");
  TF1 *line = new TF1("myLine", fline, 0, 2.E5, 2);
  TString histname;
  for (int ichan = 0; ichan < 12; ++ichan)
  {
    histname.Form("LatePeakSumChan%i", ichan);
    TH1D *hist;
    fin->GetObject(histname, hist);
    cout << ichan << "  " << hist->GetName() << endl;
    hlist.push_back(hist);
  }
  cout << " got " << hlist.size() << endl;
  TSpectrum *s = new TSpectrum();

  /*
  for(unsigned i=0; i<hlist.size();++i) {
    TH1D* hpeak = (TH1D*) hlist[i]->Clone(Form("PeaksChan%i",i));
    for(int j=0; j<nbins; ++j) {
      if( hpeak->GetBinLowEdge(j+1) < 5000. )
        hpeak->SetBinContent(j+1,0);
      else
        hpeak->SetBinContent(j+1,hlist[i]->GetBinContent(j+1));
    }
    */

  double source[nbins];
  double dest[nbins];
  for (unsigned i = 0; i < hlist.size(); ++i)
  {
    if (i == 5)
      continue;
    // int nfound = s->Search(hlist[i],4,"",5);
    for (int j = 0; j < nbins; ++j)
    {
      source[j] = hlist[i]->GetBinContent(j + 1);
    }
    /*
  source: pointer to the vector of source spectrum.
  destVector: pointer to the vector of resulting deconvolved spectrum.
  ssize: length of source spectrum.
  sigma: sigma of searched peaks, for details we refer to manual.
  threshold: threshold value in % for selected peaks, peaks with amplitude less than threshold*highest_peak/100 are ignored, see manual.
  backgroundRemove: logical variable, set if the removal of background before deconvolution is desired.
  deconIterations-number of iterations in deconvolution operation.
  markov: logical variable, if it is true, first the source spectrum is replaced by new spectrum calculated using Markov chains method.
  averWindow: averaging window of searched peaks, for details we refer to manual (applies only for Markov method).
  */

    fPositionX.clear();
    fPositionY.clear();
    fFitADC.clear();
    fFitADCY.clear();
    fFitADCError.clear();
    fSpeNumber.clear();
    fSpeNumberError.clear();
    double ymax = 0;
    /*
    Parameters:
source: pointer to the vector of source spectrum.
destVector: pointer to the vector of resulting deconvolved spectrum.
ssize: length of source spectrum.
sigma: sigma of searched peaks, for details we refer to manual.
threshold: threshold value in % for selected peaks, peaks with amplitude less than threshold*highest_peak/100 are ignored, see manual.
backgroundRemove: logical variable, set if the removal of background before deconvolution is desired.
deconIterations-number of iterations in deconvolution operation.
markov: logical variable, if it is true, first the source spectrum is replaced by new spectrum calculated using Markov chains method.
averWindow: averaging window of searched peaks, for details we refer to manual (applies only for Markov method).
    */
    int nfound = s->SearchHighRes(source, dest, nbins, 8, .1, kTRUE, 2, kTRUE, 2);
    printf("\n\npeaks for hist %s found %i \n", hlist[i]->GetName(), nfound);
    Double_t *xpeaks = s->GetPositionX();
    for (int k = 0; k < nfound; k++)
    {
      int bin = 1 + Int_t(xpeaks[k] + 0.5);
      double x = hlist[i]->GetBinCenter(bin);
      double y = hlist[i]->GetBinContent(bin);
      // cut close peaks
      /*
      if (fPositionX.size() > 0)
      {
        double diff = abs(x - fPositionX[fPositionX.size() - 1]);
        //printf("  check  ..  %i %f %f diff %f \n", k, x, fPositionX[fPositionX.size() - 1], diff);
        if (diff < 10.E3 || y / ymax < 1.E-4  || y<1.)
        {
          //printf("  skip ..  %i %f %f \n", k, x, fPositionX[fPositionX.size() - 1]);
          continue;
        }
        printf("  skip ..  %i %f %f \n", k, x, fPositionX[fPositionX.size() - 1]);
      }
      */
      if (y > ymax)
        ymax = y;

      //printf("\t det %i peak %i  height %f x %f \n", i, k, y, x);
      //

      if(x<150)
        continue;
      if(y<10)
        continue;
      
      if (fPositionX.size()>0){
        if (x - fPositionX[fPositionX.size()-1] < 100)
          continue;
      }

      fPositionX.push_back(x);
      fPositionY.push_back(y);
      double previous = 0;
      if (fPositionX.size() > 1)
        previous = fPositionX[fPositionX.size() - 2];

      printf(" add peak #  %lu at x = %f previous x %f y= %f \n", fPositionX.size(), x, previous,y);

      if (fPositionX.size() > 4)
        break;
    }


    // fit
    fFitADC.resize(fPositionX.size());
    fFitADCY.resize(fPositionX.size());
    fFitADCError.resize(fPositionX.size());
    fSpeNumber.resize(fPositionX.size());
    fSpeNumberError.resize(fPositionX.size());
    printf(" peaks channel %i n= %lu \n", i, fPositionX.size());
    for (unsigned long j = 0; j < fPositionX.size(); ++j)
    {
      fSpeNumberError[j] = double(0);
      fSpeNumber[j] = double(j);

      // skip fit each peak to gaussian
      // TH1D *hFit = (TH1D *)hlist[i]->Clone(Form("gausFitPeak%luHist%i", j, i));
      if(0){
      hFit->Reset("ICES");
      printf("\t fit to hist %i point %lu  x %f \n", i, j, fPositionX[j]);
      for (int ibin = 0; ibin < hlist[i]->GetNbinsX(); ++ibin)
        hFit->SetBinContent(ibin, hlist[i]->GetBinContent(ibin));
      // TH1D *hFit = (TH1D *)hlist[i]->Clone(Form("gausFitPeak%luHist%i", j, i));
      hFit->Fit("gaus","","", fPositionX[j]-100, fPositionX[j] + 100);
      TF1 *gFit = (TF1 *)hFit->GetListOfFunctions()->FindObject("gaus");
      double mean = 0;
      double meanError = 0;
      double val = 0;
      if (gFit != nullptr)
      {
        mean = gFit->GetParameter(1);
        meanError = gFit->GetParError(1);
        int bin = hFit->FindBin(mean);
        val = hFit->GetBinContent(bin);
        fFitADC[j] = mean;
        fFitADCY[j] = val;
        fFitADCError[j] = meanError;
        printf("\t fit to hist %i point %lu  (%f,%f) bin %i x %f xadc %f  \n", i, j,fPositionX[j],fPositionY[j],bin,
        fFitADC[j],fFitADCY[j]);
      }
      else
      {
        printf("\n\n!!!!!fit fails hist %i point %lu \n\n\n", i, j);
      }
    } else {
      fFitADC[j] = fPositionX[j];
      fFitADCY[j] = fPositionY[j];
      fFitADCError[j] =0;
    }
    } // loop over points

    printf(" fit for channel %i n= %lu \n", i, fFitADC.size());
    for (unsigned long j = 0; j < fFitADC.size(); ++j)
      printf("  point %lu Position %.2f ADC %.2f +/- %.2f  y %.2f\n", 
        j, fPositionX[j], fFitADC[j], fFitADCError[j],fFitADCY[j]);

    if (fFitADC.size()<2)
      continue;
    TGraphErrors *g = new TGraphErrors(fFitADC.size(), &fSpeNumber[0], &fFitADC[0], &fSpeNumberError[0], &fFitADCError[0]);
    //g->Print();
    g->SetMarkerStyle(23);
    g->SetMarkerColor(kRed);
    g->SetMarkerSize(1.3);

    // add marker
    TPolyMarker *pmold = (TPolyMarker *)hlist[i]->GetListOfFunctions()->FindObject("TPolyMarker");
    if (pmold)
    {
      hlist[i]->GetListOfFunctions()->Remove(pmold);
      delete pmold;
    }
    TPolyMarker *pm = new TPolyMarker(fFitADC.size(), &fFitADC[0], &fFitADCY[0]);
    hlist[i]->GetListOfFunctions()->Add(pm);
    double *xp = pm->GetX();
    double *yp = pm->GetY();
    for (int ip = 0; ip < pm->GetN(); ++ip)
      printf("det %i poly %i %f %f \n", i, ip, xp[ip], yp[ip]);
    pm->SetMarkerStyle(23);
    pm->SetMarkerColor(kRed);
    pm->SetMarkerSize(1.5);

    // make graph

    TString gname;
    gname.Form("gainChan%i", i);
    TString gtitle;
    gtitle.Form(" gain channel  %i ;  number SPE ; ADC/SPE", i);
    g->SetName(gname);
    g->SetTitle(gtitle.Data());
    line->SetParameters(0.5, 0);
    g->Fit("myLine");
    // g->GetHistogram()->GetListOfFunctions()->ls();
    TF1 *gFit = g->GetFunction("myLine");
    if (gFit == nullptr)
      continue;
    TCanvas *gcan = new TCanvas(Form("GainChan%i", i), Form("chan%i", i));
    gPad->SetLogy(0);
    gStyle->SetOptFit();
    g->GetHistogram()->GetXaxis()->SetRangeUser(0, 4);
    g->GetHistogram()->GetYaxis()->SetRangeUser(0, 1.5E5);
    gFit->SetLineStyle(5);
    gFit->SetLineWidth(5);
    g->Draw("APE1");
    gFit->Draw("same");
    gPad->SetGrid();
    gcan->Print(".png");
    fout->Add(g);
    if (gFit->GetParameter(1)<0)
      continue;
    sipmGain.push_back(gFit->GetParameter(1));
    sipmGainError.push_back(gFit->GetParError(1));
    sipmNumber.push_back(i);
    sipmNumberError.push_back(0);

    

    // plot
    TCanvas *can = new TCanvas(Form("PeaksChan%i", i), Form("chan%i", i));
    gPad->SetLogy(1);
    hlist[i]->Draw();
    fout->Add(hlist[i]);
    can->Print(".png");
    //fout->Add(can);
  } // channel number loop
  // plot
  TGraphErrors *gGain = new TGraphErrors(sipmGain.size(), &sipmNumber[0], &sipmGain[0], &sipmNumberError[0], &sipmGainError[0]);
  gGain->SetName("gGain");
  gGain->SetTitle("sipm gain");
  gGain->SetMarkerStyle(23);
  gGain->SetMarkerColor(kRed);
  gGain->SetMarkerSize(1.3);

  TCanvas *canGain = new TCanvas("sipmGain", "sipm gain");
  gPad->SetGrid();
  gGain->Draw("AP");
  canGain->Print(".png");
  fout->Add(gGain);
  fout->Write();
}
