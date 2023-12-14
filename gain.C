std::vector<double> fPositionX;
std::vector<double> fPositionY;
std::vector<double> fFitX;
std::vector<double> fFitXSigma;
std::vector<double> fFitY;
std::vector<double> fFitYError;

std::vector<TH1D *> hlist;
enum
{
  nbins = 4000
};

void gain()
{
  TH1D *hFit = new TH1D("GausFit","GausFit",4000, -20.E3, 200.E3);
  TFile *fin = new TFile("post-11_26_2023-553878.root", "readonly");
  TFile *fout = new TFile("gain", "recreate");
  TString histname;
  for (int ichan = 0; ichan < 12; ++ichan)
  {
    histname.Form("TotSumChan%i", ichan);
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
    // int nfound = s->Search(hlist[i],4,"",5);
    for (int j = 0; j < nbins; ++j)
    {
      if (hlist[i]->GetBinLowEdge(j + 1) < 6000.)
        source[j] = 0;
      else
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
    fFitX.clear();
    fFitY.clear();
    fFitXSigma.clear();
    fFitYError.clear();
    double ymax = 0;
    int nfound = s->SearchHighRes(source, dest, nbins, 8, 3, kTRUE, 2, kTRUE, 2);
    printf("peaks for hist %s found %i \n", hlist[i]->GetName(), nfound);
    Double_t *xpeaks = s->GetPositionX();
    for (int k = 0; k < nfound; k++)
    {
      int bin = 1 + Int_t(xpeaks[k] + 0.5);
      double x = hlist[i]->GetBinCenter(bin);
      double y = hlist[i]->GetBinContent(bin);
      // cut close peaks
      if (fPositionX.size() > 0)
      {
        double diff = abs(x - fPositionX[fPositionX.size() - 1]);
        printf("  check  ..  %i %f %f diff %f \n", k, x, fPositionX[fPositionX.size() - 1], diff);
        if (diff < 10.E3 || y < 10 || y / ymax < 1.E-3)
        {
          printf("  skip ..  %i %f %f \n", k, x, fPositionX[fPositionX.size() - 1]);
          continue;
        }
      }
      if (y > ymax)
        ymax = y;
      //
      fPositionX.push_back(x);
      fPositionY.push_back(y);
      //printf(" add peak %lu at x = %f val %f \n", fPositionX.size(), fPositionX[fPositionX.size() - 1], fPositionY[fPositionY.size() - 1]);
    }

    

    // fit
    fFitX.resize(fPositionX.size());
    fFitY.resize(fPositionX.size());
    fFitXSigma.resize(fPositionX.size());
    fFitYError.resize(fPositionX.size());
    printf(" peaks channel %i n= %lu \n", i, fPositionX.size());
    for (unsigned long j = 0; j < fPositionX.size(); ++j)
    {

      // fit each peak to gaussian
      //TH1D *hFit = (TH1D *)hlist[i]->Clone(Form("gausFitPeak%luHist%i", j, i));
      
      hFit->Reset("ICES");
      for (int ibin = 0; ibin < hlist[i]->GetNbinsX(); ++ibin)
        hFit->SetBinContent(ibin, hlist[i]->GetBinContent(ibin));
      // TH1D *hFit = (TH1D *)hlist[i]->Clone(Form("gausFitPeak%luHist%i", j, i));
      hFit->Fit("gaus", "QO", "", fPositionX[j], fPositionX[j] + 1000);
      TF1 *gfit = (TF1 *)hFit->GetListOfFunctions()->FindObject("gaus");
      double mean = 0;
      double sigma = 0;
      double val = 0;
      if (gfit != nullptr)
      {
        mean = gfit->GetParameter(1);
        sigma = gfit->GetParameter(2);
        int bin = hFit->FindBin(mean);
        val = hFit->GetBinContent(bin);
        fFitX[j]=mean;
        fFitXSigma[j]=sigma;
        fFitY[j]=val;
        fFitYError[j]=sqrt(val);
      } else {
        printf("\n\n!!!!!fit fails hist %i point %lu \n\n\n",i,j);
      }
    }

      printf(" \n \n fit for channel %i n= %lu \n", i, fFitX.size());
      for (unsigned long j = 0; j < fFitX.size(); ++j)
        printf(" peak %lu x (%f,%f) sigma %f y (%f,%f) \n", j, fPositionX[j],fFitX[j], fFitXSigma[j], fPositionY[j],fFitY[j]);

      // make graph
      auto g = new TGraphErrors(fFitX.size(), &fFitX[0], &fFitY[0], &fFitXSigma[0], &fFitYError[0]);
      g->SetMarkerStyle(23);
      g->SetMarkerColor(kRed);
      g->SetMarkerSize(1.3);

      

      TString gname;
      gname.Form("gainChan%i", i);
      TString gtitle;
      gtitle.Form(" gain channel  %i ;  ADC ; gain", i);
      g->SetName(gname);
      g->SetTitle(gtitle.Data());
      g->Fit("pol0", "");
      g->GetHistogram()->GetListOfFunctions()->ls();
      TF1 *gfit = g->GetFunction("pol0");
      TCanvas *gcan = new TCanvas(Form("GainChan%i", i), Form("chan%i", i));
      gPad->SetLogy(0);
      gStyle->SetOptFit();
      g->Draw("AP");
      gfit->SetLineStyle(5);
      gfit->SetLineWidth(5);
      gfit->Draw("same");
      gcan->Print(".png");
      fout->Add(g);

      // add marker
      TPolyMarker *pm = (TPolyMarker *)hlist[i]->GetListOfFunctions()->FindObject("TPolyMarker");
      if (pm)
      {
        hlist[i]->GetListOfFunctions()->Remove(pm);
        delete pm;
      }

      // pm = new TPolyMarker(fFitX.size(), &fFitX[0], &fFitY[0]);
      pm = new TPolyMarker(fFitX.size(), &fFitX[0], &fFitY[0]);
      hlist[i]->GetListOfFunctions()->Add(pm);
      pm->SetMarkerStyle(23);
      pm->SetMarkerColor(kRed);
      pm->SetMarkerSize(1.3);

      // plot
      TCanvas *can = new TCanvas(Form("PeaksChan%i", i), Form("chan%i", i));
      gPad->SetLogy(1);
      hlist[i]->Draw();
      can->Print(".png");
      fout->Add(hlist[i]);
    }

    fout->Write();
}
