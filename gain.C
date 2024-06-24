/* read and fit gains from summary June 24 2024 */
#include <ctime>
#include "TDirectory.h"
double nominalGain = 227.4; // average
std::vector< std::vector<TF1*> > vfit;
std::vector<double> fPositionX;
std::vector<double> fPositionY;
std::vector<double> fFitADCY;
std::vector<int> fFitBin;
std::vector<double> fFitADC;
std::vector<double> fFitADCError;
std::vector<double> fSpeNumber;
std::vector<double> fSpeNumberError;
std::vector<double> sipmGain;
std::vector<double> sipmGainError;
std::vector<double> sipmNumber;
std::vector<double> sipmNumberError;

enum {NPEAKS=4};

std::vector<TH1D *> hlist;

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
  //TH1D *hFit = new TH1D("GausFit", "GausFit", 4000, -20.E3, 200.E3);
  TFile *fin = new TFile("summary-12_28_2023-nfiles-2-dir-caenData-2024-06-24-13-18.root", "readonly");
  std::string sdate = currentDate();
  TFile *fout = new TFile(Form("gains-%s.root", sdate.c_str()), "recreate");
  TF1 *line = new TF1("myLine", fline, 0, 2.E5, 2);
  // get histos from file
  TDirectory *gainSumDir=nullptr;
  fin->GetObject("gainSumDir",gainSumDir);

 //
  TIter next(gainSumDir->GetListOfKeys());
  TKey *key;
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());

    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    TString hname(h->GetName());
    int ichan = TString(hname(hname.Last('n') + 1, hname.Length())).Atoi();
    // struggle to get only last hist cycle
    TString lastName;
    if(hlist.size()>0) lastName=hlist[hlist.size()-1]->GetName();
    //cout << " xxx " << h->GetName() << " " << lastName << endl;
    if(TString(h->GetName())==lastName) continue;
    if(ichan>11) continue;
    cout << "hist name " << h->GetName() << " chan " << ichan << " cycle " << key->GetCycle() << endl;
    hlist.push_back(h);
  }
   printf(" have %lu gain hists \n",hlist.size());

   // fit peaks 
   vfit.resize(hlist.size());
   double width = 100.;
   unsigned ihist=4;
   for(unsigned ihist=0; ihist<hlist.size(); ++ihist) {
   for (unsigned j = 0; j < NPEAKS; ++j)
    {
      TH1D* h=hlist[ihist];
      double nominalPeak =nominalGain*double(j+1);
      double fitStart = nominalGain*double(j+1) -width;
      double fitEnd   = nominalGain*double(j+1) + width;
      int bin =h->FindBin(nominalPeak);
      double xbin = h->GetBinLowEdge(bin);
      double peakVal=h->GetBinContent(bin);
      if(peakVal<3) continue;
      // fit here
      h->GetListOfFunctions()->Clear();
      //printf("fit to hist %s peak # %u bin %i xbin %f val %f from %f to %f \n",
      //    h->GetName(), j+1, bin,xbin,peakVal,fitStart,fitEnd);
      h->Fit("gaus","QQSS","", fitStart,fitEnd);
      TF1 *gFit = (TF1 *)h->GetListOfFunctions()->FindObject("gaus");
      if (gFit != nullptr)
      {
        vfit[ihist].push_back(gFit);
        double mean = gFit->GetParameter(1);
        double meanError = gFit->GetParError(1);
        int jbin = h->FindBin(mean);
        double val = h->GetBinContent(jbin);
        printf("fit to hist %s point %u  nominal bin %i nominal x %.2f from %.0f to %.0f  mean %.2f error %.2f peak bin %i val %.2f  \n", 
          h->GetName(),j+1,bin, xbin, fitStart,fitEnd, mean,meanError, jbin, val);
      }
    }
   }

   for(unsigned ihist=0; ihist< hlist.size(); ++ihist) printf(" good fits to %s = %lu \n",
    hlist[ihist]->GetName(),vfit[ihist].size());

    // plot
    for(unsigned ihist=0; ihist< hlist.size(); ++ihist) {
      hlist[ihist]->GetListOfFunctions()->Clear();
      TString cname;
      cname.Form("Gain%s",hlist[ihist]->GetName());
      TCanvas *c = new TCanvas(cname,cname);
      hlist[ihist]->Draw();
      for(unsigned ifit=0; ifit< vfit[ihist].size(); ++ifit) vfit[ihist][ifit]->Draw("same"); 
    }
}