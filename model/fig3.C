// times in ns
static const double lar3=126.9;
static const double FWHMar3=2.*3.6;
static const double war3 = 2.45; //FWHM/2.355
//
static const double ltrap=126.5;
static const double l1 = 11.6;



using namespace TMath;

TF1* flor;
TF1* fatt;
TF1* far;
TF1* fcon;

std::vector<double> vwave;
std::vector<double> vamp;
std::vector<double> vline;
std::vector<double> wline;

std::vector<double> xvec;
std::vector<double> avec;


static double Lorentz(Double_t *x, Double_t *par) // max, mean
{
  double ghalf = 1./Pi()/par[0];
  double x0 = par[1];
  double l = 1./Pi()*ghalf/( pow(x[0]-x0,2.) +  pow(ghalf,2.));
  return l;
}

static double MLorentz(Double_t *x, Double_t *par) // max, mean
{
  double ghalf = 1./Pi()/par[0];
  double x0 = par[1];
  double l = 1./Pi()*ghalf/( pow(x[0]-x0,2.) +  pow(ghalf,2.));
  return 1.-l;
}

unsigned long getSpectrum( TString fileName ){
  ifstream* in = new ifstream(fileName,std::ios::in);
  char line[1024];

  while( !(in->eof())  ) {
    in->getline(line,1024);
    TString tline(line);
    TObjArray* tokenArray = tline.Tokenize(',');
    if(tokenArray->GetEntries()<2) break;
    TObjString *sv = (TObjString*) tokenArray->At(0);
    if(sv->GetString().Contains("wavelength")) {
      cout << " comment " << tline << endl;
      continue;
    } else {
      TObjString *sw = (TObjString*) tokenArray->At(0);
      TObjString *sa = (TObjString*) tokenArray->At(1);
      double a = sa->GetString().Atof();
      double w = sw->GetString().Atof();
      if(a==0) continue;
      vwave.push_back(w);
      vamp.push_back(a);
      //cout << vwave.size()-1 << " w " << w << " a " << a << endl;
    }
  }
  return vwave.size();
}

void fig3()
{
  TString fileName1("0.1ppm_transmission.csv");
  TString fileName2("1ppm_2mm_trans.txt");
  TFile *fout = new TFile("fig3Out.root","recreate");
  TString canTitle;

  unsigned long nlines = getSpectrum(fileName2);
  printf(" number of points in spectrum is %lu \n",nlines);

  double wlow=120;
  double whigh=136;


  canTitle.Form("VUVEmission1PPM");
  TGraph *arSpect = new TGraph(vwave.size(),&vwave[0],&vamp[0]);
  TCanvas *canSpect = new TCanvas(canTitle,canTitle);
  gStyle->SetOptFit();
  arSpect->SetTitle("VUV Emission 1PPM");
  arSpect->GetYaxis()->SetTitle("amplitude");
  arSpect->GetXaxis()->SetTitle("wavelength (nm) ");
  arSpect->SetMarkerStyle(22);
  arSpect->SetMarkerSize(0.4);
  //arSpect->Fit("gaus");
  arSpect->Draw("alp");
  canSpect->Print(".png");
}
