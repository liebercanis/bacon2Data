// times in ns
#include <fstream>
#include <iostream>

static const double lar3=126.9;
static const double FWHMar3=2.*3.6;
static const double war3 = 2.45; //FWHM/2.355
//
static const double ltrap=126.5;
static const double l1 = 11.6;
ofstream textFile;


using namespace TMath;

TF1* flor;
TF1* fatt;
TF1* far;
TF1* fcon;
TF1* fishida;

std::vector<double> vwave;
std::vector<double> vamp;
std::vector<double> vtrans;
std::vector<double> xline;
std::vector<double> wline;

std::vector<double> logIshida;
std::vector<double> vishida;
std::vector<double> wishida;


std::vector<double> xvec;
std::vector<double> avec;
std::vector<double> vabs;

TGraph *arTrans;
TGraph *IshidaAtten;
TH1D* hTrans;

void writeGraph(TGraph *gr)
{
    textFile  << gr->GetTitle() << "   x=wavelength [nm], y = transmission "  << endl;

    double x,y;

    for (int i = 0; i < gr->GetN(); ++i){
        gr->GetPoint(i,x,y);
        textFile << i << "  " << x << "  " << y << endl;
    }
}



static double TwoExp(Double_t *xx, Double_t *par) // max, mean
{
  double x = xx[0];
  double c = par[0];
  double l1 = par[1];
  double l2 = par[2];
  double f=  c*Exp(-x/l1) + (1-c)*Exp(-x/l2);
  return f;
}


static double Lorentz(Double_t *x, Double_t *par) // max, mean
{
  double ghalf = 1./Pi()/par[0];
  double x0 = par[1];
  double l = 1./Pi()*ghalf/( pow(x[0]-x0,2.) +  pow(ghalf,2.));
  return l;
}


static double fitAtten(Double_t *xx, Double_t *par) // max, mean
{
  double x = xx[0];
  double tau = par[0];
  double A = par[1];
  return A*(1.- Exp(-x/tau));
}





static double MLorentz(Double_t *x, Double_t *par) // max, mean
{
  double ghalf = 1./Pi()/par[0];
  double x0 = par[1];
  double l =1.- 1./Pi()*ghalf/( pow(x[0]-x0,2.) + pow(ghalf,2.));
  if(x[0]<ltrap) l=hTrans->Interpolate(x[0]);
  return l;
}


void getSpectrum(){
  TString fileName("0.1ppm_transmission.csv");
  ifstream* in = new ifstream(fileName,std::ios::in);
  char line[1024];

  unsigned ilines=0;
  while( !(in->eof())  ) {
    in->getline(line,1024);
    TString tline(line);
    TObjArray* tokenArray = tline.Tokenize(',');
    int narray = tokenArray->GetEntries();
    //printf(" %i %s \n", narray, tline.Data() );
    if(narray<2) break;
    TObjString *sv = (TObjString*) tokenArray->At(0);
    if(ilines++==0) {
      cout << tline << endl;
      continue;
    } else {
      TObjString *sw = (TObjString*) tokenArray->At(0);
      TObjString *st = (TObjString*) tokenArray->At(1);
      TObjString *sa = (TObjString*) tokenArray->At(4);
      double a = sa->GetString().Atof();
      double w = sw->GetString().Atof();
      double t = st->GetString().Atof();
      if(a==0) continue;
      vwave.push_back(w);
      vamp.push_back(a);
      vtrans.push_back(t);
      //cout << vwave.size()-1 << " w " << w << " a " << a << endl;
    }
  }
  printf(" number of points in spectrum is %lu \n",vwave.size());
  return ;
}

void fillIshida()
{
  wishida.push_back(132.6); logIshida.push_back(1.0);//0.957
  wishida.push_back(237.8); logIshida.push_back(0.93); //0.858
  wishida.push_back(311); logIshida.push_back(0.87); //0.736
  wishida.push_back(417); logIshida.push_back(0.82);
  wishida.push_back(525); logIshida.push_back(0.76);
  wishida.push_back(632); logIshida.push_back(0.69);
  wishida.push_back(740); logIshida.push_back(0.64);

  double xzero=wishida[0];
  for(unsigned i=0; i< wishida.size(); ++i) {
    wishida[i] = (wishida[i]-xzero)/10.;
    vishida.push_back( pow(10.,logIshida[i])/10.0);
  }

  IshidaAtten = new TGraph(wishida.size(),&wishida[0],&vishida[0]);

}

void atten()
{
 
  TFile *fout = new TFile("attenOut.root","recreate");
  TString canTitle;

  getSpectrum();

  fillIshida();
  canTitle.Form("IshidaArgon");
  TCanvas *canIshida = new TCanvas(canTitle,canTitle);
  gStyle->SetOptFit();
  gPad->SetLogy();
  IshidaAtten->SetTitle("Ishida Ar attenuation ");
  IshidaAtten->GetYaxis()->SetTitle("attenuation");
  IshidaAtten->GetXaxis()->SetTitle("distance (cm) ");
  IshidaAtten->SetMarkerStyle(21);
  IshidaAtten->SetMarkerSize(.5);
  //IshidaAtten->Fit("expo");
  IshidaAtten->Draw("ap");
  canIshida->Print(".pdf");


  double wlow=120;
  double whigh=136;


  canTitle.Form("ArSpectrum");
  TGraph *arSpect = new TGraph(vwave.size(),&vwave[0],&vamp[0]);
  TCanvas *canSpect = new TCanvas(canTitle,canTitle);
  gStyle->SetOptFit();
  arSpect->SetTitle("Ar spectrum");
  arSpect->GetYaxis()->SetTitle("amplitude");
  arSpect->GetXaxis()->SetTitle("wavelength (nm) ");
  arSpect->SetMarkerStyle(22);
  arSpect->SetMarkerSize(0.4);
  arSpect->Fit("gaus");
  TF1 *gausSpectFit = (TF1*)arSpect->GetListOfFunctions()->FindObject("gaus");
  gausSpectFit->Print();

  for(int j=0; j<3; ++j) printf(" %i %f \n",j,gausSpectFit->GetParameter(j)); 
  arSpect->Draw("acp");
  canSpect->Print(".pdf");

  
  canTitle.Form("Trans");
  arTrans = new TGraph(vwave.size(),&vwave[0],&vtrans[0]);
  TCanvas *canTrans = new TCanvas(canTitle,canTitle);
  gStyle->SetOptFit();
  arTrans->SetTitle("Ar trans");
  arTrans->GetYaxis()->SetTitle("transmission");
  arTrans->GetXaxis()->SetTitle("wavelength (nm) ");
  arTrans->SetMarkerStyle(22);
  arTrans->SetMarkerSize(0.4);
  arTrans->Draw("acp"); 
  canTrans->Print(".pdf");

  // make trans hist
  hTrans = new TH1D("hTranFunc"," transmission ",vwave.size(),vwave[0], vwave[vwave.size()-1 ]);
  for(unsigned ibin=0; ibin< vwave.size(); ++ibin) hTrans->SetBinContent(ibin+1,vtrans[ibin]);

  TF1* arSpectNorm= new TF1("arSpectNorm","gaus",vwave[0],  vwave[vwave.size()-1 ]);
  arSpectNorm->SetParameter(0,1.);
  arSpectNorm->SetParameter(1,gausSpectFit->GetParameter(1));
  arSpectNorm->SetParameter(2,gausSpectFit->GetParameter(2));
  arSpectNorm->SetTitle("ArEmission-And-Transmission");
  arSpectNorm->Print();
  for(int j=0; j<3; ++j) printf(" %i %f \n",j,arSpectNorm->GetParameter(j)); 


  canTitle.Form("ArTransFunction");
  TCanvas *canTransFunc = new TCanvas(canTitle,canTitle);
  canTransFunc->SetGrid();
  arSpectNorm->GetYaxis()->SetTitle("Transmission");
  arSpectNorm->GetXaxis()->SetTitle("wavelength nm");
  arSpectNorm->Draw("");
  hTrans->Draw("same");
  arSpectNorm->SetLineColor(kRed);


  // original line 
  far =  new TF1("far","[0]/sqrt(2.*Pi())/[2]*exp(-0.5*((x-[1])/[2])**2)",wlow,whigh);
  //far = new TF1("far",Lorentz,wlow,whigh,2);
  double arMax = 2./Pi()/FWHMar3;  // peak value related to FWHM
  far->SetParameters(1.,lar3,war3);
  far->SetLineColor(kBlue);


  double wtrap=1./Pi();

  double att1 = 1./Pi()/wtrap; // 11.6  
  flor = new TF1("flor",Lorentz,wlow,whigh,2);
  flor->SetParameters(att1,ltrap);
  flor->SetParNames("max","mean");
  flor->SetLineColor(kGreen);

  fatt = new TF1("fatt",MLorentz,wlow,whigh,2);
  fatt->SetParameters(att1,ltrap);
  fatt->SetParNames("max","mean");
  fatt->SetLineColor(kBlue);

  printf(" >>>>>>>  Lor integral %f Atten integral %f <<<<<<< \n",flor->Integral(wlow,whigh),fatt->Integral(wlow,whigh));

  
  canTitle.Form("trans-norm-%.4f",att1);
  TCanvas *canlor = new TCanvas(canTitle,canTitle);
  fatt->Draw("C");
  flor->Draw("Csame");

  fout->Append(flor);
  fout->Append(far);
  fout->Append(fatt);

  printf(" trap width %f max %f \n",wtrap,fatt->Eval(ltrap));

  // make iterated line shape
  enum {NPOINT=100};
  xline.resize(NPOINT);
  wline.resize(NPOINT);
  double deltaLine = (whigh-wlow)/double(wline.size());
  enum {NTRY=100};

  TGraph *gtry[NTRY];

  // zeroth iter
  for(unsigned il=0; il<wline.size(); ++il) {
    wline[il]=wlow+double(il)*deltaLine;
    xline[il]=far->Eval(wline[il]);
  }

  //for(unsigned il=0; il<wline.size(); ++il) printf(" %i %f %f \n",il,wline[il],xline[il]);
  
  for(int itry = 0; itry<NTRY; ++itry) {
    double sum=0;
    for(unsigned il=0; il<wline.size(); ++il) { 
      if(itry>0) xline[il] *= fatt->Eval(wline[il]); 
      sum+=xline[il];
    }
    gtry[itry] = new TGraph(wline.size(),&wline[0],&xline[0]);
    gtry[itry]->SetTitle(Form("line-iter-%i",itry));
    gtry[itry]->SetMarkerSize(.7);
    gtry[itry]->SetMarkerColor(kBlack+itry);
    xvec.push_back(double(itry)*l1);
    avec.push_back(sum);
  }

   // normalize sum
  double anorm = avec[0];
  for(unsigned it=0; it<  avec.size(); ++it) avec[it]=avec[it]/anorm;

  for(unsigned it=0; it<  avec.size(); ++it) vabs.push_back(1-avec[it]);


  canTitle.Form("atten-lor-norm-%.4f",att1);
  TCanvas *can = new TCanvas(canTitle,canTitle);
  //gPad->SetLogy();
  //far->GetYaxis()->SetRangeUser(0,1);
  //fatt->GetYaxis()->SetRangeUser(0,1);
  gtry[1]->Draw("AP");
  far->Draw("Csame");
  fatt->Draw("Csame");
  can->Print(".pdf");


  canTitle.Form("transByStep-lor-norm-%.1f-%.1f-steps-%i",att1,l1,NTRY);
  TCanvas *cans = new TCanvas(canTitle,canTitle);
  TString gtitle;
  gPad->SetLogy(0);
  for(int jt=0; jt<NTRY; ++jt) {
    int it= jt;
    gtitle.Form("step length %.1f step number %i",l1,jt);
    gtry[it]->SetTitle(gtitle);
    gtry[it]->GetXaxis()->SetTitle("wavelength [nm] ");
    gtry[it]->GetYaxis()->SetTitle("transmission");
    gtry[it]->GetXaxis()->SetRangeUser(120,140);
    gtry[it]->SetMarkerColor(kBlack+it);
    gtry[it]->SetLineColor(kBlack+it);
    gtry[it]->SetLineStyle(it);
    if(jt==0) gtry[it]->Draw("APC");
    else gtry[it]->Draw("PCsame");
  }
  cans->Print(".pdf");


  // write out waveforms
  TString fileName; fileName.Form("transByStep-lor-norm-%.1f-%.1f-steps-%i.txt",att1,l1,NTRY);
  textFile.open(fileName);
  for(int jt=0; jt<NTRY; ++jt) writeGraph(gtry[jt]);
  textFile.close();


  canTitle.Form("transByStep-lor-far-norm-%.2f",att1);
  TCanvas *cans2 = new TCanvas(canTitle,canTitle);
  gPad->SetLogy(0);
  for(int jt=0; jt<NTRY; ++jt) {
    int it= jt;//NTRY-100+jt;
    gtry[it]->GetXaxis()->SetTitle("wavelength [nm] ");
    gtry[it]->GetYaxis()->SetTitle("transmission");
    gtry[it]->GetXaxis()->SetRangeUser(120,140);
    gtry[it]->SetMarkerColor(kBlack+it);
    gtry[it]->SetLineColor(kBlack+it);
    gtry[it]->SetLineStyle(it);
    if(jt==0) gtry[it]->Draw("APC");
    else gtry[it]->Draw("PCsame");
  }
  cans2->Print(".pdf");


   
  TF1* fita = new TF1("fita","expo",400,1000);
  TF1* fitb = new TF1("fitb","expo",l1*double(NTRY-100),l1*double(NTRY));

    
  fishida= new TF1("fishida","[0]*Exp(-x/[1])",xvec[0],xvec[xvec.size()-1]);
  fishida->SetParameters(1.,66);
  fishida->SetLineColor(kGreen);


  canTitle.Form("transVdist-lor-norm-%.2f",att1);
  TGraph *adist = new TGraph(xvec.size(),&xvec[0],&avec[0]);
  double pointx,pointy;
  adist->GetPoint(0,pointx,pointy);
  printf(" \t point 0 sum is %f %f \n",avec[0],pointy);
  TCanvas *cana = new TCanvas(canTitle,canTitle);
  gPad->SetLogy();
  gStyle->SetOptFit();
  adist->SetTitle("transmission versus distance");
  adist->GetXaxis()->SetTitle("dist [cm]");
  adist->GetYaxis()->SetTitle("transmission");
  adist->SetMarkerStyle(22);
  adist->SetMarkerSize(1);
  adist->Fit(fita,"R");
  adist->Draw("ac");
  fishida->SetTitle("Isida");
  fishida->Draw("same");
  IshidaAtten->SetTitle("IsidaAtten");
  IshidaAtten->Draw("psame");
  cana->BuildLegend();
  cana->Print(".pdf");

  adist->GetPoint(0,pointx,pointy);
  printf(" \t point 0 sum is %f x %f y %f \n",avec[0],pointx,pointy);

  TF1* fitDist =  new TF1("fitDist",fitAtten,l1*double(NTRY-100),l1*double(NTRY),2);
  fitDist->SetParNames("norm","tau");
  fitDist->SetParameters(1,1);


  TF1* fitA =  new TF1("fitA",fitAtten,l1*double(NTRY-100),l1*double(NTRY),2);
  fitA->SetParNames("norm","tau");
  fitA->SetParameters(1,1);


 
  canTitle.Form("attVdist-lor-norm-%.2f",att1);
  TGraph *babs = new TGraph(xvec.size(),&xvec[0],&vabs[0]);
  TCanvas *canc = new TCanvas(canTitle,canTitle);
  //gPad->SetLogy();
  gStyle->SetOptFit();
  babs->SetTitle("attenuation versus distance");
  babs->GetXaxis()->SetTitle("dist [cm]");
  babs->GetYaxis()->SetTitle("attenuation");
  babs->GetYaxis()->SetRangeUser(0.1,1.0);
  babs->SetMarkerStyle(22);
  babs->SetMarkerSize(0.7);
  //babs->Fit(fitA,"","",40,1000);
  canc->SetGridx();
  canc->SetGridy();

  babs->Draw("ac");
  //canb->BuildLegend();
  canc->Print(".pdf");



  TF1* fitd = new TF1("TwoExp",TwoExp,0,1000,3);
  fitd->SetParName(0,"A");
  fitd->SetParName(1,"lambda1");
  fitd->SetParName(2,"lambda2");
  fitd->SetParameter(0,.5);
  fitd->SetParLimits(0,0,1);
  fitd->SetParameter(1,30);
  fitd->SetParameter(2,60);
  fitd->SetLineColor(kRed);
  fitd->SetNpx(10000);



  canTitle.Form("transVdistb-lor-norm-%.2f",att1);
  TGraph *bdist = new TGraph(xvec.size(),&xvec[0],&avec[0]);
  TCanvas *canb = new TCanvas(canTitle,canTitle);
  //gPad->SetLogy();
  gStyle->SetOptFit();
  bdist->SetTitle("transmission versus distance");
  bdist->GetXaxis()->SetTitle("dist [cm]");
  bdist->GetYaxis()->SetTitle("transmission");
  bdist->GetYaxis()->SetRangeUser(0.1,1.0);
  bdist->SetMarkerStyle(22);
  bdist->SetMarkerSize(0.7);
  bdist->SetLineWidth(3);
  //bdist->Fit(fitd,"R");
  //bdist->Fit(fitDist,"","",40,100);
  canb->SetGridx();
  canb->SetGridy();
  bdist->Fit(fitd,"RM");
  bdist->Draw("ap");
  fitd->Draw("same");
  //canb->BuildLegend();
  canb->Print(".pdf");



  for(unsigned j=0; j<10; ++j) {
    double lamb = l1/(avec[j]-avec[j+1]);
    double sigma = 476./lamb;
    printf(" step %i dist %.2f atten %.2E lambda %.3f sigma %.3f [Mb] \n",j,xvec[j],avec[j],lamb,sigma);
  }


  double c1 =       fitd->GetParameter(0);
  double lambda1 =  fitd->GetParameter(1);
  double lambda2 =  fitd->GetParameter(2);

  double frac1 = c1*lambda1/(c1*lambda1 + (1-c1)*lambda2);

  printf(" fitd integral %f  frac1 %f \n",fitd->Integral(0,10000),frac1);


  

  //fout->ls();
  fout->Write();

}
