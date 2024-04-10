using namespace TMath;
void voltage(int NP=1000)
{
  double tunit=1;
  double omega = 2.*Pi()/tunit*2./NP;
  TComplex Z(0.1,10.);
  double RC = 10.;
  TComplex A(0,omega*RC);
  TComplex C=A*Z;
  TComplex b = 1./(1./C+1.);
  cout << b.Re() << "  " << b.Im()  << endl;
  vector<double> time;
  vector<double> vout;
  vector<double> vin;
  for(unsigned int i=0; i< NP; ++i) {
    double t = double(i)*tunit;
    time.push_back(t);
    vin.push_back(Cos(omega*t));
    vout.push_back(b.Rho()*Cos(omega*t+b.Theta()));
    //printf(" t %f wt %f vin %f vout %f \n",t,omega*t,vin[vin.size()-1],vout[vout.size()-1]);
  }

  TGraph* gin = new TGraph(time.size(),&time[0],&vin[0]);
  TGraph* gout = new TGraph(time.size(),&time[0],&vout[0]);
  gin->SetLineStyle(2);
  gout->SetLineColor(kBlue);

  TCanvas *can = new TCanvas("voltage","voltage");
  gout->Draw("ac");
  gin->Draw("cSAME");
}
