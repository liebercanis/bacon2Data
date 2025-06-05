TRandom3 *ran;

void minuitTest(double slope =1 , double intercept = 0)
{
  TFile *fout = new TFile("minuitTest.root","recreate");
  int ngen = 1.E2;
  double xlow=-1.;
  double xhigh=1.;
  double sigma = 0.1;
  TNtuple *nt= new TNtuple("ntLine","ntLine","x:y");
  ran = new TRandom3();
  TF1 *f1 = new TF1("myLine","[0]*x+[1]",xlow,xhigh);
  f1->SetParName(0, "slope");
  f1->SetParameter(0,slope);
  f1->SetParName(1, "intercept");
  f1->SetParameter(1, intercept);

  std::vector<double> xval;
  std::vector<double> yval;
  std::vector<double> yerr;

  // generate points about line

  for(int igen=0; igen<ngen; ++igen)
  {
    double x = (xhigh-xlow)*ran->Rndm() + xlow;
    double ymean = f1->Eval(x);
    double y = ran->Gaus(ymean,sigma);
    xval.push_back(x);
    yval.push_back(y);
    yerr.push_back(sigma);

    nt->Fill(x,y);
  }
  TGraphErrors *gr =new TGraphErrors(xval.size(), &xval[0],&yval[0], nullptr,&yerr[0]);
  gr->SetName("ranLine");
  gr->SetTitle("random points about a Line");
  gr->Fit(f1,"myLine","R+",xlow,xhigh);



  TCanvas *c = new TCanvas("minuitTest", "minuitTest");
  gr->Draw("ap");
  f1->Draw("same");
  //gStyle->SetOptStat(1001101);
  gStyle->SetOptFit(1011);
  c->Update();
  
  fout->Add(f1);
  fout->Add(gr);
  fout->Write();
  
}
