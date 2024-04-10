// path is average path length in cm
void absFact(double path=50.0) {
  enum {NP=6};
  double PPM[NP] = {0.1,1,2,5,10,100};
  double a[NP];

  double l0 = 12.7; // at 0.1 ppm
  for(int i=0; i<NP; ++i) {
    double li = l0*PPM[i]/0.1;
    a[i] = TMath::Exp(-path/li);
    printf(" %i ppm %f li %f a %f \n",i,PPM[i],li,a[i]);
  }

  TString ctitle;
  ctitle.Form("absFact-path-%0.f-cm",path);
  TCanvas *can = new TCanvas(ctitle,ctitle);
  can->SetGridx();
  can->SetGridy();
  can->SetLogx();

  TGraph *gr = new TGraph(NP,&PPM[0],&a[0]);
  gr->SetTitle(ctitle);
  gr->GetXaxis()->SetTitle("PPM");
  gr->GetYaxis()->SetTitle("absorption constant A");
  gr->SetMarkerStyle(22);
  gr->SetMarkerSize(1.0);
  gr->Draw("ap");

}
