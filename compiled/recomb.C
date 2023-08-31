{
  using namespace TMath;
  double Avagadro = 6.02e23;
  double rhoAr = 1.3954;              //% Density of liquid Argon g/cm^3
  double mrhoAr = rhoAr/40;           //% Molar density of LAr
  double NTotal = 30000;

  TFile *fout = new TFile("recomb.root","recreate");

// recombination rate from https://doi.org/10.1016/S0009-2614(02)01177-6 Mariusz Wojcik, M.Tachiya 
  double k_rec = 1.6e-5 * 1e-6 * mrhoAr * Avagadro * 1e-9;  //%(s^-1)
  printf("recombination k %E tau %E \n",k_rec,1./k_rec);


  TF1 *fs = new TF1("fsinglet","[0]*Exp(-x/[1])",0,400);
  fs->SetParameters(0.15,7.);
  fs->SetTitle("singlet");
  
  TF1 *ft = new TF1("ftriplet","[0]*Exp(-x/[1])",0,400);
  ft->SetParameters(0.85,1300.);
  ft->SetTitle("triplet");


  double ar = 2.0;
  TF1 *fr = new TF1("frecomb","[0]*pow(1+x/[1],-2)",0,400);
  fr->SetParameters(ar,37.);
  fr->SetTitle("recomb");


  TF1* fsum = new TF1("fsum","[0]*Exp(-x/[1]) + [2]*Exp(-x/[3])+[4]*pow(1+x/[5],-2)",0,400);
  fsum->SetParameters(0.15,7,0.85,1300,ar,37);
  fsum->SetTitle("sum");


  fs->SetLineColor(kRed);
  ft->SetLineColor(kBlue);
  fr->SetLineColor(kGreen);
  fsum->SetLineColor(kBlack);



  TCanvas *can = new TCanvas("reco","reco");
  fsum->Draw();
  fsum->GetXaxis()->SetTitle("time (ns)");
  fsum->GetYaxis()->SetRangeUser(1.E-2,2.);
  fsum->Draw();
  fs->Draw("sames");
  ft->Draw("sames");
  fr->Draw("sames");
  can->BuildLegend();

  fout->Add(fsum);
  fout->Add(fs);
  fout->Add(ft);
  fout->Add(fr);

  fout->Write();
 
} 
