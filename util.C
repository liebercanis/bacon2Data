void browse() 
{
 TBrowser* browser = new TBrowser();
 gSystem->Load("libTreeViewer");
 printf(" browse \n");
}

void dir() 
{
  gDirectory->ls("-m");
}
//create hint1 filled with the bins integral of h1
 TH1F *integralHist(TH1F* hf) {
   TString hname = TString(hf->GetName())+ TString("Integral");
   TString htitle = hf->GetTitle();
   Int_t nbins = hf->GetNbinsX();
   Float_t first =  hf->GetBinLowEdge(hf->GetXaxis()->GetFirst());
   Float_t last =   hf->GetBinLowEdge(hf->GetXaxis()->GetLast()+1);
   printf(" integralHist entries %i making %s %s %i %f %f \n",int(hf->GetEntries()),hname.Data(),htitle.Data(),nbins,first,last);

   TH1F *hint1 = new TH1F(hname,htitle,nbins,first,last);
   Float_t sum = 0;
   for (Int_t i=1; i<=nbins ; i++) {
      sum += hf->GetBinContent(i);
      hint1->SetBinContent(i,sum);
   }
   return hint1;
 }

