// Macro myfunc.C
TF1* f1;  // file scope 
double mySinFunction(double *x, double *par)  // x is f(x) with paramters passed in
{
   double xx =x[0];
   double f = TMath::Abs(par[0]*sin(par[1]*xx)/xx);
   return f;
}

// same name as file name
void myFunction()
{
   f1 = new TF1("myfunc",mySinFunction,0,10,2);  // "mySinFunction is the function pointer range of function begin to end with 2 parameters
   f1->SetParameters(2,1);
   f1->SetParNames("constant","coefficient");
   TCanvas *canf = new TCanvas("func","func");
   f1->Draw();
   auto h1 = new TH1F("h1","test",100,0,10);
   cout << f1->GetParameter(0) << " " <<  f1->GetParameter(0) << endl;

   h1->FillRandom("myfunc",20000);
   cout << " first bin "  << h1->GetBinContent(1) << endl;
   // stating for fit 
   f1->SetParameters( 1.1*h1->GetBinContent(1),1.2);
   // call minimization MINUIT
   h1->Fit("myfunc");
   TCanvas *canh = new TCanvas("hist","hist");
   h1->Draw();
}
