
vector<double> vtime;
vector<double> va;
vector<double> vb;

unsigned maxline=10000;

void getSpectrum(){
  ifstream in("delta.dat",std::ios::in);
  unsigned ilines=0;
  double a;
  while( !(in.eof())  ) {
    in >> a;
    //cout << ++ilines << " " << a <<  endl;
    if(a==0) break;
    //a+=1;
    if(ilines>maxline) break;
    vtime.push_back(double(ilines++));
    va.push_back(double(a));
  }
  in.close();
  return ;
}


void getSpectrum2(){
  ifstream in("circuit.dat",std::ios::in);
  unsigned ilines=0;
  double b;
  while( !(in.eof())  ) {
    in >> b;
    //cout << ++ilines << " " << a <<  endl;
    if(b==0) break;
    //a+=1;
    if(ilines>maxline) break;
    //cout << "\t" << ilines << " " << b << endl;
    vb.push_back(double(b));
  }
  in.close();
  return ;
}




void response(){
  getSpectrum();
  getSpectrum2();

  double sum=0;
  double sumb=0;
  double minval=1E9;
  double maxval=0;
  for(unsigned i=0; i< vtime.size() ; ++i) {
    sum += va[i];
    sumb += vb[i];
    if(va[i]<minval) minval=va[i];
    if(va[i]>maxval & i>0) maxval=va[i];
  }
  cout << " delta " << va[0] << " sum " << sum << " b " << sumb << " to time  " << vtime[vtime.size()-1] <<  endl;
  printf("  va sum  %f vb sum %f \n",std::accumulate(va.begin(), va.end(),0.0),std::accumulate(vb.begin(), vb.end(),0.0));
  printf(" va[0,1]= [%f.%f]  vb[0,1]= [%f.%f] \n",va[0],va[1],vb[0],vb[1]);

  for(unsigned i=0; i< vtime.size() ; ++i) {
    va[i]+=1; 
    vb[i]+=1; 
  }

  printf("  starting %f %f \n",va[0],vb[0]);


  TGraph *gr = new TGraph(vtime.size(),&vtime[0],&va[0]);
  TGraph *gb = new TGraph(vtime.size(),&vtime[0],&vb[0]);
  TCanvas *can = new TCanvas("delta","delta");
  
  gPad->SetLogy();
  //gPad->SetLogx();
  gr->SetTitle("delta");
  //gr->GetYaxis()->SetRangeUser(.99*minval,1.01*maxval);
  gr->GetYaxis()->SetTitle("reponse");
  gr->GetXaxis()->SetTitle("time");
  gr->SetMarkerStyle(22);
  gr->SetMarkerSize(.2);
  gr->Draw("ap");
   //gPad->SetLogx();
  gb->SetTitle("deltab");
  //gb->GetYaxis()->SetRangeUser(.99*minval,1.01*maxval);
  gb->GetYaxis()->SetTitle("reponse");
  gb->GetXaxis()->SetTitle("time");
  gb->SetMarkerStyle(22);
  gb->SetMarkerSize(.2);
  gb->SetMarkerColor(kRed);
  gb->Draw("psame");


}
