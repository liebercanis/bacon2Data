// times in ns
const double alpha3=0.86;
double t3=2100;
double ts=5;
double kxe = 8.8E-5;
double kp = 0;// 1.3E-4;

/*I fix all of the lifetimes: microseconds*/
/*I fix all of the lifetimes: microseconds*/
static double pmodtS = 7E-3;
//static double pmodtT = 1.4;  //1.96;
static double pmodtG = 3.48;
static double pmodtX = 20E-3;
static double xstart=1.05;
static double xstop=10;
static double tailStart=7.0;
static double tailStop=10.0;
static double tiny = 1.0E-9;

using namespace TMath;

// fmodelFit[ifit]->SetParNames("binwidth","T0","taut","k1","tauX","Am","tauM","ppm0","ppm","type");
static double fmodel(double *xx, double *par) 
{
  double  t= xx[0]-xstart;
  if(t<0) return tiny;
  int type = int(par[9]);
  double binwidth = par[0];
  double ta = par[2];
  double ppmTrue = (par[8]-par[7]);
  double k = par[3]*ppmTrue;
  double tr=1.0;
  //if(ppmTrue>.001) tr= (35. - 2.54*TMath::Log(5.8*ppmTrue))/48.0;
  double tx = par[4];
  double tg = par[6];
  double tap = 1./( 1./ta + k);
  double cnorm = 1;//Exp(-xstart/tap);
  double gnorm = 1.-Exp(-(xstop-xstart)/tg);
  double A0 = par[1]*binwidth/cnorm;
  double A = A0*Exp(-t/tap)/ta;
  double G0 = par[5]*binwidth/gnorm;
  double t1=1./( 1./tx - k);
  double t2=1./( 1./tx- 1./tap);
  double C = A0*ta*k;
  double D = C*( Exp(-t*k) - Exp(-t/tap) );
  double X = C/tx*k*(t1*Exp(-t*k) - t2*Exp(-t/tap)) + C/tx*k*(t2-t1)*Exp(-t/tx);
  double G = G0/tg*Exp(-t/tg);
  double f;
  if(type ==0) f=A;
  else if(type ==1) f=D;
  else if(type ==2) f=X;
  else if(type ==3) f=G;
  else if(type ==4) f=A+X;  
  else f=A+X+G; 
  return f;
}



static Double_t light(Double_t *xx, Double_t *par)
{ 
  double kx = kxe*par[0];
  double td = 1/kx;
  double tq = 1./(1./t3+kp);
  double tr = 1./(kx + 1./tq);
  double x=xx[0];

  double fs = (1-alpha3)* Exp( -x/ts)/ts;
  double f3 = alpha3/t3*Exp(-x/tr);
  double fx =  alpha3*pow(kx,2)*tq*( Exp( -x/td) - Exp( -x/tr) );
  return fs+f3+fx;
}

static Double_t yield(Double_t ppm)
{ 
  double kx = kxe*ppm;
  double td = 1/kx;
  double tq = 1./(1./t3+kp);
  double tr = 1./(kx + 1./tq);

  double y = (1-alpha3) + alpha3*tr*(kx+1./t3);
  return y;
}



void model()
{

  double PPM[3]={0,10.,100.};

  double ppm=PPM[1];
  double kx = kxe*ppm;
  double td = 1./kxe*ppm;
  double tq = 1./(1./t3+kp);
  double tr = 1./(kx + 1./tq);

  double cut = 100;
  double max = 3000;


  printf(" ts %.2E t3 %.2E tr %.2E tq %.2E td %.2E @ %.1f PPM \n\n",ts,t3,tr,tq,td,ppm);



  TF1* f10  = new TF1("ften",light,0,max,1);
  f10->SetParName(0,"PPM");
  f10->SetParameter(0,10.);
  f10->SetTitle(Form("light%.0f-ppm",PPM[1]));

  
  TF1* f20  = new TF1("ftwenty",light,0,max,1);
  f20->SetParName(0,"PPM");
  f20->SetParameter(0,20.);
  f20->SetTitle(Form("light%.0f-ppm",PPM[1]));


  TF1* fp  = new TF1("fpure",light,0,max,1);
  fp->SetParName(0,"PPM");
  fp->SetParameter(0,0.001);
  fp->SetTitle(Form("light%.0f-ppm",PPM[1]));

  fp->SetLineColor(kBlack);
  f20->SetLineColor(kRed);
  f10->SetLineColor(kGreen);


  TString canTitle;
  canTitle.Form("light%.0f-ppm",ppm);
  TCanvas *can = new TCanvas(canTitle,canTitle);
  gPad->SetLogy();
  fp->Draw();
  f10->Draw("same");
  f20->Draw("same");

  
  double y10 = f10->Integral(0,max);
  double y10fast = f10->Integral(0,cut);
  double y10long = f10->Integral(cut,max);

  double y20 = f20->Integral(0,max);
  double y20fast = f20->Integral(0,cut);
  double y20long = f20->Integral(cut,max);


  double yp = fp->Integral(0,max);
  double ypfast = fp->Integral(0,100);
  double yplong = fp->Integral(100,max);
  printf(" int time %.0f (>%.0f cut) yp %.2E  ratio %.0f PPM %.2E  %.0f PPM  %.2E \n",max,cut,yplong,PPM[1],y10long/yplong,PPM[2],y20long/yplong);
  printf(" int time %.0f total  yp %.2E ratio %.0f PPM %.2E  %.0f PPM  %.2E \n", max, yp, PPM[1],y10/yp,PPM[2],y20/yp);

  printf(" L from eq 23 : 0 %.3f (%.3f)  10PPM %.3f  100PPM  %.3f \n",yield(0),yp, yield(10),yield(100));


  printf(" \n\n >>> setFit %i  %.f PPM fit tauG %.3f integral %.3E  G0Par %.3E T0 %.3E IFIT %i \n\n\n",iset,ppm,tauG,integral,G0Par,A0,ifit);
  TF1 *fp= new TF1(Form("LfFit-%i-%0.f-PPM",iset,ppm),fmodel,xstart,xstop,npars);
  fp->SetParNames("binwidth","T0","taut","k1","tauX","G0","tauG","ppm0","ppm","type");
  fp->SetParameters(binwidth,A0, pmodtT, k1 ,pmodtX,G0Par,tauG, ppm0,  ppm , ifit );
  fp->SetNpx(1000); // numb points for function
  fp->FixParameter(0,binwidth);
  if(fixTriplet) fp->FixParameter(2,pmodtT);


}
