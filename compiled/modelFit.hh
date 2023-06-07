// time is in microseconds
// model with absorption
//static double tres = 7.86;
static double tres = 5.4;
static double tTriplet = 1600.0; //2100.0;
static double tSinglet = 5.0;
static double tMix = 4700.;
static double tXe = 20.0;
static double kxe = 8.8E-5;
static double kplusZero = 1.3E-4;
static double xTrigger = 1000.;
static double xMax=10000.;

enum {NSETS =12};
enum {NTYPES=5};



static double expGaus(double x, double tau)
{
  x = x + 0.33*tres;
  double arg1 = (tres * tres / tau - 2. * x) / 2. / tau;
  double arg2 = (tres * tres / tau - x) / sqrt(2) / tres;
  double f = 0.5 * TMath::Exp(arg1) * TMath::Erfc(arg2);
  return f;
}

static double lightModel(Double_t *xx, Double_t *par)
{
  double x = xx[0] - xTrigger;
  double bw =  par[0];
  double norm = par[1];
  double ppm=  par[2];
  double tTrip = par[3];
  double kp = par[4]*kplusZero;
  double tmixPar = par[9];
  double lmix = 1./tmixPar;
  //double alpha1 = bw*norm;
  //double alpha3 = (1.-par[5])/par[5]*alpha1;
  double alpha1 = par[5]*bw*norm;
  double alpha3 = (1.-par[5])*bw*norm;
  double ab= par[6];
  double k1Zero = kxe * 131. / 40.;
  double kx = k1Zero*ppm;
  double kPrime = lmix+kx+kplusZero*par[7];
  int type = int(par[8]);

  double lS = 1./tSinglet;
  double lT = 1/ tTrip;
  double l1 = 1./tSinglet + kp + kx;
  double l3 = 1./tTrip    + kp + kx;
  double lX = 1./tXe;
  double c1 = kx + ab*lS;
  double c3 = kx + ab*lT;
  double t1 = 1./l1;
  double t3 = 1./l3;
  double tkPrime = 1./kPrime;

  double f=0;
  
  // model
  double fs = (1.-ab)*alpha1 / tSinglet * expGaus(x, t1);
  double ft = (1.-ab)*alpha3/tTrip * expGaus(x, t3);
  double fm = alpha1*c1/(l1-kPrime)*(  expGaus(x, tkPrime)  - expGaus(x, t1) ) + alpha3*c3/(l3-kPrime)*(  expGaus(x, tkPrime)  - expGaus(x, t3) );
  fm /= tmixPar;

  double x1 = c1*kx*alpha1/(l1-kPrime)/tXe * (  (expGaus(x, tkPrime) - expGaus(x, tXe))/(lX - kPrime)  -  (expGaus(x, t1)  -  expGaus(x, tXe))/(lX - l1));
  double x3 = c3*kx*alpha3/(l3-kPrime)/tXe * (  (expGaus(x, tkPrime) - expGaus(x, tXe))/(lX - kPrime)  -  (expGaus(x, t3)  -  expGaus(x, tXe))/(lX - l3));
  double fx = x1 + x3;
 
  /*segreto 
  double td = 1 / kx;
  double tq = 1. / (1. / tTrip + kp);
  double tr = 1. / (kx +  1. / tq);
  double fs = alpha1 / tSinglet * expGaus(x, tSinglet);
  double ft = alpha3 /tTrip * expGaus(x, tr);
  double fx = alpha3 * pow(kx, 2) * tq * (expGaus(x, td) - expGaus(x, tr));
  double fm=fx;
  */

  
  
  
  
  /*exp
  //double fs =  alpha1/tSinglet * TMath::Exp(-1.*x/tSinglet);
  //double ft = alpha3/tTrip * TMath::Exp(-1.*x/tTrip);
  double fs = alpha1 / tSinglet * expGaus(x, tSinglet);
  double ft = alpha3 /tTrip * expGaus(x, tTrip);
  double fm=0.;
  double fx=0.;
  */


  //printf(" %E %E %E %E %E  \n"  ,fs,ft,fm,x1,x3);

  if (type == 0)
    f = fs;
  else if (type == 1)
    f = ft;
  else if (type == 2)
    f = fm;
  else if (type == 3)
    f = fx;
  else if (type == 4)
    f = fs+ft+fx + fm;

  return f;

}

class modelFit
{
  public:
    //double taut = 2100.;
    modelFit(int ifit, int iset); 
    virtual ~modelFit() { ; }
    TF1 *fp;
    enum { NPARS =10 };
    double XPPM[NSETS]={.01,.05,.07,0.1,.2,.5,1,2,5,10,100,1000};
    double PPM[NSETS]={.1,1,2,5,10,100,1000,1000,1000,1000,0,0};
    double binwidth = 1.0; //0.04;
    double norm = 3.0E4;
    double sFrac = 0.2;
    int setColor[NSETS];
    int setStyle[NSETS];
    int theFit;
    int theSet;
    void show();
};

void modelFit::show() {
  {
    printf(" \n\n >>> modelFit fitted parameters fit %i set %i \n",theFit,theSet);
    for (int ii = 0; ii < NPARS; ++ii)
    {
      printf("\t  param %i %s %.3E +/- %.3E \n", ii, fp->GetParName(ii), fp->GetParameter(ii), fp->GetParError(ii));
    }
    double tTrip =  fp->GetParameter(3);
    double kp = fp->GetParameter(4)*kplusZero;
    double ppm=  fp->GetParameter(2);
    double ab=  fp->GetParameter(6);
    double tmixPar=  fp->GetParameter(9);

    double k1Zero = kxe * 131. / 40.;
    double kx = k1Zero*ppm;
    double kPrime = 1./tMix + kx+kplusZero*fp->GetParameter(7);
    double lS = 1./tSinglet;
    double lT = 1/tTrip;
    double l1 = 1./tSinglet + kp + kx;
    double l3 = 1./tTrip    + kp + kx;
    double lX = 1./tXe;
    double c1 = kx + ab*lS;
    double c3 = kx + ab*lT;
    double t1 = 1./l1;
    double t3 = 1./l3;
    double tkPrime = 1./kPrime;
    double sfrac = fp->GetParameter(5);
    double bw = fp->GetParameter(0);
    double norm = fp->GetParameter(1);
    double alpha1 = sfrac*bw*norm;
    double alpha3 = (1.-sfrac)*bw*norm; 
    double snorm = alpha1*c1*kx/(l1-kPrime)/tXe;
    double tnorm = alpha3*c3*kx/(l3-kPrime)/tXe;
    printf( "\t ls %.3E c1 %.3E  lt %.3E  c3 %.3E  tMix %.3E  kPrime %.3E  l1-kPrime %.3E  l3-kPrime %.3E \n",lS,c1,lT,c3,tmixPar, kPrime,  l1-kPrime, l3-kPrime);
    printf( "\t S norm %.3E T norm %.3E\n",snorm,tnorm);
  }
}
  
modelFit::modelFit(int ifit, int iset) 
{
  theFit=ifit; theSet = iset;
  setColor[0] = kOrange + 4;
  setColor[1] = kBlue;
  setColor[2] = kRed;
  setColor[3] = kGreen - 2;
  setColor[4] = kBlack;
  setColor[5] = kCyan+2;
  setColor[6] = kMagenta;

  for(int i=0; i<NSETS; ++i ) setStyle[i]=1;
  double ab=0;
  double ppm=PPM[iset];
  double kplus = 1;
  double kPrime = 1;

  fp = new TF1(Form("LifeFit-%.2f-type-%i-set-%i",ab,ifit,iset), lightModel, xTrigger-20.,xMax, NPARS);
  printf(" modelFit: set %i fit range %f to %f  \n", iset,xTrigger-20.,xMax);
  fp->SetParName(0,"binw");
  fp->SetParName(1,"norm");
  fp->SetParName(2,"PPM");
  fp->SetParName(3,"tau3");
  fp->SetParName(4,"kp");
  fp->SetParName(5,"sfrac");
  fp->SetParName(6,"ab");
  fp->SetParName(7,"kprime");
  fp->SetParName(8,"type");
  fp->SetParName(9,"tmix");

  fp->FixParameter(0, binwidth);
  fp->SetParameter(1, norm);
  fp->FixParameter(2, ppm);
  fp->FixParameter(3, tTriplet);
  fp->FixParameter(4,kplus);
  fp->SetParameter(5,sFrac);
  fp->SetParLimits(5,.01,.5);
  fp->SetParameter(6,ab);
  fp->SetParLimits(6,1.E-9,1.);
  fp->FixParameter(7,2*kPrime);
  fp->FixParameter(8,ifit);
  fp->FixParameter(9,tMix);
  //fp->SetParLimits(9,1.E3,20.E3);
  fp->SetTitle(Form("LfFit-type-%i-set-%i-%0.1f-PPM-ab-%.2f", ifit, iset, ppm,ab));
  fp->SetNpx(1000); // numb points for function
  fp->Print();
  fp->SetLineColor(setColor[iset]);
  fp->SetMarkerStyle(setStyle[iset]);

   cout << "defined function " << fp->GetName() << "  " <<  fp->GetTitle() << endl;
}
