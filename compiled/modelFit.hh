// MG added additions from Doug and geometry Jun 27, 2023
// time is in microseconds
// model with absorption
// static double tres = 7.86;
// static double tres = 5.4;
static double tres = 10;
static double tTriplet = 1600.0; // 2100.0;
static double tSinglet = 5.0;
static double tMix = 4700.;
static double tXe = 20.0;
static double kxe = 8.8E-5; //diffusion time 
static double kqZero = 1.3E-4; // kq in the paper quenching rate
static double xTrigger = 1411.; // 1200.;
static double xMin = 1300.;     // 15000.;
static double xMax = 9000.;     // 15000.;
static double nPhotons = 50.E3 * 5.486; // 274300.00
static double distanceLevel[4];
static double effGeo[13];
static double trec = 37.1; // ns Eur. Phys. J. C (2013) 73:2618

enum
{
  NTYPES = 5,
  NCHAN = 13,
  NPARS = 13 
};

/*
 Calculate absorption as a function of distance and xenon concentration.
// Taken from fits to Neumeier data.
        A = 0.615;
        lambda1 = 12.7*0.1/ppmx(i);
        lambda2 = 740*0.1/ppmx(i);
        Tr128 = A*exp(-davg/lambda1)+(1-A)*exp(-davg/lambda2);
        Absorb = 1-Tr128;
        cond1 = S(0) == Sfrac*NTotal;
        cond2 = T(0) == (1-Sfrac-Ifrac)*NTotal;
        cond3 = I(0) == Ifrac*NTotal;
        cond4 = M(0) == 0;
        cond5 = X(0) == 0;
        conds = [cond1, cond2, cond3, cond4, cond5];
// With nitrogen
        ode1 = diff(S,t) == -(1/taus + 1/taudx + 1/taudn + 1/tauquench)*S;  % One taud for the transfer to the mixed state, one for transfer to contaminants (N)
        ode2 = diff(T,t) == -(1/taut + 1/taudx + 1/taudn + 1/tauquench)*T;
        ode3 = diff(I,t) == -2*Ifrac*NTotal*tau3^2/(tau3+t)^3;
        ode4 = diff(M,t) == (1/taus)*Absorb*S + (1/taut)*Absorb*T + (1/taudx)*(S+T) + (1/tau3)*Absorb*I- (1/taum + 1/taudx + 1/taudn + 1/tauquenchmix)*M;
        ode5 = diff(X,t) == (1/taudx)*M - (1/taux + 1/taudn + 1/tauquench)*X;
        odes = [ode1, ode2, ode3, ode4, ode5];
        Solve = dsolve(odes, conds);
        SSol = Solve.S;
        TSol = Solve.T;
        MSol = Solve.M;
        XSol = Solve.X;
        ISol = Solve.I;
// Hamamatsu PMT or SiPM QE corrections
        SiPMQE150 = 0.238;
        SiPMQE175 = 0.238;
        PMTQE175 = 0.38;
*/

double effGeoFunc(int ichan)
{
  int ilevel = -1;
  if (ichan == 6 || ichan == 7 || ichan == 8)
    ilevel = 0;
  else if (ichan == 3 || ichan == 4 || ichan == 5)
    ilevel = 1;
  else if (ichan == 0 || ichan == 1 || ichan == 2)
    ilevel = 2;
  else if(ichan==12)
    ilevel = 3;

  double e = 1.0;
  if (ilevel < 0)
    return e;
  /*
  Area of SiPMs is 6.0mm x 6.0mm

      Channels 6, 7, and 8 are at 11.6 cm
      from the source Channels 3, 4, and 5 are at 23.2 cm
      from the source Channels 0, 1, and 2 are at 34.8 cm from the source Channel 12 is at 36 cm from the source.
      */
  double a = pow(0.6, 2.);
  double b = 4.0 * TMath::Pi();
  double distance2[4];
  distance2[0] = pow(distanceLevel[0], 2.);
  distance2[1] = pow(distanceLevel[1], 2.);
  distance2[2] = pow(distanceLevel[2], 2.);
  distance2[3] = pow(distanceLevel[3], 2.);

  e = a / b / distance2[ilevel];
  return e;
}

static int level(int ichan)
{
  int ilevel = -1;
  if (ichan == 6 || ichan == 7 || ichan == 8)
    ilevel = 0;
  else if (ichan == 3 || ichan == 4 || ichan == 5)
    ilevel = 1;
  else if (ichan == 0 || ichan == 1 || ichan == 2)
    ilevel = 2;
  else if (ichan == 12)
    ilevel = 3;
  return ilevel;
}

static double QEff128(double ppm, double dist)
{
  return 0.16 - (10 - ppm * dist) * 0.002;
}

static double Absorption(double ppm, double dist)
{
  // Calculate absorption as a function of distance and xenon concentration.% Taken from fits to Neumeier data.
  double A = 0.615;
  double lambda1 = 12.7 * 0.1 / ppm;
  double lambda2 = 740 * 0.1 / ppm;
  double Tr128 = A * exp(-dist / lambda1) + (1 - A) * exp(-dist / lambda2);
  return 1. - Tr128;
}


static double expGaus(double x, double tau)
{
  double arg1 = (tres * tres / tau - 2. * x) / 2. / tau;
  double arg2 = (tres * tres / tau - x) / sqrt(2) / tres;
  double f = 0.5 * TMath::Exp(arg1) * TMath::Erfc(arg2);
  return f;
}

static double lightModel(Double_t *xx, Double_t *par)
{
  int ichan = int(par[10]);
  double bw = par[11];
  int ifit = int(par[12]);
  double x = xx[0] - xTrigger;
  double norm = par[0];
  double ppm = par[1];
  double tTrip = par[2];
  double kp = par[3] * kqZero;
  double tmixPar = par[8];
  double lmix = 1. / tmixPar;
  double bkg = par[9];
  // double alpha1 = bw*norm;
  // double alpha3 = (1.-par[5])/par[5]*alpha1;
  double sfrac = par[4];
  double rfrac = par[5];
  double alpha1 = sfrac * bw * norm * effGeo[ichan];
  //printf(" %f %f %f %f %f .....", sfrac, bw, norm, effGeo[ichan],alpha1);
  double frec = bw * norm * effGeo[ichan] * pow(1 + x / trec, -2.);
  double alpha3 = (1. - sfrac -rfrac) * bw * norm * effGeo[ichan];
  double ab = par[6]; 
  double k1Zero = kxe * 131. / 40.;
  double kx = k1Zero * ppm;
  double kxPrime = lmix + kx + kqZero * par[7];

  double lS = 1. / tSinglet;
  double lT = 1 / tTrip;
  double l1 = 1. / tSinglet + kp + kx;
  double l3 = 1. / tTrip + kp + kx;
  double lX = 1. / tXe;
  double c1 = kx + ab * lS;
  double c3 = kx + ab * lT;
  double t1 = 1. / l1;
  double t3 = 1. / l3;
  double tkxPrime = 1. / kxPrime;

  // model
  double fs = (1. - ab) * alpha1 / tSinglet * expGaus(x, t1);
  double ft = (1. - ab) * alpha3 / tTrip * expGaus(x, t3);
  double fm = alpha1 * c1 / (l1 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t1)) + alpha3 * c3 / (l3 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t3));
  fm /= tmixPar;

  //printf("  %E %E  %E  \n",fs,x,expGaus(x,t1));

  double x1 = c1 * kx * alpha1 / (l1 - kxPrime) / tXe * ((expGaus(x, tkxPrime) - expGaus(x, tXe)) / (lX - kxPrime) - (expGaus(x, t1) - expGaus(x, tXe)) / (lX - l1));
  double x3 = c3 * kx * alpha3 / (l3 - kxPrime) / tXe * ((expGaus(x, tkxPrime) - expGaus(x, tXe)) / (lX - kxPrime) - (expGaus(x, t3) - expGaus(x, tXe)) / (lX - l3));
  double fx = x1 + x3;

  // printf(" %E %E %E %E %E  \n"  ,fs,ft,fm,x1,x3);
  // QE eff factors
  // Tot = (SiPMQE128*SLight + SiPMQE128*TLight + SiPMQE128*ILight + SiPMQE175*XLight + SiPMQE150*MLight);
  int ilevel = level(ichan);
  double dist = distanceLevel[ilevel];
  double SiPMQ128 = QEff128(ppm, dist);
  double SiPMQE150 = 0.238;
  double SiPMQE175 = 0.238;
  double PMTQE175 = 0.38;

  fs = fs * SiPMQ128;
  ft = ft * SiPMQ128;
  fm = fm * SiPMQE150;

  // PMT sees only  175
  if (ichan == 5)
  {
    fx = fx * SiPMQE175;
    fs = 0;
    ft = 0;
    fm = 0;
    frec = 0;
  }
  else if (ichan == 12)
  {
    fx = fx * PMTQE175;
    fs = 0;
    ft = 0;
    fm = 0;
    frec = 0;
  }
  else
  { // sipms do not see 175
    fx = fx * SiPMQE175;
  }
  double f = 0;
  if (ifit == 0)
    f = fs;
  else if (ifit == 1)
    f = ft;
  else if (ifit == 2)
    f = fm;
  else if (ifit == 3)
    f = fx;
  else if (ifit == 4)
    f = fs + ft + fx + fm;
    //f = fs + ft + fx + fm +frec;

  f = f + bkg;

  return f;
}

class modelFit
{
public:
  // double taut = 2100.;
  modelFit(int thefit, int ichan, double ppm);
  virtual ~modelFit() { ; }
  TF1 *fp;
  double binwidth = 2.0; // ns;
  double norm = 3.0E4;
  double sFrac = 0.2;
  int theChan;
  double thePPM;
  void show();
  void showEff();
};

void modelFit::showEff() {
  cout << " modelFit::showEff " << endl;
  for (int ichan = 0; ichan < 13; ++ichan)
  {
    int ilevel = level(ichan);
    double d = distanceLevel[ilevel];
    double e = effGeo[ilevel];
    printf("chan %i level %i dist %.3f e %.3E \n", ichan, ilevel, d, e);
    //printf("effGeo[%i]=%f ; \n", ichan, e);
  }
}

void modelFit::show()
{
  printf(" \n\n >>> modelFit fitted parameters fit chan %i ppm %f \n", (int)fp->GetParameter(10), fp->GetParameter(1));
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  param %i %s %.4E +/- %.4E \n", ii, fp->GetParName(ii), fp->GetParameter(ii), fp->GetParError(ii));
  }
  double tTrip = fp->GetParameter(3);
  double kp = fp->GetParameter(4) * kqZero;
  double ppm = fp->GetParameter(2);
  double ab = fp->GetParameter(6);
  double tmixPar = fp->GetParameter(9);

  double k1Zero = kxe * 131. / 40.;
  double kx = k1Zero * ppm;
  double kxPrime = 1. / tMix + kx + kqZero * fp->GetParameter(7);
  double lS = 1. / tSinglet;
  double lT = 1 / tTrip;
  double l1 = 1. / tSinglet + kp + kx;
  double l3 = 1. / tTrip + kp + kx;
  double lX = 1. / tXe;
  double c1 = kx + ab * lS;
  double c3 = kx + ab * lT;
  double t1 = 1. / l1;
  double t3 = 1. / l3;
  double tkxPrime = 1. / kxPrime;
  double sfrac = fp->GetParameter(4);
  double rfrac = fp->GetParameter(5);
  double bw = fp->GetParameter(11);
  double norm = fp->GetParameter(0);
  double alpha1 = sfrac * bw * norm;
  double alpha3 = (1. - sfrac-rfrac) * bw * norm;
  double snorm = alpha1 * c1 * kx / (l1 - kxPrime) / tXe;
  double tnorm = alpha3 * c3 * kx / (l3 - kxPrime) / tXe;
  printf("\t ls %.3E c1 %.3E  lt %.3E  c3 %.3E  tMix %.3E  kxPrime %.3E  l1-kxPrime %.3E  l3-kxPrime %.3E  ab %.3E\n", lS, c1, lT, c3, tmixPar, kxPrime, l1 - kxPrime, l3 - kxPrime,ab);
  printf("\t sfrac %.3f  tfrac %.4E alpha1 %.3E alpha3 %.3E snorm %.3E tnorm %.3E \n", sfrac, 1.-sfrac-rfrac, alpha1, alpha3, snorm, tnorm);
}

modelFit::modelFit(int thefit, int ichan, double ppm)
{
  TString names[NTYPES];
  names[0] = TString("singlet");
  names[1] = TString("triplet");
  names[2] = TString("mixed");
  names[3] = TString("xenon");
  names[4] = TString("total");

  int ilevel = level(ichan);
  double dist = distanceLevel[ilevel];
  double ab = Absorption(ppm, dist);
  double kplus = 1;
  double kxPrime = 1;
  theChan = ichan;
  thePPM = ppm;

  distanceLevel[0] = 11.6;
  distanceLevel[1] = 23.2;
  distanceLevel[2] = 34.8;
  distanceLevel[3] = 36.0;

  for (int ichan = 0; ichan < 13; ++ichan)
  {
    effGeo[ichan] = effGeoFunc(ichan);
  }

  double xlow = 0.;
  double xhigh = 15000; //15000 full waveform
  fp = new TF1(Form("ModelFit-%.2f-type-%s-chan-%i", ab, names[thefit].Data(), ichan), lightModel, xlow, xhigh, NPARS);
  printf(" modelFit: set %i fit range %f to %f  \n", thefit, 0., xMax);
 
  fp->SetParName(0, "norm");
  fp->SetParName(1, "PPM");
  fp->SetParName(2, "tau3");
  fp->SetParName(3, "kp");
  fp->SetParName(4, "sfrac");
  fp->SetParName(5, "rfrac");
  fp->SetParName(6, "ab");
  fp->SetParName(7, "kxprime");
  fp->SetParName(8, "tmix");
  fp->SetParName(9, "bgk");
  fp->SetParName(10, "chan");
  fp->SetParName(11, "binw");
  fp->SetParName(12, "type");

  fp->SetParameter(0, nPhotons);
  fp->SetParLimits(0, .01 * nPhotons, 10 * nPhotons);
  fp->SetParameter(1, ppm);
  fp->SetParameter(2, tTriplet);
  fp->SetParameter(3, kplus);
  fp->SetParameter(4, 0.2);
  fp->SetParLimits(4, .01, 1.);
  fp->SetParameter(5, 1./1000.);
  fp->SetParLimits(5, 0,0.1);
  fp->SetParameter(6, ab);
  fp->SetParLimits(6, 1.E-9, 1.);
  fp->SetParameter(7, 2 * kxPrime);
  fp->SetParameter(8, tMix);
  fp->SetParameter(9, 50);

  //these are always fixed
  //fp->FixParameter(5, sFrac);
  fp->FixParameter(10, ichan);
  fp->FixParameter(11, binwidth);
  fp->FixParameter(12, thefit);

  // fp->SetParLimits(9,1.E3,20.E3);
  fp->SetTitle(Form("ModelFit-type-%s-chan-%i-%0.1f-PPM-ab-%.2f", names[thefit].Data(), ichan, ppm, ab));
  fp->SetNpx(1000); // numb points for function
  fp->Print();

  cout << "defined function " << fp->GetName() << "  " << fp->GetTitle() << endl;
}
