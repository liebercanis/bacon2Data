// MG added additions from Doug and geometry Jun 27, 2023
// time is in microseconds
// model with absorption
// static double tres = 7.86;
// static double tres = 5.4;
static double tres = 5.4;
static double tTriplet = 1600.0; // 2100.0;
static double tSinglet = 5.0;
static double tMix = 4700.;
static double tXe = 20.0;
static double kxe = 8.8E-5;
static double kplusZero = 1.3E-4;
static double xTrigger = 1406.09; // 1200.;
static double xMax = 10000.;
static double nPhotons = 50.E3 * 5.486;

enum
{
  NTYPES = 5,
  NCHAN = 13
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

static double getDistance(int ilevel)
{
  if (ilevel < 0)
    return 1.0;
  double distance[4];
  distance[0] = 11.6;
  distance[1] = 23.2;
  distance[2] = 34.8;
  distance[3] = 36.0;
  return distance[ilevel];
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

static double effGeo[13];

static double expGaus(double x, double tau)
{
  x = x + 0.33 * tres;
  double arg1 = (tres * tres / tau - 2. * x) / 2. / tau;
  double arg2 = (tres * tres / tau - x) / sqrt(2) / tres;
  double f = 0.5 * TMath::Exp(arg1) * TMath::Erfc(arg2);
  return f;
}

static double lightModel(Double_t *xx, Double_t *par)
{
  int ichan = int(par[10]);
  int ifit = int(par[8]);
  double x = xx[0] - xTrigger;
  double bw = par[0];
  double norm = par[1];
  double ppm = par[2];
  double tTrip = par[3];
  double kp = par[4] * kplusZero;
  double tmixPar = par[9];
  double lmix = 1. / tmixPar;
  // double alpha1 = bw*norm;
  // double alpha3 = (1.-par[5])/par[5]*alpha1;
  double alpha1 = par[5] * bw * norm * effGeo[ichan];
  double alpha3 = (1. - par[5]) * bw * norm * effGeo[ichan];
  double ab = par[6];
  double k1Zero = kxe * 131. / 40.;
  double kx = k1Zero * ppm;
  double kPrime = lmix + kx + kplusZero * par[7];

  double lS = 1. / tSinglet;
  double lT = 1 / tTrip;
  double l1 = 1. / tSinglet + kp + kx;
  double l3 = 1. / tTrip + kp + kx;
  double lX = 1. / tXe;
  double c1 = kx + ab * lS;
  double c3 = kx + ab * lT;
  double t1 = 1. / l1;
  double t3 = 1. / l3;
  double tkPrime = 1. / kPrime;

  // model
  double fs = (1. - ab) * alpha1 / tSinglet * expGaus(x, t1);
  double ft = (1. - ab) * alpha3 / tTrip * expGaus(x, t3);
  double fm = alpha1 * c1 / (l1 - kPrime) * (expGaus(x, tkPrime) - expGaus(x, t1)) + alpha3 * c3 / (l3 - kPrime) * (expGaus(x, tkPrime) - expGaus(x, t3));
  fm /= tmixPar;

  double x1 = c1 * kx * alpha1 / (l1 - kPrime) / tXe * ((expGaus(x, tkPrime) - expGaus(x, tXe)) / (lX - kPrime) - (expGaus(x, t1) - expGaus(x, tXe)) / (lX - l1));
  double x3 = c3 * kx * alpha3 / (l3 - kPrime) / tXe * ((expGaus(x, tkPrime) - expGaus(x, tXe)) / (lX - kPrime) - (expGaus(x, t3) - expGaus(x, tXe)) / (lX - l3));
  double fx = x1 + x3;
  

  // printf(" %E %E %E %E %E  \n"  ,fs,ft,fm,x1,x3);
  // QE eff factors
  // Tot = (SiPMQE128*SLight + SiPMQE128*TLight + SiPMQE128*ILight + SiPMQE175*XLight + SiPMQE150*MLight);
  int ilevel = level(ichan);
  double distance[4];
  distance[0] = 11.6;
  distance[1] = 23.2;
  distance[2] = 34.8;
  distance[3] = 36.0;
  double dist = distance[ilevel];
  double SiPMQ128 = QEff128(ppm, dist);
  double SiPMQE150 = 0.238;
  double SiPMQE175 = 0.238;
  double PMTQE175 = 0.38;

  fs = fs * SiPMQ128;
  ft = ft * SiPMQ128;
  fm = fm * SiPMQE150;
  
  // PMT sees only  175
  if (ichan == 5) {
    fx = fx * SiPMQE175;
    fs = 0;
    ft = 0;
    fm = 0;
  } else if (ichan == 12) {
    fx = fx * PMTQE175;
    fs = 0;
    ft = 0;
    fm = 0;
  } else { // sipms do not see 175
    fx = 0;
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

  return f;
}

class modelFit
{
public:
  // double taut = 2100.;
  modelFit(int thefit, int ichan, double ppm);
  virtual ~modelFit() { ; }
  TF1 *fp;
  enum
  {
    NPARS = 11
  };
  double binwidth = 2.0; // ns;
  double norm = 3.0E4;
  double sFrac = 0.2;
  int theChan;
  double thePPM;
  void show();
};

void modelFit::show()
{
  printf(" \n\n >>> modelFit fitted parameters fit chan %i ppm %f \n", (int)fp->GetParameter(11), fp->GetParameter(2));
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  param %i %s %.3E +/- %.3E \n", ii, fp->GetParName(ii), fp->GetParameter(ii), fp->GetParError(ii));
  }
  double tTrip = fp->GetParameter(3);
  double kp = fp->GetParameter(4) * kplusZero;
  double ppm = fp->GetParameter(2);
  double ab = fp->GetParameter(6);
  double tmixPar = fp->GetParameter(9);

  double k1Zero = kxe * 131. / 40.;
  double kx = k1Zero * ppm;
  double kPrime = 1. / tMix + kx + kplusZero * fp->GetParameter(7);
  double lS = 1. / tSinglet;
  double lT = 1 / tTrip;
  double l1 = 1. / tSinglet + kp + kx;
  double l3 = 1. / tTrip + kp + kx;
  double lX = 1. / tXe;
  double c1 = kx + ab * lS;
  double c3 = kx + ab * lT;
  double t1 = 1. / l1;
  double t3 = 1. / l3;
  double tkPrime = 1. / kPrime;
  double sfrac = fp->GetParameter(5);
  double bw = fp->GetParameter(0);
  double norm = fp->GetParameter(1);
  double alpha1 = sfrac * bw * norm;
  double alpha3 = (1. - sfrac) * bw * norm;
  double snorm = alpha1 * c1 * kx / (l1 - kPrime) / tXe;
  double tnorm = alpha3 * c3 * kx / (l3 - kPrime) / tXe;
  printf("\t ls %.3E c1 %.3E  lt %.3E  c3 %.3E  tMix %.3E  kPrime %.3E  l1-kPrime %.3E  l3-kPrime %.3E \n", lS, c1, lT, c3, tmixPar, kPrime, l1 - kPrime, l3 - kPrime);
  printf("\t S norm %.3E T norm %.3E\n", snorm, tnorm);
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
  double dist = getDistance(ilevel);
  double ab = Absorption(ppm, dist);
  double kplus = 1;
  double kPrime = 1;
  theChan = ichan;
  thePPM = ppm;
  effGeo[0] = 0.000024;
  effGeo[1] = 0.000024;
  effGeo[2] = 0.000024;
  effGeo[3] = 0.000053;
  effGeo[4] = 0.000053;
  effGeo[5] = 0.000053;
  effGeo[6] = 0.000213;
  effGeo[7] = 0.000213;
  effGeo[8] = 0.000213;
  effGeo[9] = 0.028648;
  effGeo[10] = 0.028648;
  effGeo[11] = 0.028648;
  effGeo[12] = 0.000022;

  fp = new TF1(Form("ModelFit-%.2f-type-%s-chan-%i", ab, names[thefit].Data(), ichan), lightModel, xTrigger - 2., xMax, NPARS);
  printf(" modelFit: set %i fit range %f to %f  \n", thefit, xTrigger - 10., xMax);
  fp->SetParName(0, "binw");
  fp->SetParName(1, "norm");
  fp->SetParName(2, "PPM");
  fp->SetParName(3, "tau3");
  fp->SetParName(4, "kp");
  fp->SetParName(5, "sfrac");
  fp->SetParName(6, "ab");
  fp->SetParName(7, "kprime");
  fp->SetParName(8, "type");
  fp->SetParName(9, "tmix");
  fp->SetParName(10, "chan");

  fp->FixParameter(0, binwidth);
  fp->SetParameter(1, nPhotons);
  fp->FixParameter(2, ppm);
  fp->SetParameter(3, tTriplet);
  fp->FixParameter(4, kplus);
  fp->SetParameter(5, sFrac);
  fp->SetParLimits(5, .01, .5);
  fp->FixParameter(6, ab);
  fp->SetParLimits(6, 1.E-9, 1.);
  fp->FixParameter(7, 2 * kPrime);
  fp->FixParameter(8, thefit);
  fp->FixParameter(9, tMix);
  fp->FixParameter(10, ichan);

  // fp->SetParLimits(9,1.E3,20.E3);
  fp->SetTitle(Form("ModelFit-type-%s-chan-%i-%0.1f-PPM-ab-%.2f", names[thefit].Data(), ichan, ppm, ab));
  fp->SetNpx(1000); // numb points for function
  fp->Print();

  cout << "defined function " << fp->GetName() << "  " << fp->GetTitle() << endl;
}
