#include "TString.h"
// MG added additions from Doug and geometry Jun 27, 2023
// time is in microseconds
// model with absorption
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
        SiPMQE150 = 0.238;:%:
        SiPMQE175 = 0.238;
        PMTQE175 = 0.38;
*/

static int theBinWidth = 2;
static double nPhotons = 4.1E4; // 7.6*5.4E3
static double fillFactor = 0.7;
static double nominalGain = 227.4;
static double nonZeroFraction = 1. - 0.61;
static double tTriplet0 = 1600.0; // 2100.0;
static double tSinglet0 = 7.0;
static double tauM0 = 4700.;
static double tXe0 = 20.0;
static double kUnit0 = 1.0E-4; // kq in the paper quenching rate  multiplied by binwidth
static double kq0 = 1.3;
static double kxe0 = 0.88;
static double distanceLevel[4];
static double effGeo[13];
// static double trec = 37.1; // ns Eur. Phys. J. C (2013) 73:2618
static double kdiffusion0 = 0.88;
; // this is starting value Eur. Phys. J. C (2013) 73:2618
// static double effGeo4=5.322512E-05;
static double landauSigma0 = 6.0;
static double gausSigma0 = 5.0;
static double landauFrac0 = 0.12;
static double resolution0 = 5.0; // hit time resolution
static int nominalTriggerBin = 688;
static double nominalTriggerTime = 2. * double(nominalTriggerBin);
static double MPV0 = 2. * 698.0;
static double sfrac0 = 0.8; // nominal singlet frac
// times in ns
static int startBin = 660;
static int singletEndBin = 725;
static int endBin = 7500;
static double singletEndTime = 2. * double(singletEndBin);
static double endTime = 2. * double(endBin);
static double startTime = 2. * double(startBin); // hWave->GetBinLowEdge(maxBin) + hWave->GetBinWidth(maxBin) / 2.;
static double fitStart = 1500.;                  // singletEndTime
static double fitEnd = endTime;                  // fit range end
double massRatio = 131. / 40.;                   // ratio of Xe/Ar molecular weight
static double SiPMQE128Ham = 0.15;
static double SiPMQE150 = 0.238;
static double SiPMQE175 = 0.238;
// static double PMTQE175 = 0.38;
static double PMTQE150 = 0.01;
static double PMTQE175 = 0.38;
static double PMTQE400 = 0.35;
double SiPMQ128 = SiPMQE128Ham;

enum
{
  NTYPES = 5,
  NCHAN = 13
};

// fit parameters
enum
{
  NORM = 0,
  SFRAC,
  PPM,
  MPV,
  GAUSSIGMA,
  LANDAUSIGMA,
  LANDAUFRAC,
  TAU3,
  KQ,
  AB,
  KDIFFUSION,
  TAUM,
  BKGCONST,
  BKGTAU,
  CHAN,
  TYPE,
  NPARS
};

double effGeoFunc(int ichan)
{
  int ilevel = -1;
  if (ichan == 6 || ichan == 7 || ichan == 8 || ichan == 9 || ichan == 10 || ichan == 11)
    ilevel = 0;
  else if (ichan == 3 || ichan == 4 || ichan == 5)
    ilevel = 1;
  else if (ichan == 0 || ichan == 1 || ichan == 2)
    ilevel = 2;
  else if (ichan == 12)
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
  double aPmt = TMath::Pi() / 4.0 * pow(6.4, 2); // R11410-20  Effective area : 64 mm dia
  double a = pow(0.6, 2.);
  if (ichan == 12)
    a = aPmt;

  double b = 4.0 * TMath::Pi();
  double distance2[4];
  distance2[0] = pow(distanceLevel[0], 2.);
  distance2[1] = pow(distanceLevel[1], 2.);
  distance2[2] = pow(distanceLevel[2], 2.);
  distance2[3] = pow(distanceLevel[3], 2.);

  e = a / b / distance2[ilevel];
  return e;
}

static double singletPeak(double *xx, double *par)
{
  int ichan = int(par[0]);
  double norm = par[1];
  double mpv = par[2];
  double lsigma = par[3];
  double gsigma = par[4];
  double lfrac = par[5];
  double t = xx[0];
  double g = TMath::Landau(t, mpv, lsigma);
  if (ichan == 4)
    return norm * lfrac * g;
  double f = TMath::Gaus(t, mpv, gsigma, 1); // normalized gaus
  return norm * (1. - lfrac) * f + norm * lfrac * g;
}

static int level(int ichan)
{
  int ilevel = -1;
  if (ichan > 5)
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
  double arg1 = (resolution0 * resolution0 / tau - 2. * x) / 2. / tau;
  double arg2 = (resolution0 * resolution0 / tau - x) / sqrt(2) / resolution0;
  double f = 0.5 * TMath::Exp(arg1) * TMath::Erfc(arg2);
  return f;
}

static double lightModel(Double_t *xx, Double_t *par)
{
  double f = 0; // function return value
  int ichan = int(par[CHAN]);
  int ifit = int(par[TYPE]);
  double binwidth = theBinWidth;         // 2 ns
  double x = xx[0] - nominalTriggerTime; // xx[0]  is bin center
  // if (xx[0] < 1376)
  //   return 0;
  //      printf("xx[0] %f .............. ",xx[0]);
  double mpv = par[MPV] - nominalTriggerTime;
  double norm = par[NORM];
  double sfrac = par[SFRAC];
  double ppm = par[PPM];
  // *** *************** do not scale all lifetimes by binwidth of 2 ns
  double tTriplet = par[TAU3];
  double tSinglet = tSinglet0;
  double tauM = par[TAUM]; // fp->SetParName(8, "tauM");
  double tXe = tXe0;
  double lsigma = par[LANDAUSIGMA];
  double gsigma = par[GAUSSIGMA];
  // **** scale kUnit multiply by binwidth
  double kUnit = kUnit0;
  // *******************************************
  double kq = par[KQ] * kUnit; // fp->SetParName(3, "kq");
  double lmix = 1. / tauM;
  double bkg = par[BKGCONST] * TMath::Exp(-x / par[BKGTAU]);
  double lfrac = par[LANDAUFRAC];
  double transmission = 1. - par[AB];
  double alpha1 = binwidth * norm * sfrac * effGeo[ichan] * fillFactor * transmission * SiPMQ128;
  // printf(" %f %f %f %f %f .....", sfrac, bw, norm, effGeo[ichan],alpha1);
  double alpha3 = binwidth * norm * (1. - sfrac) * effGeo[ichan] * fillFactor * transmission * SiPMQ128;
  double k1Zero = kUnit * massRatio * par[KDIFFUSION]; // p->SetParName(7, "kdiffusion";
  double kx = k1Zero * ppm;
  // double kxPrime = lmix + kx + kUnit * par[7]; //fp->SetParName(7, "kx");
  double kxPrime = lmix + kx + kq; // fp->SetParName(7, "kx");

  double lS = 1. / tSinglet;
  double lT = 1 / tTriplet;
  double l1 = 1. / tSinglet + kq + kx;
  double l3 = 1. / tTriplet + kq + kx;
  double lX = 1. / tXe;
  double c1 = kx + par[AB] * lS;
  double c3 = kx + par[AB] * lT;
  double t1 = 1. / l1;
  double t3 = 1. / l3;
  double tkxPrime = 1. / kxPrime;

  // singlet parameters defined from start of waveform
  // double fs = (1. - ab) * alpha1 / tSinglet * expGaus(x, t1);
  double fs = 0;

  double gsinglet = TMath::Landau(x, mpv, lsigma, 1); // normalized
  double fsinglet = TMath::Gaus(x, mpv, gsigma, 1);   // normalized gaus
  if (ichan == 4 || ichan == 1)
    fs = alpha1 * lfrac * gsinglet;
  else
    fs = alpha1 * ((1 - lfrac) * gsinglet + lfrac * gsinglet);
  // add in sipm efficiency
  // printf("%f alpha1 %E alpha3 %E \n", x, alpha1, alpha3);
  // fs = max(fs, 1.E-20);
  if (isnan(fs))
  {
    printf("NAN at x=%f \n", x);
    fs = 0;
  }

  // return if just fitting singlet
  // if (ifit == 0)
  //  return fs;
  // end of singlet

  double ft = 0;
  double fm = 0;
  double fx = 0;

  /* divide by lifetimes to convert into light yield */
  if (x > 0)
  {
    ft = alpha3 / tTriplet * expGaus(x, t3);
    fm = alpha1 * c1 / (l1 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t1)) + alpha3 * c3 / (l3 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t3));
    fm /= tauM;

    // printf("  %E %E  %E  \n",fs,x,expGaus(x,t1));
    double x1 = c1 * kx * alpha1 / (l1 - kxPrime) / tXe * ((expGaus(x, tkxPrime) - expGaus(x, tXe)) / (lX - kxPrime) - (expGaus(x, t1) - expGaus(x, tXe)) / (lX - l1));
    double x3 = c3 * kx * alpha3 / (l3 - kxPrime) / tXe * ((expGaus(x, tkxPrime) - expGaus(x, tXe)) / (lX - kxPrime) - (expGaus(x, t3) - expGaus(x, tXe)) / (lX - l3));

    /*
       x1 = c1 * kx * alpha1 / (l1 - kxPrime) / tXe * ((TMath::Exp(-x/tkxPrime) - TMath::Exp(-x/tXe)) / (lX - kxPrime) - (TMath::Exp(-x/t1) - TMath::Exp(-x/tXe)) / (lX - l1));
      x3 = c3 * kx * alpha3 / (l3 - kxPrime) / tXe * ((TMath::Exp(-x/tkxPrime) - TMath::Exp(-x/tXe)) / (lX - kxPrime) - (TMath::Exp(-x/t3) - TMath::Exp(-x/tXe)) / (lX - l3));
      */

    fx = x1 + x3;
    // printf(" f1 %E f2 %E x1 %E x3 %E fx = %E expa %E expb %E t1 %E t2 %E \n", (TMath::Exp(-x/tkxPrime) - TMath::Exp(-x/tXe)) / (lX - kxPrime) ,  (TMath::Exp(-x/t3) - TMath::Exp(-x/tXe)) / (lX - l1)  , x1, x3 , fx,TMath::Exp(-x/t3), TMath::Exp(-x/tXe),x/t3,x/tXe);
  }

  // QE eff factors
  // Tot = (SiPMQE128*SLight + SiPMQE128*TLight + SiPMQE128*ILight + SiPMQE175*XLight + SiPMQE150*MLight);
  int ilevel = level(ichan);
  double dist = distanceLevel[ilevel];
  double SiPMQ128 = QEff128(ppm, dist);
  double SiPMQE150 = 0.238;
  double SiPMQE175 = 0.238;
  double PMTQE175 = 0.38;

  // ft = ft;
  fm = fm * SiPMQE150 / SiPMQ128;

  // PMT sees only  175
  if (ichan == 5)
  {
    fx = fx * SiPMQE175 / SiPMQ128;
    fs = 0;
    ft = 0;
    fm = 0;
  }
  else if (ichan == 12)
  {
    fx = fx * PMTQE175 / SiPMQ128;
    fs = 0;
    ft = 0;
    fm = 0;
  }
  else
  { // sipms do not see 175
    fx = fx * SiPMQE175 / SiPMQ128;
  }
  if (ifit == 0)
    f = fs;
  else if (ifit == 1)
    f = ft;
  else if (ifit == 2)
    f = fm;
  else if (ifit == 3)
    f = fx;
  else if (ifit == 4)
    f = fs + ft + fx + fm + bkg; // return all
  // f = f + bkg;
  // printf(" time = %f (%f) bkg %E %E %E %E tot %E \n"  ,xx[0],x,bkg,fs,ft,fm,f);

  return f;
}

class modelFit
{
public:
  // double taut = 2100.;
  modelFit(int theFit, int ichan, double ppm);
  virtual ~modelFit() { ; }
  TF1 *fp;
  double binwidth = 2.0; // ns;
  double norm = 3.0E4;
  double sFrac = 0.2;
  int theChan;
  double thePPM;
  void show();
  void showEff(int ichan);
  void showEff();
};

void modelFit::showEff(int ichan)
{
  cout << " modelFit::showEff " << endl;
  int ilevel = level(ichan);
  double d = distanceLevel[ilevel];
  double e = effGeo[ilevel];
  printf("chan %i level %i dist %.3f e %.3E \n", ichan, ilevel, d, e);
  // printf("effGeo[%i]=%f ; \n", ichan, e);
}
void modelFit::showEff()
{
  cout << " modelFit::showEff " << endl;
  for (int ichan = 0; ichan < 13; ++ichan)
  {
    int ilevel = level(ichan);
    double d = distanceLevel[ilevel];
    double e = effGeo[ilevel];
    printf("chan %i level %i dist %.3f e %.3E \n", ichan, ilevel, d, e);
    // printf("effGeo[%i]=%f ; \n", ichan, e);
  }
}

void modelFit::show()
{
  int ichan = fp->GetParameter(CHAN);
  double k1Zero = kUnit0 * fp->GetParameter(KDIFFUSION) * massRatio;
  int ilevel = level(ichan);
  printf(">>> modelFit static starting parameters chan %i \n", ichan);
  printf("\t level %i distance %.2f geometric efficiency %.2E \n", ilevel, distanceLevel[ilevel], effGeoFunc(ichan));
  printf("\t norm photon yield %.3E \n", nPhotons);                 // nominalTriggerTime
  printf("\t nominal trigger time %.2f ns \n", nominalTriggerTime); // nominalTriggerTime
  printf("\t Singlet Landau sigma  %.2f ns \n", landauSigma0);      // from fit
  printf("\t Singlet Landau fraction  %.2f ns \n", landauFrac0);    // from fit
  printf("\t Singlet Gausian  %.2f ns \n", gausSigma0);             // from fit
  printf("\t triplet lifetime %.2f ns\n", tTriplet0);
  printf("\t singlet lifetime %.2f ns\n", tSinglet0);
  printf("\t xenon lifetime %.2f ns\n", tXe0);
  printf("\t mixed state  lifetime %f ns\n", tauM0);
  printf("\t kUnit %.2E ns-1\n", kUnit0);
  printf("\t kq %.2f ns-1 \n", kq0);
  printf("\t kx  %.2f ns-1 \n", kxe0);
  printf("\t k1Zero %.2E ns-1 \n", k1Zero); // ref[12]
  printf("\t hit time resolution %.2f ns\n", resolution0);
  printf("\t fit range  %.0f  %.0f \n", fitStart, fitEnd); // nominalTriggerTime
  // J. Calvo, et al., Measurement of the attenuation length of argon scin- tillation light in the ArDM LAr TPC,

  printf(" \n >>> initial modelFit fitted parameters  \n");
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  %i & %s & %.2E  \\\\\n", ii, fp->GetParName(ii), fp->GetParameter(ii));
  }
}

modelFit::modelFit(int theFit, int ichan, double ppm)
{
  TString names[NTYPES];
  names[0] = TString("singlet");
  names[1] = TString("triplet");
  names[2] = TString("mixed");
  names[3] = TString("xenon");
  names[4] = TString("total");

  std::vector<TString> vparNames;
  vparNames.resize(int(NPARS));
  vparNames[NORM] = TString("NORM");
  vparNames[TAU3] = TString("TAU3 ns");
  vparNames[PPM] = TString("PPM");
  vparNames[MPV] = TString("MPV");
  vparNames[LANDAUSIGMA] = TString("LandauSigma");
  vparNames[GAUSSIGMA] = TString("GausSigma");
  vparNames[LANDAUFRAC] = TString("LandauFrac");
  vparNames[KQ] = TString("KQ ns-1");
  vparNames[SFRAC] = TString("SFRAC");
  vparNames[AB] = TString("AB");
  vparNames[KDIFFUSION] = TString("kdiffusion");
  vparNames[TAUM] = TString("taum ns");
  vparNames[BKGCONST] = TString("bkgconst");
  vparNames[BKGTAU] = TString("bkgtau");
  vparNames[CHAN] = TString("chan");
  vparNames[TYPE] = TString("type");

  int ilevel = level(ichan);
  double dist = distanceLevel[ilevel];
  double ab = Absorption(ppm, dist);

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

  fp = new TF1(Form("ModelFit-%.2f-type-%s-chan-%i", ab, names[theFit].Data(), ichan), lightModel, fitStart, fitEnd, NPARS);
  printf(" modelFit: set %i fit range %f to %f  \n", theFit, fitStart, fitEnd);

  for (int ipar = 0; ipar < NPARS; ++ipar)
    fp->SetParName(ipar, vparNames[ipar].Data());

  // default initial values and limits
  fp->SetParameter(NORM, nPhotons);
  fp->SetParameter(SFRAC, sfrac0);
  fp->SetParameter(PPM, ppm);
  fp->SetParameter(MPV, MPV0);
  fp->SetParameter(GAUSSIGMA, gausSigma0);
  fp->SetParameter(LANDAUSIGMA, landauSigma0);
  fp->SetParameter(LANDAUFRAC, landauFrac0);
  fp->SetParameter(TAU3, tTriplet0);
  fp->SetParameter(KQ, kq0);
  fp->SetParameter(AB, ab);
  fp->SetParameter(KDIFFUSION, kdiffusion0); // this is starting value Eur. Phys. J. C (2013) 73:2618
  fp->SetParameter(TAUM, tauM0);
  fp->SetParameter(CHAN, ichan);
  fp->SetParameter(TYPE, theFit);

  // par limits
  // fp->SetParLimits(NORM1, 1.E-4 * nPhotons, 10 * nPhotons);
  fp->SetParLimits(KDIFFUSION, 0, 0.1);
  fp->SetParLimits(LANDAUSIGMA, 1.0E-5 * landauSigma0, 10. * landauSigma0);
  fp->SetParLimits(GAUSSIGMA, 1.0E-5 * gausSigma0, 10. * gausSigma0);
  fp->SetParLimits(LANDAUFRAC, 0., 1.);
  if (ichan == 1 || ichan == 4)
  {
    fp->FixParameter(GAUSSIGMA, gausSigma0); // not fitting
    fp->FixParameter(LANDAUFRAC, 1.0);       // all Landau
  }

  // these are always fixed
  // fp->FixParameter(5, sFrac);
  fp->FixParameter(CHAN, ichan);
  fp->FixParameter(TYPE, theFit);

  // fp->SetParLimits(9,1.E3,20.E3);
  fp->SetTitle(Form("ModelFit-type-%s-chan-%i-%0.1f-PPM-ab-%.2f", names[theFit].Data(), ichan, ppm, ab));

  fp->SetNpx(1000); // numb points for function
  fp->Print();

  cout << "defined function " << fp->GetName() << "  " << fp->GetTitle() << endl;
}
