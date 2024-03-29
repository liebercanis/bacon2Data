// MG added additions from Doug and geometry Jun 27, 2023
// time is in microseconds
// model with absorption
// static double tres = 7.86;
// static double tres = 5.4;
#include "TString.h"
#include "TF1.h"
#include <TNtuple.h>

enum
{ 
  NCOMP =7,
  NCHAN = 13,
  NPARS = 10 // add another parameter for mixed state decay
  // add another parameter for recombination time.
};

TNtuple* ntScan = new TNtuple("ntScan","ntScan","ppm:fx:f");
static double tres = 10;     // maybe 13 for PMT
static double tresPmt = 13;  // maybe 13 for PMT
static double tTriplet = 1600.0; // 2100.0;
static double tSinglet = 5.0;
static double tMix = 4700.;
static double tXe = 20.0;
static double kxe = 8.8E-5;             // diffusion time
static double kqZero = 1.3E-4;          // kq in the paper quenching rate
static double xTrigger = 1411.;         // 1200.;
static double nPhotons = 50.E3 * 5.486; // 274300.00
static double distanceLevel[4];
static double effGeo[13];
/* Ion-beam excitation of liquid argon M. Hofmann et al.  Eur. Phys. J. C (2013) 73:2618 */
static double trecon = 37.1; // ns   
/*
Electron transport and electron–ion recombination in liquid argon simulation based on the Cohen–Lekner theory 
  Mariusz Wojcik a, Tachiya b doi.org/10.1016/S0009-2614(02)01177-6
*/
static double trecomb = 2.976081E-03; //ns
static double foffset = 0;
static double buff[13][7500];
static double lpar[10];       // pass parameters to light model
static TString lparNames[10]; // pass parameters to light model

static double SiPMQE128Ham = 0.15;
static double SiPMQE150 = 0.238;
static double SiPMQE175 = 0.238;
//static double PMTQE175 = 0.38;
static double PMTQE150 = 0.01;
static double PMTQE175 = 0.38;
static double PMTQE400 = 0.35;
static double background[13];
TF1 *ffit[13];
TF1 *fmodel[13];
TF1 *ffitPmt[NCOMP];
TF1 *ffitChan[NCOMP];
TString compName[NCOMP];
static double xlow = 900;
static double xhigh = 12000;


static bool goodChannel(int ic)
{
  bool val = true;
  if (ic == 5 || ic == 6 || ic==8 || ic == 3 || ic == 9 || ic==10 || ic==11)
    val = false;
  return val;
}
static void setParNames() // tousif
{
  lparNames[0] = TString("norm");
  lparNames[1] = TString("ppm");
  lparNames[2] = TString("tau3");
  lparNames[3] = TString("kp");
  lparNames[4] = TString("sfrac");
  lparNames[5] = TString("rfrac");
  lparNames[6] = TString("kxprime");
  lparNames[7] = TString("tmix");
  lparNames[8] = TString("taumix");
  lparNames[9] = TString("taurecomb");
}

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
*/

static double effGeoFunc(int ichan)
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

  double e = 1.0;
  if (ilevel < 0)
    return e;
  /*
  Area of SiPMs is 6.0mm x 6.0mm

      Channels 6, 7, and 8 are at 11.6 cm
      from the source Channels 3, 4, and 5 are at 23.2 cm
      from the source Channels 0, 1, and 2 are at 34.8 cm from the source Channel 12 is at 36 cm from the source.
      */
  double aPmt = TMath::Pi() / 4.0 * pow(6.4, 2); // R11410-20  Effective area : 64 mm dia unit is cm!
  double a = pow(0.6, 2.);
  if(ichan==12)
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

static int level(int ichan)
{
  distanceLevel[0] = 11.6;
  distanceLevel[1] = 23.2;
  distanceLevel[2] = 34.8;
  distanceLevel[3] = 36.0;
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
  // double val =  0.16 - (10 - ppm * dist) * 0.002;
  return SiPMQE128Ham;
}

static double Absorbtion(double ppm, double dist)
{
  // Calculate absorption as a function of distance and xenon concentration.%
  // Taken from fits to Neumeier data at 0.1 PPM and scaled;
  double A = 0.615;
  ppm = max(1.0E-9, ppm);
  double lambda1 = 12.7 * 0.1 / ppm;
  double lambda2 = 740 * 0.1 / ppm;
  double Tr128 = A * exp(-dist / lambda1) + (1 - A) * exp(-dist / lambda2);
  return 1. - Tr128;
}

static double expGaus(double x, double tau)
{
  x += 7.2; // compensate for shift in mean due to smearing of 10 percent
  double arg1 = (tres * tres / tau - 2. * x) / 2. / tau;
  double arg2 = (tres * tres / tau - x) / sqrt(2) / tres;
  double f = 0.5 * TMath::Exp(arg1) * TMath::Erfc(arg2);
  return f;
}

/*
npar: number of currently variable parameters
par: array of (constant and variable) parameters
flag: Indicates what is to be calculated (see example below)
grad: array of gradients Output parameters:
fval: The calculated function value.
grad: The (optional) vector of first derivatives).
*/

/* parameter definitions
  fp->setparname(0, "norm");
  fp->SetParName(1,"ppm");
  fp->SetParName(2, "tau3");
  fp->SetParName(3, "kp");
  fp->SetParName(4, "sfrac");
  fp->SetParName(5, "rfrac");
  fp->SetParName(6, "kxprime");
  fp->SetParName(7, "tmix");
  fp->SetParName(8, "bkg");
*/

/* returns fit light yield  for given model paraters */
static double model(int ichan,int ifit,  double xbin, double ab, double SiPMQ128, double effGeo)
{
  double bw = 2.;
  double x = xbin - xTrigger;
  double norm = lpar[0];
  double ppm = lpar[1];
  double tTrip = lpar[2];
  double kp = lpar[3];
  double tmixPar = lpar[7];
  double lmix = 1. / tmixPar;
  double bkg = background[ichan];

  // double alpha1 = bw*norm;
  // double alpha3 = (1.-lpar[5])/lpar[5]*alpha1;
  double sfrac = lpar[4];
  double rfrac = lpar[5];
  double alpha1 = sfrac * bw * norm * effGeo;
  double alpha3 = (1. - sfrac) * bw * norm * effGeo;
  //double alpha3 = (1. - sfrac - rfrac) * bw * norm * effGeo;
  double k1Zero = kxe * 131. / 40.;
  double kx = k1Zero * ppm;
  double kxPrime = lmix + kx + lpar[6];
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
  // printf(" fs %E ab %f alpha1 %f tsinglet %f expG %E sipm %f \n", fs, ab, alpha1, tSinglet, expGaus(x, t1), SiPMQ128);

  double ft = 0;
  double fm = 0;
  double fx = 0;
  double frec = 0;
  //double Nrec = norm / 1000.;
  /* recombination function fit to Landau
  1  Constant     8.75979e+02   1.30528e+01   9.10733e-02  -2.82338e-06
     2  MPV          1.41693e+03   2.90900e-01   2.31557e-03   1.82918e-04
     3  Sigma        1.44940e+01   1.74659e-01   5.88558e-06  -1.15829e-01
      Landau(Double_t x, Double_t mpv = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
  */
  frec = rfrac * bw * norm * effGeo * TMath::Landau(x,0, lpar[9]);
  //pow(1 + x / trecon, -2.);
  // frecomb = bw * Nrec * TMath::Exp(-x/trecomb);
  ft = (1. - ab) * alpha3 / tTrip * expGaus(x, t3);
  double fmnorm = alpha1 * c1 / (l1 - kxPrime);
  fmnorm = max(0., fmnorm);
  fm = fmnorm * (expGaus(x, tkxPrime) - expGaus(x, t1)) + alpha3 * c3 / (l3 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t3));
  fm /= tMix;

  double fxnorm1 = c1 * kx * alpha1 / (l1 - kxPrime) / tXe;
  double fxnorm3 = c3 * kx * alpha3 / (l3 - kxPrime) / tXe;
  fxnorm1 = max(0., fxnorm1);
  fxnorm3 = max(0., fxnorm3);
  if (fxnorm1 < 0)
    printf(" negative form1!!! fx = %f %f c1 %f c3 %f fxnorm1 = %E fxnorm3 %E  \n ", fx, alpha1, c1, c3, fxnorm1, fxnorm3);

  double x1 = fxnorm1 * ((expGaus(x, tkxPrime) - expGaus(x, tXe)) / (lX - kxPrime) - (expGaus(x, t1) - expGaus(x, tXe)) / (lX - l1));
  double x3 = fxnorm3 * ((expGaus(x, tkxPrime) - expGaus(x, tXe)) / (lX - kxPrime) - (expGaus(x, t3) - expGaus(x, tXe)) / (lX - l3));

  fx = x1 + x3;

  /*
  double x1 = c1 * kx * alpha1 / (l1 - kxPrime) / tXe * ((exp(-x/tkxPrime) - exp(-x/tXe)) / (lX - kxPrime) - (exp(-x/t1) - exp(-x/tXe)) / (lX - l1));
  double x3 = c3 * kx * alpha3 / (l3 - kxPrime) / tXe * ((exp(-x/tkxPrime) - exp(-x/tXe)) / (lX - kxPrime) - (exp(-x/t3) - exp(-x/tXe)) / (lX - l3));
  */

  fs = fs * SiPMQ128;
  ft = ft * SiPMQ128;
  fm = fm * SiPMQE150;

  //  glass covered  175
  if (ichan == 5)
  {
    fx = fx * SiPMQE175;
    fs = 0;
    ft = 0;
    fm = 0;
    frec = frec * SiPMQE175;
  }
  // PMT sees only  175
  else if (ichan == 12)
  {
    fx = fx * PMTQE175;
    fs = 0;
    ft = 0;
    fm = fm*PMTQE150;
    frec = frec * PMTQE400;
  }
  else
  { // sipms do not see 175
    fx = fx * SiPMQE175;
    frec = frec * SiPMQE175;
  }
  double f = fs + frec + ft + fx + fm + bkg+foffset;
  //double f = fs + frec + ft + fx + fm;
  double fcomp[6];
  fcomp[0] = f;
  fcomp[1] = fs;
  fcomp[2] = ft;
  fcomp[3] = fx;
  fcomp[4] = fm;
  fcomp[5] = frec;
  //  leave these diognostics in
  if (1)
  {
    if (isnan(fs))
      printf("xxx fs is NAN x = %f\n", x);
    if (isnan(ft))
      printf("xxx ft is NAN x = %f\n", x);
    if (isnan(fm))
      printf("xxx fm is NAN x = %f\n", x);
    if (isnan(fx))
      printf("xxx fx is NAN x = %f\n", x);
    if (isnan(bkg))
      printf("xxx bkg is NAN x = %f\n", x);
    if (isnan(frec))
      printf("xxx frec is NAN x = %f\n", x);
    if (isnan(foffset))
      printf("xxx foffset is NAN x = %f\n", x);
    if (isnan(f))
      printf("xxxxxxxxxxx f is NAN x = %f\n", x);
  }

  return fcomp[ifit];
}

// for plotting st the end
double modelFunc(double *xx, double *par)
{
  int ic = par[0]; // channel
  int ifit = par[1]; // fit component e.g. singelet , triplet ... 
  double ppm = lpar[1];
  double dist = distanceLevel[level(ic)];
  double SiPMQ128 = QEff128(ppm, dist);
  double ab = Absorbtion(ppm, dist);
  double effGeo = effGeoFunc(ic);
  return model(ic,ifit, xx[0], ab, SiPMQ128, effGeo);
}

void createFunctions()
{
  compName[0] = TString("fTot");
  compName[1] = TString("fSinglet");
  compName[2] = TString("fTriplet");
  compName[3] = TString("fXe");
  compName[4] = TString("fMixed");
  compName[5] = TString("fRecomb");
  compName[6] = TString("fBackground");
  
  /* 
    make light yield functions-- return yield at time 
    these functions have 2 parameters.  
    Parameter 0 is the channel number
    Parameter 1 is the fit component  
  */
  for (int ifit = 0; ifit < 13; ++ifit)  // we  have 13 channels
  {
    // for fitted  parameters
    ffit[ifit] = new TF1(Form("FitToModelChan%i", ifit), modelFunc, xlow, xhigh, 2);
    ffit[ifit]->SetParameter(0, ifit);
    ffit[ifit]->SetParameter(1, 0);
    ffit[ifit]->SetNumberFitPoints(17500);
    ffit[ifit]->GetYaxis()->SetTitle("yield");
    ffit[ifit]->GetXaxis()->SetTitle("time [ns]");
    // for default parameters
    fmodel[ifit] = new TF1(Form("ModelChan%i", ifit), modelFunc, xlow, xhigh, 2);
    fmodel[ifit]->SetParameter(0, ifit);
    fmodel[ifit]->SetParameter(1, 0);
    fmodel[ifit]->SetNumberFitPoints(17500);
    fmodel[ifit]->GetYaxis()->SetTitle("yield");
    fmodel[ifit]->GetXaxis()->SetTitle("time [ns]");
  }
  // make PMT functions by comp
  for (int ifit = 0; ifit < 6; ++ifit) // we have 6 components of light
  {
    ffitPmt[ifit] = new TF1(Form("FitToModelPmtComp%s", compName[ifit].Data()), modelFunc, xlow, xhigh, 2);
    ffitPmt[ifit]->SetParameter(0, 12);
    ffitPmt[ifit]->SetParameter(1, ifit);
    ffitPmt[ifit]->SetNumberFitPoints(17500);
    ffitPmt[ifit]->GetYaxis()->SetTitle("yield");
    ffitPmt[ifit]->GetXaxis()->SetTitle("time [ns]");
  }
  // make chan 7  functions by comp
  for (int ifit = 0; ifit < 6; ++ifit)
  {
    ffitChan[ifit] = new TF1(Form("FitToModelChan7Comp%s",compName[ifit].Data()), modelFunc, xlow, xhigh, 2);
    ffitChan[ifit]->SetParameter(0, 7);
    ffitChan[ifit]->SetParameter(1, ifit);
    ffitChan[ifit]->SetNumberFitPoints(17500);
    ffitChan[ifit]->GetYaxis()->SetTitle("yield");
    ffitChan[ifit]->GetXaxis()->SetTitle("time [ns]");
  }
}

void showChannel(int ichannel)
{

  double effGeo = effGeoFunc(ichannel);
  double dist = distanceLevel[level(8)];
  double ppm = lpar[1];
  double ab = Absorbtion(ppm, dist);
  double bw = 2.;
  double norm = lpar[0];
  double tTrip = lpar[2];
  double kp = lpar[3];
  double tmixPar = lpar[7];
  double lmix = 1. / tmixPar;
  double bkg = background[8];

  // double alpha1 = bw*norm;
  // double alpha3 = (1.-lpar[5])/lpar[5]*alpha1;
  double sfrac = lpar[4];
  double rfrac = lpar[5];
  double alpha1 = sfrac * bw * norm * effGeo;
  double frec = rfrac * bw * norm * effGeo * pow(1 + 2, -2.);
  double alpha3 = (1. - sfrac - rfrac) * bw * norm * effGeo;
  double k1Zero = kxe * 131. / 40.;
  double kx = k1Zero * ppm;
  double kxPrime = lmix + kx +  lpar[7];
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

  double snorm = alpha1 * c1 * kx / (l1 - kxPrime) / tXe;
  double tnorm = alpha3 * c3 * kx / (l3 - kxPrime) / tXe;

  double total = ffit[ichannel]->Integral(xlow, xhigh);

  printf("%i total %E \n", ichannel, total);

  /*
  printf("\t ls %.3E c1 %.3E  lt %.3E  c3 %.3E  tMix %.3E  kxPrime %.3E  l1-kxPrime %.3E  l3-kxPrime %.3E  ab %.3E\n", lS, c1, lT, c3, tmixPar, kxPrime, l1 - kxPrime, l3 - kxPrime, ab);
  printf("\t sfrac %.3f  tfrac %.4E alpha1 %.3E alpha3 %.3E snorm %.3E tnorm %.3E \n", sfrac, 1. - sfrac - rfrac, alpha1, alpha3, snorm, tnorm);
  */
}

void show()
{
  setParNames();
  printf(" \n\n >>> modelFit fitted parameters fit ppm %f  \n", lpar[1]);
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("\t  param %i %s %.4E  \n", ii, lparNames[ii].Data(), lpar[ii]);
  }
  showChannel(8);
  showChannel(7);
  showChannel(4);
  showChannel(2);
  showChannel(1);
  showChannel(0);
}

/* 
fcn is required by Minuit to have exactly these argements 
returns likelihood value for some set of parameters
*/
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // pack parameters into static array onto lightModel
  // for (int k = 0; 8 < npar; ++k)
  //  lpar[k] = par[k];

  f = 0;
  distanceLevel[0] = 11.6;
  distanceLevel[1] = 23.2;
  distanceLevel[2] = 34.8;
  distanceLevel[3] = 36.0;
  double ppm = par[1];
  // loop over channels
  for (int ic = 0; ic < 13; ++ic)
  {
    if (ic == 5 || ic == 6 || ic == 8 || ic == 3 || ic == 9 || ic ==10 || ic==11 )
        continue;

    // Fit only the PMT
    if(ic!=12) continue;

    // level
    int ilevel = -1;
    if (ic == 6 || ic == 7 || ic == 8)
      ilevel = 0;
    else if (ic == 3 || ic == 4 || ic == 5)
      ilevel = 1;
    else if (ic == 0 || ic == 1 || ic == 2)
      ilevel = 2;
    else if (ic == 12)
      ilevel = 3;
    /** dist */
    double dist = distanceLevel[ilevel];
    // SiPMQ128
    double SiPMQ128 = SiPMQE128Ham;
    // double ab = Absorbtion(ppm, dist);
    double lambda1 = 12.7 * 0.1 / ppm;
    double lambda2 = 740 * 0.1 / ppm;
    double Tr128 = 0.615 * exp(-dist / lambda1) + (1 - 0.615) * exp(-dist / lambda2);
    double ab = 1. - Tr128;
    // double effGeo = effGeoFunc(ic);
    double fourPi = 12.566371;
    double effGeo = pow(0.6, 2.) / fourPi / pow(dist, 2.);
    double aPmt = TMath::Pi() / 4.0 * pow(6.40, 2); // R11410-20  Effective area : 64 mm dia units here are cm
    if (ic == 12)
      effGeo = aPmt/ fourPi/ pow(dist, 2.);

    // loop over bins
    // for (int j = 0; j < 7500; ++j) // 7500 is total
    int ilow = int(xlow / 2.);
    int ihigh = int(xhigh / 2.);
    for (int j = ilow; j < ihigh; ++j) // 7500 is total
    {
      double xbin = 2. * (double(j) + 0.5); // bin center convert to ns
      // double x = model(ic, xbin, ab, SiPMQ128, effGeo); // expected
      /* *******/
      double bw = 2.;
      double x = xbin - xTrigger;
      double norm = par[0];
      double ppm = par[1];
      double tTrip = par[2];
      double kp = par[3] ;
      double tmixPar = par[7];
      double lmix = 1. / tmixPar;
      double bkg = background[ic];

      // double alpha1 = bw*norm;
      // double alpha3 = (1.-par[5])/par[5]*alpha1;
      double sfrac = par[4];
      double rfrac = par[5];
      double alpha1 = sfrac * bw * norm * effGeo;
      //double alpha3 = (1. - sfrac - rfrac) * bw * norm * effGeo;
      double alpha3 = (1. - sfrac ) * bw * norm * effGeo;
      double k1Zero = kxe * 131. / 40.;
      double kx = k1Zero * ppm;
      double kxPrime = lmix + kx +  par[6];
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

      double ft = 0;
      double fm = 0;
      double fx = 0;

      double frec = 0;
      frec = rfrac * bw * norm * effGeo * TMath::Landau(x, 0, par[9]);
      ft = (1. - ab) * alpha3 / tTrip * expGaus(x, t3);
      double fmnorm = alpha1 * c1 / (l1 - kxPrime);
      fmnorm = max(0., fmnorm);
      fm = fmnorm * (expGaus(x, tkxPrime) - expGaus(x, t1)) + alpha3 * c3 / (l3 - kxPrime) * (expGaus(x, tkxPrime) - expGaus(x, t3));
      fm /= tMix;

      double fxnorm1 = c1 * kx * alpha1 / (l1 - kxPrime) / tXe;
      double fxnorm3 = c3 * kx * alpha3 / (l3 - kxPrime) / tXe;
      fxnorm1 = max(0., fxnorm1);
      fxnorm3 = max(0., fxnorm3);

      double x1 = fxnorm1 * ((expGaus(x, tkxPrime) - expGaus(x, tXe)) / (lX - kxPrime) - (expGaus(x, t1) - expGaus(x, tXe)) / (lX - l1));
      double x3 = fxnorm3 * ((expGaus(x, tkxPrime) - expGaus(x, tXe)) / (lX - kxPrime) - (expGaus(x, t3) - expGaus(x, tXe)) / (lX - l3));
      fx = x1 + x3;


      fs = fs * SiPMQ128;
      ft = ft * SiPMQ128;
      fm = fm * SiPMQE150;

      // glass covered sees only  175
      if (ic == 5)
      {
        fx = fx * SiPMQE175;
        fs = 0;
        ft = 0;
        fm = 0;
        frec = frec * SiPMQE175;
      }
      else if (ic == 12)
      {
        fx = fx * PMTQE175;
        fs = 0;
        ft = 0;
        fm = fm * PMTQE150;
        frec = frec * PMTQE400;
      }
      else
      { // all other sipms 
        fx = fx * SiPMQE175;
        frec = frec * SiPMQE175;
      }
      double mval = fs + ft + fx + frec + fm + bkg + foffset;
      // double f = fs + frec + ft + fx + fm + bkg;


      /*******/
      if (mval < 0)
      {
        //printf("yyy negative model value ibin %i x %E set to 1 \n", j, x);
        mval = 1.;
      }
      double y = buff[ic][j]; // observed
      // we do an NLL
      double yterm = 0.0; // in this case Prob=1 so log=0
      if (y > 0)
        yterm = y - y * log(y);
      f += mval - y * log(mval) - yterm;
      ntScan->Fill(par[1],fx,f);
      // leave warnning printout
      if (isnan(f))
      {
        // printf("nnnn  ibin f is NAN %i x = %E y = %E \n", j, x, y);
        // show();
      }
  }
  /*
  printf(".... f =  %.3E ;;", f);
  for (int ii = 0; ii < NPARS; ++ii)
  {
    printf("  %i %.2E ;", ii, lpar[ii]);
  }
  printf("\n");
  */
  }
}
