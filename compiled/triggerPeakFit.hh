#ifndef TRIGGERPEAKFIT_DEFINED
#define TRIGGERPEAKFIT_DEFINED
/*
    June 4 2025
    fit to trigger interaction point and peak value
*/

#include "Math/Vector3D.h"
#include "TMath.h"

enum
{
    nTriggerSipms = 3
};
static double peakFitQsum[nTriggerSipms];  // qsum input
static double peakMeanQsum[nTriggerSipms]; // qsum input
static bool triggerPeakFitShow = false;

/* everthing must be in this one routing */
static double peakFit(double *par)
{
    // geometry
    double trigRadius = 1.;
    double trigTheta = 55.06 / 360. * TMath::TwoPi(); // 11,10,9
    double trigPhi[nTriggerSipms];
    trigPhi[0] = 0.;                           // 9
    trigPhi[1] = 240. / 360. * TMath::TwoPi(); // 10
    trigPhi[2] = 120. / 360. * TMath::TwoPi(); // 11
    double effGeo = 1.280427E-01;
    double SiPMQE128Ham = 0.15;
    double area = pow(0.6, 2.);
    double fillFac = 0.7;
    double eff = effGeo * SiPMQE128Ham * fillFac; // total efficiency

    // sipm position vectors
    ROOT::Math::XYZVector rsipm[3];
    // construct XYZ coordinates
    double cost = cos(trigTheta);
    double sint = sin(trigTheta);
    double cosp = cos(trigPhi[0]);
    double sinp = sin(trigPhi[0]);
    double x = trigRadius * sint * cosp;
    double y = trigRadius * sint * sinp;
    double z = trigRadius * cost;
    rsipm[0] = ROOT::Math::XYZVector(x, y, z);
    //
    cosp = cos(trigPhi[1]);
    sinp = sin(trigPhi[1]);
    x = trigRadius * sint * cosp;
    y = trigRadius * sint * sinp;
    z = trigRadius * cost;
    rsipm[1] = ROOT::Math::XYZVector(x, y, z);
    //
    cosp = cos(trigPhi[2]);
    sinp = sin(trigPhi[2]);
    x = trigRadius * sint * cosp;
    y = trigRadius * sint * sinp;
    z = trigRadius * cost;
    rsipm[2] = ROOT::Math::XYZVector(x, y, z);

    double totalPhotonYield = par[0]; // total expected photons
    // printf("xxxxxxx par 0 %f  %f\n", par[0], totalPhotonYield);
    //  interaction point
    ROOT::Math::XYZVector ipVector = ROOT::Math::XYZVector(par[1], par[2], par[3]);

    // compute mean values
    for (int i = 0; i < nTriggerSipms; ++i)
    {
        // relatve postion vector
        ROOT::Math::XYZVector relative = rsipm[i] - ipVector;
        // correct solid angle
        ROOT::Math::XYZVector runit = rsipm[i].Unit();
        ROOT::Math::XYZVector ounit = relative.Unit();
        double cos = runit.Dot(ounit);
        if (cos < 0)
            cos = 0;
        peakMeanQsum[i] = totalPhotonYield * cos * SiPMQE128Ham * fillFac * area / (4. * TMath::Pi() * relative.Mag2());
        if (triggerPeakFitShow)
            printf("triggerPeakFit:: trig %i cos %f total eff %E peakMeanQsum %f \n", i, cos, eff, peakMeanQsum[i]);
    }

    /* compute NLL see TH1F reference
        NLL = sum ( f_i+y_i log(y_i/f_i) -y_i )
        f_i = observed
        y_i = epected
    */
    double nLL = 0; // negative log lilelihood value
    for (int i = 0; i < nTriggerSipms; ++i)
    {
        if (peakFitQsum[i] > 0)
            nLL += peakFitQsum[i] + peakMeanQsum[i] * log(peakMeanQsum[i] / peakFitQsum[i]) - peakMeanQsum[i];
    }

    if (triggerPeakFitShow)
    {
        // shift phi for printing.
        double localPhi;
        printf("triggerPeakFit:: paramters total photons %f IP r,theta,phi =  (%f,%f,%f) \n", par[0], par[1], par[2], par[3]);
        for (int itrg = 0; itrg < 3; ++itrg)
        {
            localPhi = rsipm[itrg].Phi() * 360. / TMath::TwoPi();
            if (localPhi < 0)
                localPhi += 360.;
            printf("\t trigger sipm %i R,Theta,Phi (%f,%f,%f) expected photons %f \n", itrg, rsipm[itrg].R(), rsipm[itrg].Theta() * 360. / TMath::TwoPi(), localPhi, peakMeanQsum[itrg]);

            printf("\t peakFit peakMeanQsum %f %f %f input qsum %f %f %f  \n", peakMeanQsum[0], peakMeanQsum[1], peakMeanQsum[2], peakFitQsum[0], peakFitQsum[1], peakFitQsum[2]);
        }
    }
    return nLL;
}
#endif