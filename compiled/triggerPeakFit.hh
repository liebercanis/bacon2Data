#ifndef TRIGGERPEAKFIT_DEFINED
#define TRIGGERPEAKFIT_DEFINED
/*
    June 4 2025
    fit to trigger interaction point and peak value
*/

#include "Math/Vector3D.h"
#include "TMath.h"

static double qsumValue[3]; // qsum input

/* everthing must be in this one routing */
static double peakFit(double *par)
{
    double nLL = 0; // negative log lilelihood value
    // geometry
    double trigRadius = 1.;
    double trigTheta = 55.06 / 360. * TMath::TwoPi(); // 11,10,9
    double trigPhi[3];
    trigPhi[0] = 0.;                           // 9
    trigPhi[1] = 240. / 360. * TMath::TwoPi(); // 10
    trigPhi[2] = 120. / 360. * TMath::TwoPi(); // 11
    double effGeo = 0.028648;
    double SiPMQE128Ham = 0.15;

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

    bool show = true;

    if (show)
    {
        // shift phi for printing.
        double localPhi;
        for (int itrg = 0; itrg < 3; ++itrg)
        {
            localPhi = rsipm[itrg].Phi() * 360. / TMath::TwoPi();
            if (localPhi < 0)
                localPhi += 360.;
            printf("\t trigger sipm %i R,Theta,Phi (%f,%f,%f) \n", itrg, rsipm[itrg].R(), rsipm[itrg].Theta() * 360. / TMath::TwoPi(), localPhi);

            printf("peakFit input %f %f %f \n", qsumValue[0], qsumValue[1], qsumValue[2]);
        }
    }
    return nLL;
}
#endif