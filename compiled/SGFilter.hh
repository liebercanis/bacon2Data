/* implimentation of one-dimensional Savitzky-Golay filter in C++
from Vishal implimented by MG*/
#ifndef SAVITZKYGOLAY_H
#define SAVITZKYGOLAY_H
// #include <vector>
// using std::vector;
#include "TH1D.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstddef>

class SavitzkyGolay
{

public:
    SavitzkyGolay() {
    };
    ~SavitzkyGolay() {
    };

    // Calculates the Gram Polynomial (s=0), or its s'th derivative evaluated at i, order k, over 2m+1 points
    double GramPoly(int i, int m, int k, int s)
    {
        // Calculates the Gram Polynomial (s=0), or its s'th derivative evaluated at i, order k, over 2m+1 points
        double gp_val;
        if (k > 0)
        {
            gp_val = (4.0 * k - 2.0) / (k * (2.0 * m - k + 1.0)) * (i * GramPoly(i, m, k - 1, s) + s * GramPoly(i, m, k - 1, s - 1)) - ((k - 1.0) * (2.0 * m + k)) / (k * (2.0 * m - k + 1.0)) * GramPoly(i, m, k - 2, s);
        }
        else
        {
            if (k == 0 && s == 0)
            {
                gp_val = 1.0;
            }
            else
            {
                gp_val = 0.0;
            } // end of if k = 0 && s = 0
        } // end of if k > 0
        return gp_val;
    } // end of GramPoly function

    // Calculates the generalised factorial (a)(a-l)...(a-b+l)
    int GenFact(int a, int b)
    {
        // Calculates the generalised factorial (a)(a-l)...(a-b+l)
        int gf = 1.0;
        for (int jj = (a - b + 1); jj < a + 1; jj++)
        {
            gf = gf * jj;
        }
        return (gf);
    } // end of GenFact function

    // Calculates the weight of the i'th data point for the t'th Least-Square point of the s'th derviative, over 2m+l points, order
    double Weight(int i, int t, int m, int n, int s)
    {
        // Calculates the weight of the i'th data point for the t'th Least-Square point of the s'th derviative, over 2m+l points, order
        double sum = 0.0;
        for (int k = 0; k < n + 1; k++)
        {
            sum = sum + (2.0 * k + 1.0) * ((double)GenFact(2 * m, k) / GenFact(2 * m + k + 1, k + 1)) * GramPoly(i, m, k, 0) * GramPoly(t, m, k, s);
        }
        return sum;
    } // end of Weight function

    // The Savitzky-Golay filter
    // Default destructor
    // Savitzky Golay filter parameter mwindow (filterLength = 2*mwindow+1)
    // Savitzky Golay filter parameter npoly ( order of polynomial )
    std::vector<double> SavGolFilter(TH1D *hist, int mwindow = 3, int npoly = 3)
    {
        std::vector<std::vector<double>> weight;
        weight.resize(3);
        // double W[3][filterLength];
        // collect the weights
        for (int jj = 0; jj < weight.size(); jj++)
        {
            weight[jj].clear();
            for (int kk = -mwindow; kk < mwindow + 1; kk++)
            {
                weight[jj].push_back(Weight(kk, 0, mwindow, npoly, jj));
                // std::cout << " weight bin " << jj << " w= " << weight[jj][weight[jj].size() - 1] << std::endl;
            }
        }

        // use weights to smooth
        // makea vector for doing sum
        std::vector<double> yfilt;
        yfilt.resize(hist->GetNbinsX());
        std::fill(yfilt.begin(), yfilt.end(), 0);

        int filterLength = 2 * mwindow + 1;
        // loop over input histogram bins
        for (int i = 0; i < hist->GetNbinsX(); i++)
        {
            //  loop over filter window starting at bin i avoid end of vector yfilt
            for (int kk = 0; kk < filterLength; kk++)
            {
                int ibin = i + kk - mwindow;
                if (ibin < 0)
                    continue;
                if (ibin > hist->GetNbinsX())
                    continue;
                yfilt[i] += hist->GetBinContent(ibin) * weight[0][kk];
                // if (kk == filterLength - 1 && ibin / 100 * 100 == ibin)
                //     printf("SGSUM bin=%i kk=%i weight %f bin %f weighted sum %f \n", ibin, kk, weight[0][kk], hist->GetBinContent(i), yfilt[i]);
            }
        }
        return yfilt;
    }
};
#endif