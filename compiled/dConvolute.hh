/* implimentation of one-dimensional Savitzky-Golay filter in C++
from Vishal implimented by MG*/
#ifndef DCONVOLUTE_H
#define DCONVOLUTE_H
// #include <vector>
// using std::vector;
#include <iostream>
#include <fstream>
#include <string>
#include <cstddef>
#include <TVirtualFFT.h>

class dConvolute
{

public:
    TVirtualFFT *fFFT;
    TVirtualFFT *fInverseFFT;
    // consruct with fft initialization
    dConvolute(TVirtualFFT *theFFT, TVirtualFFT *theInverseFFT)
    {
        fFFT = theFFT;
        fInverseFFT = theInverseFFT;
    };
    ~dConvolute() {
    };

    /* method to do derivativer */
    std::vector<double> differentiate(std::vector<double> digi, double diffStep = 1)
    {
        std::vector<double> ddigi;
        ddigi.resize(digi.size());
        Double_t sump = 0;
        Double_t summ = 0;
        unsigned nsamples = digi.size();
        ddigi[0] = 0; // first entry is zero
        for (unsigned i = 1; i < nsamples; ++i)
        {
            // sum limit
            int maxSum = diffStep;
            if (i < diffStep)
                maxSum = i;
            if (nsamples - 1 - i < diffStep)
                maxSum = nsamples - 1 - i;
            //
            sump = 0;
            for (unsigned j = 0; j < maxSum; ++j)
            {
                sump += digi[i + 1 + j];
            }
            summ = 0;
            for (unsigned j = 0; j < maxSum; ++j)
            {
                summ += digi[i - 1 - j];
            }
            // if(verbose) printf(" hitFinder::differentiate bin %i maxSum %u sump %E summ %E \n",i,maxSum,sump,summ);
            ddigi[i] = sump - summ;
        }
        return ddigi;
    }
    /* fft method */
    // input is real wave to transform output is complex FFT
    // dResponse is fft of derivative of response function
    std::vector<std::complex<double>> FFT(std::vector<double> vin, std::vector<std::complex<double>> dResponse)
    {
        // do transform
        std::vector<std::complex<double>> complexVector;
        for (unsigned is = 0; is < vin.size(); ++is)
            fFFT->SetPoint(is, vin[is]);
        fFFT->Transform(); //

        // build return vector
        std::vector<Double_t> realVec, imVec;
        for (int i = 0; i < int(vin.size()); ++i)
        {
            double rl, im;
            fFFT->GetPointComplex(i, rl, im);
            std::complex<double> c(rl, im);
            // c = c / sqrt(double(nsamples)); // normalize
            complexVector.push_back(c);
        }

        // if needed, convolute with fft of derivative response
        if (dResponse.size() > 0)
            for (int i = 0; i < complexVector.size(); ++i)
                complexVector[i] *= dResponse[i];

        return complexVector;
    }

    /* inverse fft method */
    // input is complex FFT output is complex FFT
    // filter is Wiener, shift is time shift
    std::vector<double> inverseFFT(std::vector<std::complex<double>> complexInput,
                                   std::vector<std::complex<double>> dNoise,
                                   std::vector<std::complex<double>> dResponse,
                                   bool filter = true, int shift = 0)
    {
        std::vector<double> realOutput;
        int nsamples = (int)complexInput.size();
        // deconvolution is C[w]=H[w]/G[w] with Wiener from wikipedia https://en.wikipedia.org/wiki/Wiener_deconvolution
        if (dResponse.size() > 0)
        {
            // loop over frequency bins
            for (int i = 0; i < complexInput.size(); ++i)
            {
                std::complex<double> W;
                if (filter)
                {
                    double noiseRatio = double(std::norm(dNoise[i])) / double(std::norm(complexInput[i]));
                    W = std::conj(dResponse[i]) / (std::norm(dResponse[i]) + noiseRatio);
                }
                else
                {
                    W = 1. / dResponse[i];
                }
                // multiply by weight
                complexInput[i] *= W;
            }
        }
        // if applying time shift if needed
        if (shift != 0)
        {
            for (int i = 0; i < complexInput.size(); ++i)
            {
                std::complex<double> phase = std::exp(1i * double(shift) * double(i) / double(nsamples) * 2. * TMath::Pi());
                complexInput[i] *= phase;
            }
        }

        // do inverse
        for (int is = 0; is < nsamples; ++is)
            fInverseFFT->SetPoint(is, complexInput[is].real(), complexInput[is].imag());

        fInverseFFT->Transform();
        /*
        ** FFTW computes an unnormalized transform, in that there is no coefficient in front of the summation in the DFT.
        ** In other words, applying the forward and then the backward transform will multiply the input by n.
        ** */

        /* build return vector */
        for (int i = 0; i < nsamples; ++i)
        {
            double rl, im;
            fInverseFFT->GetPointComplex(i, rl, im);
            // normalize
            rl /= double(nsamples);
            realOutput.push_back(rl);
        }
        return realOutput;
    }
};
#endif