int nsamples = 1024;
TFile *fout;
TH1D *htemplate;
TH1D *htemplateFFT;
std::vector<std::complex<double>> templateTransform;
std::vector<std::complex<double>> forwardFFT(std::vector<double> rdigi);
std::vector<Double_t> backwardFFT(std::vector<std::complex<double>> VectorComplex);
std::vector<Double_t> SPEdigi;
TRandom3 *ran3;
double QPE = 3000.;

/// The fft class to take the fourier transform.
TVirtualFFT *fFFT;
/// The fft class to take the inverse fourier transform.
TVirtualFFT *fInverseFFT;

TH1D *getTemplate(int ichan)
{
    TH1D *hist = NULL;
    TString fileName = TString("../bobj/templates-2023-05-01-15-06.root");
    TFile *f1 = new TFile(fileName, "readonly");
    if (f1->IsZombie())
    {
        printf(" no  file for %s \n", fileName.Data());
        return hist;
    }
    f1->GetObject(Form("QPEShapeChan%i", ichan), hist);
    if (hist)
        printf(" got template %s  from file %s \n", hist->GetName(), fileName.Data());
    return hist;
}

void fftTest(int npoints = 10)
{

    ran3 = new TRandom3(65539);

    fout = new TFile("fftTest.root", "recreate");
    htemplate = new TH1D("template", "template", nsamples, 0, nsamples);
    htemplateFFT = new TH1D("templateFFT", "templateFFT", nsamples / 2, 0, nsamples / 2);
    TH1D *hinputWave = new TH1D("inputWave", "inputWave", nsamples, 0, nsamples);
    TH1D *hinputFFT = new TH1D("inputFFT", "inputFFT", nsamples / 2, 0, nsamples / 2);
    TH1D *hconvolutionFFT = new TH1D("convolutionFFT", "convolutionFFT", nsamples / 2, 0, nsamples / 2);
    TH1D *hsignalWave = new TH1D("signalWave", "signalWave", nsamples, 0, nsamples);
    TH1D *hdeconvWave = new TH1D("deconvlWave", "deconvWave", nsamples, 0, nsamples);

    for (int i = 0; i < npoints; ++i)
    {
        double t = 30.*TMath::Log(1. / (1. - ran3->Rndm()));
        printf(" %i %f bin %i \n", i, t, hinputWave->FindBin(t));
        hinputWave->SetBinContent(hinputWave->FindBin(t), QPE);
    }


    // initialize fft
    fFFT = TVirtualFFT::FFT(1, &nsamples, "R2C M K");
    fInverseFFT = TVirtualFFT::FFT(1, &nsamples, "C2R M K");

    TH1D *SPETemplate = getTemplate(6);
    fout->Add(SPETemplate);
    printf(" template integral %f\n ", SPETemplate->Integral());

    SPEdigi.resize(nsamples);

    for (int ibin = 0; ibin < SPETemplate->GetNbinsX(); ++ibin)
    {
        SPEdigi[ibin] = SPETemplate->GetBinContent(ibin);
    }

    for (int ibin = 0; ibin < SPETemplate->GetNbinsX(); ++ibin)
        htemplate->SetBinContent(ibin, SPEdigi[ibin]);

    // make transorm
    std::vector<std::complex<double>> templateTransform = forwardFFT(SPEdigi);

    // fill histos
    for (unsigned isample = 0; isample < SPEdigi.size(); isample++)
        htemplate->SetBinContent(isample, SPEdigi[isample]);

    for (int i = 0; i < nsamples / 2; ++i)
        htemplateFFT->SetBinContent(i, std::abs(templateTransform[i]));

    // FFT input wave
    std::vector<double> digi;
    digi.resize(nsamples);
    for (int ibin = 0; ibin < hinputWave->GetNbinsX(); ++ibin)
        digi[ibin] = hinputWave->GetBinContent(ibin);
    //
    std::vector<std::complex<double>> inputWaveTransform = forwardFFT(digi);
    for (int i = 0; i < nsamples / 2; ++i)
        hinputFFT->SetBinContent(i, std::abs(inputWaveTransform[i]));

    // convolve with shape
    unsigned maxFrequency = inputWaveTransform.size();
    std::vector<std::complex<double>> convolutionTransform;
    convolutionTransform.resize(nsamples);
    printf(" max freq = %u  \n", maxFrequency);
    // apply shape FFT convolution here
    for (unsigned iw = 1; iw < maxFrequency; ++iw)
    {
        convolutionTransform[iw] = inputWaveTransform[iw] * templateTransform[iw];
    }

    for (int i = 0; i < nsamples / 2; ++i)
        hconvolutionFFT->SetBinContent(i, std::abs(convolutionTransform[i]));

    std::vector<double> fdigi = backwardFFT(convolutionTransform);
    for (unsigned i = 0; i < fdigi.size(); ++i)
        fdigi[i] /= double(nsamples);

    for (int ibin = 0; ibin < hsignalWave->GetNbinsX(); ++ibin)
        hsignalWave->SetBinContent(ibin, fdigi[ibin]);

    // deconvolute back
    std::vector<std::complex<double>> deconvTransform;
    deconvTransform.resize(nsamples);
    for (unsigned iw = 1; iw < maxFrequency; ++iw)
    {
        deconvTransform[iw] = convolutionTransform[iw] / templateTransform[iw];
    }
    std::vector<double> gdigi = backwardFFT(deconvTransform);
    for (unsigned i = 0; i < gdigi.size(); ++i)
        gdigi[i] /= double(nsamples);
    for (int ibin = 0; ibin < hdeconvWave->GetNbinsX(); ++ibin)
        hdeconvWave->SetBinContent(ibin, gdigi[ibin]);
}

std::vector<std::complex<double>> forwardFFT(std::vector<double> rdigi)
{
    std::vector<std::complex<double>> VectorComplex;
    for (int is = 0; is < nsamples; ++is)
        fFFT->SetPoint(is, rdigi[is]);
    fFFT->Transform();

    std::vector<Double_t> realVec, imVec;
    for (int i = 0; i < nsamples; ++i)
    {
        double rl, im;
        fFFT->GetPointComplex(i, rl, im);
        std::complex<double> c(rl, im); //.real or .imag accessors
        VectorComplex.push_back(c);
    }
    return VectorComplex;
}

std::vector<Double_t> backwardFFT(std::vector<std::complex<double>> VectorComplex)
{
    std::vector<Double_t> Signal;
    for (int is = 0; is < nsamples; ++is)
    {
        fInverseFFT->SetPoint(is, VectorComplex[is].real(), VectorComplex[is].imag());
    }
    fInverseFFT->Transform();

    for (int i = 0; i < nsamples; ++i)
    {
        double rl = fInverseFFT->GetPointReal(i);
        Signal.push_back(rl);
    }
    return Signal;
}
