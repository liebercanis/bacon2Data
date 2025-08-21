#include "TReadGains.hxx"
ClassImp(TReadGains)

    TReadGains::TReadGains()
{
    /**************** define nominal gains ***************/
    nominalGain = 134.786401;     // 170.;     // was 160.0; set Jue 13 2025
    nominalTrigGain = 735.688747; //
    nominalQsumGain = 4940.503519;
    nominalQsumTrigGain = 32056.789775;
    nominalPmtGain = 502.;
    nominalQsumPmtGain = 1713;
    clear();

    // new gain file
    gainFilePeakName = TString(getenv("BOBJ")) + TString("/gainPeakCurrent.root");
    gainFileSumName = TString(getenv("BOBJ")) + TString("/gainSumCurrent.root");
    cout << "read gains from file " << gainFilePeakName << "" << gainFileSumName << endl;
    readPeakGains(gainFilePeakName);
    readSumGains(gainFileSumName);
    printGains();
}

void TReadGains::clear()
{
    sipmPeakGain.clear();
    sipmPeakGainError.clear();
    sipmSumGain.clear();
    sipmSumGainError.clear();
}

bool TReadGains::readPeakGains(TString fileName)
{

    /* define nominal */
    sipmPeakGain.clear();
    sipmPeakGainError.clear();
    sipmPeakGain.resize(NUMCHANNELS);
    sipmPeakGainError.resize(NUMCHANNELS);
    for (int i = 0; i < NUMCHANNELS; ++i)
        sipmPeakGain[i] = nominalGain;
    sipmPeakGain[9] = nominalTrigGain;
    sipmPeakGain[10] = nominalTrigGain;
    sipmPeakGain[11] = nominalTrigGain;
    sipmPeakGain[12] = nominalPmtGain;

    for (int i = 0; i < NUMCHANNELS; ++i)
        sipmPeakGainError[i] = sqrt(sipmPeakGain[i]);

    /* look for gain file */
    bool exists = false;
    FILE *aFile;
    aFile = fopen(fileName.Data(), "r");
    if (aFile)
    {
        fclose(aFile);
        exists = true;
    }

    if (!exists)
    {
        printf(" fopen couldnt open template file %s\n", fileName.Data());
        return false;
    }

    TFile *fin = new TFile(fileName, "readonly");
    if (fin->IsZombie())
    {
        std::cout << "Error opening file" << fileName << std::endl;
        return false;
    }
    cout << " opened sipm gain file " << fileName << endl;
    TGraphErrors *gGain = NULL;
    fin->GetObject("gainPeak", gGain);
    if (gGain == NULL)
    {
        cout << "no gGain in file " << endl;
        return false;
    }
    cout << "found graph named " << gGain->GetName() << " in file " << fileName << endl;

    for (int i = 0; i < gGain->GetN(); ++i)
    {
        int index = int(gGain->GetPointX(i));
        sipmPeakGain[index] = gGain->GetPointY(i);
        sipmPeakGainError[index] = gGain->GetErrorY(i);
    }

    return true;
}

bool TReadGains::readSumGains(TString fileName)
{
    sipmSumGain.clear();
    sipmSumGainError.clear();
    sipmSumGain.resize(NUMCHANNELS);
    sipmSumGainError.resize(NUMCHANNELS);
    for (int i = 0; i < NUMCHANNELS; ++i)
        sipmSumGain[i] = nominalQsumGain;
    sipmSumGain[9] = nominalQsumTrigGain;
    sipmSumGain[10] = nominalQsumTrigGain;
    sipmSumGain[11] = nominalQsumTrigGain;
    sipmSumGain[12] = nominalQsumPmtGain;

    for (int i = 0; i < NUMCHANNELS; ++i)
        sipmSumGainError[i] = sqrt(sipmSumGain[i]);

    /* look for gain file */
    bool exists = false;
    FILE *aFile;
    aFile = fopen(fileName.Data(), "r");
    if (aFile)
    {
        fclose(aFile);
        exists = true;
    }

    if (!exists)
    {
        printf(" fopen couldnt open template file %s\n", fileName.Data());
        return false;
    }

    TFile *fin = new TFile(fileName, "readonly");
    if (fin->IsZombie())
    {
        std::cout << "Error opening file" << fileName << std::endl;
        return false;
    }
    cout << " opened sipm gain file " << fileName << endl;
    TGraphErrors *gGain = NULL;
    fin->GetObject("gainSum", gGain);
    if (gGain == NULL)
    {
        cout << "no gGain in file " << endl;
        return false;
    }
    cout << "found graph named " << gGain->GetName() << " in file " << fileName << endl;

    for (int i = 0; i < gGain->GetN(); ++i)
    {
        int index = int(gGain->GetPointX(i));
        sipmSumGain[index] = gGain->GetPointY(i);
        sipmSumGainError[index] = gGain->GetErrorY(i);
    }

    return true;
}

void TReadGains::printGains()
{
    printf("TReadGains:: %lu gains \n", sipmPeakGain.size());
    for (unsigned long j = 0; j < sipmPeakGain.size(); ++j)
    {
        printf(" chan %lu  peak gain %.4f error %.4f  sum gain %.4f error %.4f \n", j, sipmPeakGain[j], sipmPeakGainError[j],
               sipmSumGain[j], sipmSumGainError[j]);
    }
}
