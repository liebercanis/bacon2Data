#include "TReadGains.hxx"
ClassImp(TReadGains)

    TReadGains::TReadGains()
{
    /**************** define nominal gains ***************/
    nominalGain = 134.786401;     // 170.;     // was 160.0; set Jue 13 2025
    nominalTrigGain = 735.688747; //
    nominalQsumGain = 4940.503519;
    nominalQsumTrigGain = 32056.789775;
    nominalPmtGain = 165.;      // changed for run 5 was 502.;
    nominalQsumPmtGain = 1900.; // changed for run 5 was 1713;
    clear();
    gPeakGraph = NULL;
    gSumGraph = NULL;

    // new gain file
    gainFileName = TString(getenv("BOBJ")) + TString("/gainsCurrent.root");
    bool gotFile = openFile();
    if (gotFile)
        cout << "TReadGains:: read gains from file " << gainFileName << "" << gainFileName << endl;
    else
        cout << "TReadGains:: error failed to open file  " << gainFileName << "" << gainFileName << endl;

    readPeakGains();
    readSumGains();
    printGains();
}

bool TReadGains::openFile()
{
    /* look for gain file */
    bool exists = false;
    FILE *aFile;
    aFile = fopen(gainFileName.Data(), "r");
    if (aFile)
    {
        fclose(aFile);
        exists = true;
    }

    if (!exists)
    {
        printf(" fopen couldnt open template file %s\n", gainFileName.Data());
        return false;
    }

    TFile *fin = new TFile(gainFileName, "readonly");
    if (fin->IsZombie())
    {
        std::cout << "Error opening file" << gainFileName << std::endl;
        return false;
    }
    cout << " opened sipm gain file " << gainFileName << endl;

    gPeakGraph = NULL;
    fin->GetObject("gainPeak", gPeakGraph);
    if (gPeakGraph == NULL)
    {
        cout << "no gainPeak in file " << endl;
        return false;
    }
    cout << "found graph named " << gPeakGraph->GetName() << " in file " << gainFileName << endl;

    gSumGraph = NULL;
    fin->GetObject("gainSum", gSumGraph);
    if (gSumGraph == NULL)
    {
        cout << "no gainSum in file " << endl;
        return false;
    }
    cout << "found graph named " << gSumGraph->GetName() << " in file " << gainFileName << endl;
    return true;
}

void TReadGains::clear()
{
    sipmPeakGain.clear();
    sipmPeakGainError.clear();
    sipmSumGain.clear();
    sipmSumGainError.clear();
}

bool TReadGains::readPeakGains()
{

    TGraphErrors *gGain = gPeakGraph;
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
    if (!gGain)
        return false;

    printf("peak gains with graph %s \n", gGain->GetName());
    for (int i = 0; i < NUMCHANNELS; ++i)
        sipmPeakGainError[i] = sqrt(sipmPeakGain[i]);

    for (int i = 0; i < gGain->GetN(); ++i)
    {
        int index = int(gGain->GetPointX(i));
        sipmPeakGain[index] = gGain->GetPointY(i);
        sipmPeakGainError[index] = gGain->GetErrorY(i);
    }

    return true;
}

bool TReadGains::readSumGains()
{
    TGraphErrors *gGain = gSumGraph;
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
    if (!gGain)
        return false;

    printf("sum gains with graph %s \n", gGain->GetName());

    for (int i = 0; i < NUMCHANNELS; ++i)
        sipmSumGainError[i] = sqrt(sipmSumGain[i]);

    // protect against bad gain values
    for (int i = 0; i < gGain->GetN(); ++i)
    {
        int index = int(gGain->GetPointX(i));
        // protect against bad values
        double val = gGain->GetPointY(i);
        if (i < 9)
        {
            if (val > nominalQsumGain / 2.)
            {
                sipmSumGain[index] = gGain->GetPointY(i);
                sipmSumGainError[index] = gGain->GetErrorY(i);
            }
            else
                printf("TReadGain WARNING gain val %f too small compared to %f\n", val, sipmSumGain[index]);
        }
        else if (i > 8 && i < 12)
        {
            if (val > nominalQsumTrigGain / 2.)
            {
                sipmSumGain[index] = gGain->GetPointY(i);
                sipmSumGainError[index] = gGain->GetErrorY(i);
            }
            else
                printf("TReadGain WARNING gain val %f too small compared to %f\n", val, sipmSumGain[index]);
        }
        else if (i == 12)
        {
            if (val > nominalQsumPmtGain / 2.)
            {
                sipmSumGain[index] = gGain->GetPointY(i);
                sipmSumGainError[index] = gGain->GetErrorY(i);
            }
            else
                printf("TReadGain WARNING gain val %f too small compared to %f\n", val, sipmSumGain[index]);
        }
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
