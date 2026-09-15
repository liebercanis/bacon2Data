#include "TReadGains.hxx"
ClassImp(TReadGains)

    TReadGains::TReadGains(bool useFile)
{
    readFromFile = useFile;
    /**************** define nominal gains ***************/
    nominalGain = 134.786401;     // 170.;     // was 160.0; set Jue 13 2025
    nominalTrigGain = 735.688747; //
    nominalQsumGain = 4940.503519;
    nominalQsumTrigGain = 32056.789775;
    nominalPmtGain = 165.;       // changed for run 5 was 502.;
    nominalQsumPmtGain = 760.00; // changed for 04_16_2026 changed for run 5 was 1713;
    clear();
    gPeakGraph = NULL;
    gSumGraph = NULL;

    // new gain file
    gainFileName = TString(getenv("BOBJ")) + TString("/gainsCurrent.root");
    bool gotFile = openFile();
    if (gotFile && readFromFile)
        cout << "TReadGains:: read gains from file " << gainFileName << "" << gainFileName << endl;
    else if (!gotFile)
        cout << "TReadGains:: error failed to open file  " << gainFileName << "" << gainFileName << endl;
    else
        cout << "TReadGains:: using nominal gains " << endl;

    readPeakGains();
    readSumGains();
    relativeLEDEff.resize(NUMCHANNELS);
    // from LED run as eff_i - mean
    // devide to correct
    relativeLEDEff[0] = 1. + 0.245877;
    relativeLEDEff[1] = 1. + 0.0567434;
    relativeLEDEff[2] = 1. + 0.21173;
    relativeLEDEff[3] = 1. - 0.217802;
    relativeLEDEff[4] = 1. + 0.274232;
    relativeLEDEff[5] = 1. + 0.00278381;
    relativeLEDEff[6] = 1. - 0.359765;
    relativeLEDEff[7] = 1. - 0.184547;
    relativeLEDEff[8] = 1. - 0.0292534;
    relativeLEDEff[9] = 1.0;
    relativeLEDEff[10] = 1.0;
    relativeLEDEff[11] = 1.0;
    relativeLEDEff[12] = 1.0;

    // relativeEff from zero PPM run
    relativeEff.resize(NUMCHANNELS);
    for (unsigned i = 0; i < relativeEff.size(); ++i)
        relativeEff[i] = 1;

    // with 0.15 PPM
    relativeEff[1] = 1.0 + -0.320668;
    relativeEff[2] = 1.0 + 0.127058;
    relativeEff[3] = 1.0 + 0.034291;
    relativeEff[4] = 1.0 + 0.081117;
    relativeEff[5] = 1.0 + -0.022577;
    relativeEff[6] = 1.0 + 0.051182;
    relativeEff[7] = 1.0 + 0.049597;
    relativeEff[9] = 1.0 + 2.050647;
    relativeEff[10] = 1.0 + 2.223320;
    relativeEff[11] = 1.0 + 1.981977;
    relativeEff[12] = 1.0 + -0.853489;

    // with 0.05 PPM
    /*
    relativeEff[1] = 1. + -0.312058;
    relativeEff[2] = 1. + 0.141342;
    relativeEff[3] = 1. + 0.002901;
    relativeEff[4] = 1. + 0.048306;
    relativeEff[5] = 1. + -0.052241;
    relativeEff[6] = 1. + 0.086694;
    relativeEff[7] = 1. + 0.085056;
    relativeEff[9] = 1. + 3.337998;
    relativeEff[10] = 1. + 3.583538;
    relativeEff[11] = 1. + 3.240349;
    relativeEff[12] = 1. + -0.839874;
    */

    /*
        0 PPM
        //relativeEff[0] = 1. + -0.832702;
        relativeEff[1] = 1. + -0.421601;
        relativeEff[2] = 1. + -0.040397;
        relativeEff[3] = 1. + -0.022906;
        relativeEff[4] = 1. + 0.021330;
        relativeEff[5] = 1. + -0.076629;
        relativeEff[6] = 1. + 0.271060;
        relativeEff[7] = 1. + 0.269144;
        //relativeEff[8] = 1. + -0.922346;
        relativeEff[9] = 1. + 5.393614;
        relativeEff[10] = 1. + 5.755505;
        relativeEff[11] = 1. + 5.249692;
        relativeEff[12] = 1. + -0.881690;
        */

    relativeNorm.resize(NUMCHANNELS);
    for (unsigned i = 0; i < relativeNorm.size(); ++i)
        relativeNorm[i] = 1;

    /** from limit range fit  old
    relativeNorm[0] = 1.0 + 2.203244;
    relativeNorm[2] = 1.0 + -0.045537;
    relativeNorm[3] = 1.0 + -0.106738;
    relativeNorm[4] = 1.0 + -0.137989;
    relativeNorm[5] = 1.0 + -0.204397;
    relativeNorm[6] = 1.0 + 0.151084;
    relativeNorm[7] = 1.0 + 0.152126;
    */
    // after baseline mode subtraction sept 15 2026 skip channels 0, 8

    relativeEff[1] = 1. + -0.446354;
    relativeEff[2] = 1. + -0.076883;
    relativeEff[3] = 1. + -0.030821;
    relativeEff[4] = 1. + 0.017561;
    relativeEff[5] = 1. + -0.090557;
    relativeEff[6] = 1. + 0.316456;
    relativeEff[7] = 1. + 0.310597;

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
    if (!gGain || !readFromFile)
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
    if (!gGain || !readFromFile)
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
            if (val > nominalQsumGain / 10.)
            {
                sipmSumGain[index] = gGain->GetPointY(i);
                sipmSumGainError[index] = gGain->GetErrorY(i);
            }
            else
                printf("TReadGain WARNING gain channel %i val %f too small compared to %f\n", i, val, nominalQsumGain);
        }
        else if (i > 8 && i < 12)
        {
            if (val > nominalQsumTrigGain / 10.)
            {
                sipmSumGain[index] = gGain->GetPointY(i);
                sipmSumGainError[index] = gGain->GetErrorY(i);
            }
            else
                printf("TReadGain WARNING gain channel %i val %f too small compared to %f\n", i, val, nominalQsumTrigGain);
        }
        else if (i == 12)
        {
            if (val > nominalQsumPmtGain / 2.)
            {
                sipmSumGain[index] = gGain->GetPointY(i);
                sipmSumGainError[index] = gGain->GetErrorY(i);
            }
            else
                printf("TReadGain WARNING gain channel %i val %f too small compared to %f\n", i, val, sipmSumGain[index]);
        }
    }

    return true;
}

double TReadGains::getNominalPeak(unsigned j)
{
    double nominal = 0;
    if (j >= 0 && j < 9)
    {
        nominal = nominalGain;
    }
    if (j > 8 && j < 12)
    {
        nominal = nominalTrigGain;
    }

    if (j == 12)
    {
        nominal = nominalPmtGain;
    }
    return nominal;
}

double TReadGains::getNominalSum(unsigned j)
{
    double nominal = 0;
    if (j >= 0 && j < 9)
    {
        nominal = nominalQsumGain;
    }
    if (j > 8 && j < 12)
    {
        nominal = nominalQsumTrigGain;
    }

    if (j == 12)
    {
        nominal = nominalQsumPmtGain;
    }
    return nominal;
}

void TReadGains::printGains()
{
    printf("TReadGains:: %lu gains \n", sipmPeakGain.size());
    double peakToNominal = 1;
    double sumToNominal = 1;
    for (unsigned long j = 0; j < sipmPeakGain.size(); ++j)
    {
        if (j >= 0 && j < 9)
        {
            peakToNominal = sipmPeakGain[j] / nominalGain;
            sumToNominal = sipmSumGain[j] / nominalQsumGain;
        }
        if (j > 8 && j < 12)
        {
            peakToNominal = sipmPeakGain[j] / nominalTrigGain;
            sumToNominal = sipmSumGain[j] / nominalQsumTrigGain;
        }

        if (j == 12)
        {
            peakToNominal = sipmPeakGain[j] / nominalPmtGain;
            sumToNominal = sipmSumGain[j] / nominalQsumPmtGain;
        }

        printf(" chan %lu  peak gain %.4f error %.4f  sum gain %.4f error %.4f peak to nominal %.4f sum to nominal %.4f  relative %.4f  \n", j, sipmPeakGain[j], sipmPeakGainError[j],
               sipmSumGain[j], sipmSumGainError[j], peakToNominal, sumToNominal, relativeEff[j]);
    }
}