#include "modelFit.hh"

void modelFitTest(int ichan=0 )
{
    TString names[NTYPES];
    names[0] = TString("singlet");
    names[1] = TString("triplet");
    names[2] = TString("mixed");
    names[3] = TString("xenon");
    names[4] = TString("total");
    vector<modelFit *> mfits;

    double ppm = 0.05;
    for (int itype = 0; itype < NTYPES; ++itype)
        mfits.push_back(new modelFit(itype, ichan, ppm));

    mfits[NTYPES - 1]->show();
    mfits[0]->fp->SetLineColor(kRed);
    mfits[1]->fp->SetLineColor(kBlue);
    mfits[2]->fp->SetLineColor(kGreen);
    mfits[3]->fp->SetLineColor(kOrange);
    mfits[4]->fp->SetLineColor(kBlack);

    vector<double> sum;
    sum.resize(mfits.size());
    printf(" SUMS for channel %i \n",ichan);
    for (int j = 0; j < sum.size(); ++j)
    {
        sum[j] = mfits[j]->fp->Integral(1400, 6000);
        printf(" type %i %s sum %.3E \n", j, names[j].Data(), sum[j]);
    }


    // model fit has pointer to underlying TF1
    TCanvas *can = new TCanvas(Form("modelChan%iPPM%f", ichan, ppm), Form("modelChan%iPPM%f", ichan,ppm));
    can->SetLogy();

    for (int itype = 0; itype < NTYPES; ++itype)
    {
        mfits[itype]->fp->GetYaxis()->SetRangeUser(1E-7,0.5);
        if (itype == 0)
            mfits[itype]->fp->Draw();
        else
            mfits[itype]->fp->Draw("same");
    }
    can->BuildLegend();
    can->Print(".png");
}
