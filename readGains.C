void readGains(TString fileName = "gains-2023-12-21-13-16.root")
{
    std::vector<double> sipmGain;
    std::vector<double> sipmGainError;
    TFile *fin = new TFile(fileName, "readonly");
    cout << " opened " << fileName << endl;
    //fin->ls();
    TGraphErrors *gGain = NULL;
    fin->GetObject("gGain",gGain);
    if(gGain==NULL) {
        cout << "no gGain in file " <<  endl;
        return;
    }
    cout << "found graph gGain" << endl;
    sipmGain.clear();
    sipmGainError.clear();
    int np = gGain->GetN();
    for (int i = 0; i < np; ++i) {
        double x, y;
        gGain->GetPoint(i, x, y);
        double ey = gGain->GetErrorY(i);
        //printf(" %i x %.0f y %.3f ey %.3f ey/y %.4f \n",i,x,y,ey,ey/y);
        sipmGain.push_back(y);
        sipmGainError.push_back(ey);
    }
    
    printf("stored gains %lu \n", sipmGain.size());
    for (unsigned long j = 0; j < sipmGain.size(); ++j)
        printf(" %lu  gain %.4f error %.4f   \n",j,  sipmGain[j],sipmGainError[j]);
     
}