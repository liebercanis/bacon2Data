#include "modelFit.hh"
void modelFitTest(){
    int itype = 4;
    modelFit *mf = new modelFit(itype,0);
    mf->show();

    // model fit has pointer to underlying TF1
    TCanvas *can = new TCanvas(Form("model%i",itype) ,Form("model%i",itype));
    mf->fp->Draw();

}
