#include <ctime>
#include <iostream>
#include <iterator>
#include <locale>
#include <TString.h>
#include <TROOT.h>
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TDatime.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TF1.h"

//
#include "anaRun.cc"

std::vector<TString> fileList;
std::vector<double> filenum;
std::vector<double> efilenum;
std::vector<double> fileTime;
std::vector<TDatime> fileDatime;
TBFile *bf;
std::vector<vector<double>> vecQsum;
std::vector<vector<double>> vecEQsum;
std::vector<vector<double>> vecQPE;
std::vector<vector<double>> vecEQPE;
std::vector<vector<double>> vSlope;
std::vector<vector<double>> vESlope;
std::vector<TH1D *> runQSum;
// summed waves
std::vector<TH1D *> hSumWave;
std::vector<TH1D *> hSumHitWave;
std::map<int, TH1D *> histMap;
TDirectory *sumDir;
TDirectory *foutSumDir;
TFile *fin;
TFile *fout;
bool first = true;
TH1D *hQPEChan;
TH1D *eventCount;

void setTimeGraph(TMultiGraph *mg, TString ylabel)
{
  mg->GetXaxis()->SetTimeDisplay(1);
  mg->GetXaxis()->SetNdivisions(-220);
  mg->GetXaxis()->SetTimeFormat("%d:%H");
  mg->GetXaxis()->SetTimeOffset(0, "gmt");
  mg->GetYaxis()->SetTitle(ylabel);
}

Int_t get_month_index(TString name)
{
  int imonth;
  if (name.EqualTo("Jan"))
    imonth = 1;
  if (name.EqualTo("Feb"))
    imonth = 2;
  if (name.EqualTo("Mar"))
    imonth = 3;
  if (name.EqualTo("Apr"))
    imonth = 4;
  if (name.EqualTo("May"))
    imonth = 5;
  if (name.EqualTo("Jun"))
    imonth = 6;
  if (name.EqualTo("Jul"))
    imonth = 7;
  if (name.EqualTo("Aug"))
    imonth = 8;
  if (name.EqualTo("Sep"))
    imonth = 9;
  if (name.EqualTo("Oct"))
    imonth = 10;
  if (name.EqualTo("Nov"))
    imonth = 11;
  if (name.EqualTo("Dec"))
    imonth = 12;
  return imonth;
}

TDatime getTime()
{
  TString strfileTime = bf->modified;
  TString tempstring = TString(strfileTime(20, 4));
  int fileYear = tempstring.Atoi();
  tempstring = TString(strfileTime(4, 3));
  int fileMonth = get_month_index(tempstring);
  tempstring = TString(strfileTime(8, 2));
  int fileDay = tempstring.Atoi();
  tempstring = TString(strfileTime(11, 2));
  int fileHour = tempstring.Atoi();
  tempstring = TString(strfileTime(14, 2));
  int fileMin = tempstring.Atoi();
  tempstring = TString(strfileTime(17, 2));
  int fileSec = tempstring.Atoi();
  TDatime datime;
  datime.Set(fileYear, fileMonth, fileDay, fileHour, fileMin, fileSec);
  printf("FileYear = %u , FileMonth = %u , FileDay = %u , FileHour = %u , FileMin = %u , FileSec = %u \n", fileYear, fileMonth, fileDay, fileHour, fileMin, fileSec);
  return datime;
}
void QPEFits(unsigned ifile)
{
  if (!sumDir)
  {
    printf("addSumHistos no sumDir!!! \n");
    return;
  }

  // sumDir->ls();

  TList *sumList = sumDir->GetListOfKeys();
  TIter next(sumList);
  TKey *key;
  printf("QPEFIT %u \n", sumList->GetEntries());
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    std::string name = string(h->GetName());

    // summed waves
    TH1D *hAdd;
    if (name.find("sumWave") != std::string::npos && name.find("Bad") == std::string::npos)
    {
      hAdd = (TH1D *)h->Clone(Form("%s-file%i", h->GetName(), ifile));
      cout << " add to foutSumDir " << hAdd->GetName() << endl;
      foutSumDir->Add(hAdd);
      hSumWave.push_back(hAdd);
    }
    if(name.find("sumHitWave") != std::string::npos){
      hAdd = (TH1D *)h->Clone(Form("%s-file%i", h->GetName(), ifile));
      hSumHitWave.push_back(hAdd);
      foutSumDir->Add(hAdd);
    }
    fout->cd();
    // channel qsum
    if (name.find("QSumChan") == std::string::npos)
      continue;
    string chan = name.substr(name.find_last_of("n") + 1);
    int ichan = stoi(chan);
    cout << " addSumHistos clone " << name << " chan " << ichan << endl;
    TH1D *hclone = (TH1D *)h->Clone(Form("fitQSumCh%i", ichan));

    // fit QPE
    double xlow = 2000.;
    double xhigh = 5000;
    if (ichan == 9 || ichan == 10 || ichan == 11)
    {
      xlow = 4.E4;
      xhigh = 8.E4;
    }
    hclone->Fit("gaus", " ", " ", xlow, xhigh);
    TF1 *gfit = (TF1 *)hclone->GetListOfFunctions()->FindObject("gaus");
    if (gfit)
    {
      hQPEChan->SetBinContent(ichan, gfit->GetParameter(1));
      hQPEChan->SetBinError(ichan, gfit->GetParError(1));
      printf(" ****** chan %i low %E hight %E QPE %E\n", ichan, xlow, xhigh, gfit->GetParameter(1));
      vecQPE[ichan].push_back(gfit->GetParameter(1));
      vecEQPE[ichan].push_back(gfit->GetParError(1));
    }
    else
    {
      vecQPE[ichan].push_back(0);
      vecEQPE[ichan].push_back(0);
    }
    // fout->Add(hclone);
  }
}

void addSumHistos()
{
  if (!sumDir)
  {
    printf("addSumHistos no sumDir!!! \n");
    return;
  }

  // sumDir->ls();

  TList *sumList = sumDir->GetListOfKeys();
  TIter next(sumList);
  TKey *key;
  printf("addSumHistos %u \n", sumList->GetEntries());
  while (TKey *key = (TKey *)next())
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    std::string name = string(h->GetName());
    if (name.find("QSumChan") == std::string::npos)
      continue;
    string chan = name.substr(name.find_last_of("n") + 1);
    int ichan = stoi(chan);
    if (first)
    {
      cout << " addSumHistos clone " << name << " chan " << ichan << endl;
      TH1D *hclone = (TH1D *)h->Clone(Form("runQSumCh%i", ichan));
      runQSum.push_back(hclone);
      histMap.insert(std::pair<int, TH1D *>(ichan, hclone));
      fout->Add(hclone);
    }
    else
    {
      cout << " addSumHistos add  " << name << " chan " << ichan << endl;
      histMap.at(ichan)->Add(h);
    }
    // name.Form("SumWave-%s-%s", h->GetName(), tag.Data());
    // TH1D *hsave = (TH1D *)h->Clone(name);
    //  rawSumDir->Add(hsave);
  }
  first = false;
}

string currentDate()
{
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  char output[30];
  strftime(output, 30, "%Y-%m-%d-%H-%M", timeinfo);
  return string(output);
}
// count subruns and channels
unsigned long countFiles()
{
  TString dirName = TString("myData/");
  TSystemDirectory dir("mydata", dirName); // TSystemDirectory
  TList *files = dir.GetListOfFiles();     //
  TIter next(files);
  TSystemFile *file;
  while ((file = (TSystemFile *)next()))
  {
    string name = string(file->GetName());
    cout << name << endl;
    string exten = name.substr(name.find_last_of(".") + 1);
    if (exten != string("root"))
      continue;
    if (name.find("anaRun-run") != std::string::npos)
      fileList.push_back(TString(name.c_str()));
  }
  return fileList.size();
}

int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << " make summary plots  " << endl;
  printf(" usage: summary  <max files 0=all>  \n ");
  if (argc < 1)
    exit(0);
  TString tag("run");

  TDatime dopeTime(2023, 3, 9, 22, 0, 0);

  countFiles();
  unsigned nchan = 13;
  vecQsum.resize(nchan);
  vecEQsum.resize(nchan);
  vecQPE.resize(nchan);
  vecEQPE.resize(nchan);
  vSlope.resize(nchan);
  vESlope.resize(nchan);
  // for (unsigned ic = 0; ic < nchan; ++ic)
  //   printf(" chan %u vecQsum %lu \n",ic,vecQsum[ic].size());

  printf(" for %s found %lu files \n", tag.Data(), fileList.size());
  Long64_t maxFiles = fileList.size();
  if (argc > 1)
  {
    maxFiles = atoi(argv[2]);
  }

  std::string sdate = currentDate();
  cout << " starting summary for  " << maxFiles << " files on " << sdate << endl;

  fout = new TFile(Form("summary-%s.root", sdate.c_str()), "recreate");
  foutSumDir = fout->mkdir("sumDir");
  fout->cd();

  hQPEChan = new TH1D("QPEChan", "QPE  by channel", 12, 0, 12);
  hQPEChan->Sumw2();

  for (unsigned ifile = 0; ifile < maxFiles; ++ifile)
  {
    cout << " starting anaRunFile " << fileList[ifile] << endl;
    TString fullName = TString("myData/") + fileList[ifile];
    fin = new TFile(fullName);
    TGraphErrors *gslope;
    fin->GetObject("tbfile", bf);
    fin->GetObject("slope-graph", gslope);
    if (!bf)
    {
      // printf(" no timestamp in file %s \n",fileList[i].Data());
      continue;
    }
    if (!gslope)
    {
      // printf(" no timestamp in file %s \n",fileList[i].Data());
      continue;
    }
    cout << bf->GetTitle() << " modified " << bf->modified << " slope graph " << gslope->GetTitle() << endl;

    // get slopes from graph
    unsigned nslope = gslope->GetN();
    printf(" nslope %u \n", nslope);
    // the graph starts with ichn = 3
    for (unsigned ic = 0; ic < nchan;  ++ic)
    {
      vSlope[ic].push_back(0.);
      vESlope[ic].push_back(0);
    }
    for (unsigned ic = 0; ic < nslope ; ++ic)
    {
      unsigned ilast = vSlope[ic].size()-1;
      unsigned ichan = unsigned(gslope->GetPointX(ic));
      double yval = gslope->GetPointY(ic);
      double yerr = gslope->GetErrorY(ic);
      printf("filling chan %u file %u %f %f \n", ichan, ilast , yval, yerr);

      vSlope[ichan][ilast] = gslope->GetPointY(ic);
      vESlope[ichan][ilast] = gslope->GetErrorY(ic);
    }

    printf("\t slopes %lu \n", vSlope.size());
    for (unsigned j = 0; j < vSlope.size(); ++j)
    {
      for (unsigned k = 0; k < vSlope[j].size(); ++k)
        printf("chan %u file %u %f %f \n", j, k, vSlope[j][k], vESlope[j][k]);
    }

    fin->GetObject("sumDir", sumDir);
    if (!sumDir)
    {
      printf(" sumDir not found in file %s \n", fin->GetName());
      continue;
    }
    printf(" sumDir for file %s \n", fin->GetName());

    // use bf->modified a std string
    /* migrate fit from post.C:        TF1 *g = (TF1*)hHitSum[i]->GetListOfFunctions()->FindObject("expo");*/

    TH1D *hqsum;
    TH1D *hqprompt;
    fin->GetObject("histQsum", hqsum);
    fin->GetObject("histQprompt", hqprompt);

    TH1D *hq = (TH1D *)hqsum->Clone(Form("hqsum%i", ifile));
    TH1D *hp = (TH1D *)hqsum->Clone(Form("hqprompt%i", ifile));
    fout->Add(hq);
    fout->Add(hp);
    eventCount = NULL;
    fin->GetObject("eventCount", eventCount);
    if (eventCount)
    {
      printf(" eventCount for file %s  entries %E \n", fin->GetName(), eventCount->GetEntries());
      TH1D *hevcount = (TH1D *)eventCount->Clone(Form("eventCount%i", ifile));
      fout->Add(hevcount);
    }

    // addSumHistos();
    QPEFits(ifile);

    if (!hqsum || !hqprompt)
      continue;

    cout << "for file  " << ifile << " " << fileList[ifile] << " " << hqsum->GetName() << " " << hqsum->GetNbinsX() << endl;
    cout << "ssssssss  summed hits    " << hSumWave.size() << " , " << hSumHitWave.size() << endl;
    for (unsigned ih = 0; ih<hSumWave.size(); ++ih)
      printf(" %u %s %s \n", ih, hSumWave[ih]->GetName(), hSumHitWave[ih]->GetName());

    filenum.push_back(double(ifile));
    efilenum.push_back(0);
    TDatime datetime = getTime();
    fileDatime.push_back(datetime);
    fileTime.push_back(datetime.Convert());

    cout << " TIME FILE  " << ifile << " " << fileList[ifile] << " modified " << bf->modified << " TDatime " << datetime.AsString() << " as int " << datetime.Convert() << endl;

    if (eventCount)
    {
      for (int i = 0; i <= eventCount->GetNbinsX(); ++i)
      {
        double norm = eventCount->GetBinContent(i + 1) / eventCount->GetBinContent(1);
        printf(" %i count %0f norm %f \n", i - 1, eventCount->GetBinContent(i + 1), norm);
      }
    }

    for (int i = 0; i < hqsum->GetNbinsX() - 1; ++i)
    {
      double norm = 1.0;
      // now overflow is every event and qsum bin is
      if (eventCount)
        norm = eventCount->GetBinContent(i) / eventCount->GetBinContent(0);
      double val = 0;
      double eval = 0;
      if (!isnan(hqsum->GetBinContent(i + 1)) && !isinf(hqsum->GetBinContent(i + 1)))
      {
        val = hqsum->GetBinContent(i + 1) / norm;
        eval = hqsum->GetBinError(i + 1) / norm;
      }
      vecQsum[i].push_back(val);
      vecEQsum[i].push_back(eval);
      cout << "chan " << i + 1 << " qsum " << hqsum->GetBinContent(i + 1) << " norm " << norm << " size " << vecQsum[i].size() << endl;
    }
  } // end loop over files
  printf(" files %lu \n", filenum.size());
  // for (unsigned ic = 0; ic < nchan; ++ic)
  //  printf(" chan %u vecQsum %lu  \n", ic, vecQsum[ic].size());
  for (unsigned jfile = 0; jfile < filenum.size(); ++jfile)
  {
    int ifile = int(filenum[jfile]);
    printf("  summary file %u  %s chan6 %f  \n", ifile, fileList[ifile].Data(), vecQsum[7][ifile]);
  }
  int myColor[13] = {41, 42, 43, 44, 45, 46, 2, 3, 4, 31, 32, 33, 34};
  int myStyle[13] = {21, 22, 23, 24, 25, 26, 21, 22, 23, 31, 32, 33, 34};

  // normalize to first file
  vector<double> normQsum;
  normQsum.resize(vecQsum.size());
  for (unsigned ic = 0; ic < vecQsum.size(); ++ic)
  {
    printf(" vecQsum %i %lu \n", ic, vecQsum[ic].size());
    normQsum[ic] = 1;
    if (vecQsum[ic].size() > 0)
    { // ave over before doping
      double beforeSum = 0;
      int normCount = 0;
      if(vecQsum[ic].size()>20){
        for (unsigned jt = 0; jt < 20; ++jt)
        {
          if (!isinf(vecQsum[ic][jt]) && vecQsum[ic][jt] > 0 && fileDatime[jt].Convert() < dopeTime.Convert())
          {
            beforeSum += vecQsum[ic][jt];
            ++normCount;
          }
      }
      }
      if(normCount>0) normQsum[ic] = beforeSum / double(normCount);
    }
    printf("\t  normQsum =  %f  \n", normQsum[ic]);
  }

  vector<double> normQPE;
  normQPE.resize(vecQPE.size());
  for (unsigned ic = 0; ic < vecQPE.size(); ++ic)
  {
    printf(" vecQPE %i %lu \n", ic, vecQPE[ic].size());
    normQPE[ic] = 1;
    if (vecQPE[ic].size() > 0)
    {
      printf("\t vecQPE =  %f  \n", vecQPE[ic][0]);
      if (!isinf(vecQPE[ic][0]) && vecQPE[ic][0] > 0)
        normQPE[ic] = vecQPE[ic][0];
    }
    printf("\t normQPE =  %f  \n", normQPE[ic]);
  }
  // graphs without norm  one graph per channel
  TString ylabel;
  vector<TGraphErrors *> gqsumUn;
  for (unsigned ic = 0; ic < nchan; ++ic)
  {
    // cout << " add " << ic << endl;
    gqsumUn.push_back(new TGraphErrors(filenum.size(), &fileTime[0], &(vecQsum[ic][0]), &efilenum[0], &(vecEQsum[ic][0])));
    gqsumUn[ic]->SetName(Form("qsumChanUn%i", ic));
    gqsumUn[ic]->SetTitle(Form("qsum-unnormalized-chan-%i", ic));
    gqsumUn[ic]->SetMarkerSize(1);
    gqsumUn[ic]->SetMarkerColor(myColor[ic]);
    gqsumUn[ic]->SetMarkerStyle(myStyle[ic]);
    fout->Add(gqsumUn[ic]);
  }
  ylabel.Form("integrated charge");
  TMultiGraph *mgL1 = new TMultiGraph();
  mgL1->Add(gqsumUn[6]);
  mgL1->Add(gqsumUn[7]);
  mgL1->Add(gqsumUn[8]);
  setTimeGraph(mgL1, ylabel);
  TCanvas *canL1 = new TCanvas(Form("QsummaryL1-%s", sdate.c_str()), Form("QsummaryL1-%s", sdate.c_str()));
  mgL1->Draw("ap");
  gPad->Update();
  canL1->BuildLegend();
  canL1->SetGrid();
  fout->Append(canL1);

  TMultiGraph *mgL2 = new TMultiGraph();
  mgL2->Add(gqsumUn[3]);
  // mgL2->Add(gqsumUn[4]);
  mgL2->Add(gqsumUn[5]);
  setTimeGraph(mgL2, ylabel);
  TCanvas *canL2 = new TCanvas(Form("QsummaryL2-%s", sdate.c_str()), Form("QsummaryL2-%s", sdate.c_str()));
  mgL2->Draw("ap");
  gPad->Update();
  canL2->BuildLegend();
  canL2->SetGrid();
  fout->Append(canL2);

  TMultiGraph *mgL3 = new TMultiGraph();
  mgL3->Add(gqsumUn[0]);
  mgL3->Add(gqsumUn[1]);
  mgL3->Add(gqsumUn[2]);
  setTimeGraph(mgL3, ylabel);
  TCanvas *canL3 = new TCanvas(Form("QsummaryL3-%s", sdate.c_str()), Form("QsummaryL3-%s", sdate.c_str()));
  mgL3->Draw("ap");
  gPad->Update();
  canL3->BuildLegend();
  canL3->SetGrid();
  fout->Append(canL3);

  // apply norms
  for (unsigned ic = 0; ic < vecQsum.size(); ++ic)
  {
    for (unsigned jf = 0; jf < vecQsum[ic].size(); ++jf)
    {
      vecQsum[ic][jf] /= normQsum[ic];
      vecEQsum[ic][jf] /= normQsum[ic];
    }
  }
  // apply norms
  for (unsigned ic = 0; ic < vecQPE.size(); ++ic)
  {
    for (unsigned jf = 0; jf < vecQPE[ic].size(); ++jf)
    {
      vecQPE[ic][jf] /= normQPE[ic];
      vecEQPE[ic][jf] /= normQPE[ic];
    }
  }

  // one graph per channel
  vector<TGraphErrors *> gqsum;
  TMultiGraph *mg = new TMultiGraph();
  for (unsigned ic = 0; ic < nchan; ++ic)
  {
    // cout << " add " << ic << endl;
    gqsum.push_back(new TGraphErrors(filenum.size(), &fileTime[0], &(vecQsum[ic][0]), &efilenum[0], &(vecEQsum[ic][0])));
    gqsum[ic]->SetName(Form("qsumChan%i", ic));
    gqsum[ic]->SetTitle(Form("qsum-chan-%i", ic));
    gqsum[ic]->SetMarkerSize(1);
    gqsum[ic]->SetMarkerColor(myColor[ic]);
    gqsum[ic]->SetMarkerStyle(myStyle[ic]);
    fout->Add(gqsum[ic]);
    if (ic == 6 || ic == 7 || ic == 8)
      mg->Add(gqsum[ic]);
  }
  // overlay all channel graphs on canvas

  // QPE graphs one graph per channel
  vector<TGraphErrors *> gqpe;
  TMultiGraph *mgQPE = new TMultiGraph();
  for (unsigned ic = 0; ic < nchan; ++ic)
  {
    // cout << " add " << ic << endl;
    gqpe.push_back(new TGraphErrors(filenum.size(), &fileTime[0], &(vecQPE[ic][0]), &efilenum[0], &(vecEQPE[ic][0])));
    gqpe[ic]->SetName(Form("QPEChan%i", ic));
    gqpe[ic]->SetTitle(Form("QPE-chan-%i", ic));
    gqpe[ic]->SetMarkerSize(1);
    gqpe[ic]->SetMarkerColor(myColor[ic]);
    gqpe[ic]->SetMarkerStyle(myStyle[ic]);
    fout->Add(gqpe[ic]);
    if (ic == 6 || ic == 7 || ic == 8)
      mgQPE->Add(gqpe[ic]);
  }
  // overlay all channel graphs on canvas
  // mg->GetXaxis()->SetNdivisions(1010);
  /*
    n = n1 + 100 * n2 + 10000 * n3 Where n1 is the number of primary divisions, n2 is the number of second order divisions and n3 is the number of third order divisions. n < 0, the axis will be forced to use exactly n divisions.
  */
  int ndiv = 10 + 100 * 5 + 10000 * 3;
  ylabel.Form("integrated charge (normed) ");
  setTimeGraph(mg, ylabel);
  TCanvas *can = new TCanvas(Form("Qsummary-%s", sdate.c_str()), Form("Qsummary-%s", sdate.c_str()));
  mg->Draw("ap");
  gPad->Update();
  can->BuildLegend();
  can->SetGrid();
  can->Print(".png");
  fout->Append(can);

  ylabel.Form("single photon charge (normed)");
  setTimeGraph(mgQPE, ylabel);
  TCanvas *canqpe = new TCanvas("QPE", "QPE");
  mgQPE->Draw("ap");
  canqpe->BuildLegend();
  canqpe->SetGrid();
  fout->Append(canqpe);

  // slope graphs
  printf(" \t\t make slope graph %lu \n", filenum.size());
  // one graph per channel
  vector<TGraphErrors *> graphSlope;
  TMultiGraph *mgslope = new TMultiGraph();
  for (unsigned ic = 6; ic < 9; ++ic)
  {
    cout << " add " << ic << " size " << vSlope[ic].size() << " size " << vESlope[ic].size() << endl;
    cout << "     " << ic << " size " << fileTime.size() << " size " << efilenum.size() << endl;
    graphSlope.push_back(new TGraphErrors(filenum.size(), &fileTime[0], &(vSlope[ic][0]), &efilenum[0], &(vESlope[ic][0])));
    unsigned ilast = graphSlope.size() - 1;
    graphSlope[ilast]->SetName(Form("slopehan%i", ic));
    graphSlope[ilast]->SetTitle(Form("qslope-chan-%i", ic));
    graphSlope[ilast]->SetMarkerSize(1);
    graphSlope[ilast]->SetMarkerColor(myColor[ic]);
    graphSlope[ilast]->SetMarkerStyle(myStyle[ic]);
    fout->Add(graphSlope[ilast]);
    mgslope->Add(graphSlope[ilast]);
  }
  ylabel.Form(" slope fit value (?) ");
  setTimeGraph(mgslope, ylabel);
  TCanvas *canSlope = new TCanvas(Form("Slope-%s", sdate.c_str()), Form("Slope-%s", sdate.c_str()));
  mgslope->Draw("ap");
  gPad->Update();
  canSlope->BuildLegend();
  canSlope->SetGrid();
  fout->Append(canSlope);

  // finish
  //fout->ls();
  fout->Write();

  // report
  for (unsigned it = 0; it < fileTime.size(); ++it)
    if (fileDatime[it].Convert() < dopeTime.Convert())
      cout << " before " << fileDatime[it].AsString() << "  " << fileList[it] << endl;
    else
      cout << " after  " << fileDatime[it].AsString() << "  " << fileList[it] << endl;
  cout << "summary finished " << maxFiles << " " << fout->GetName() << endl;
  exit(0);
}