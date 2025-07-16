#define Steps_cxx
#include <iostream>
#include <fstream>
#include <numeric>
#include "TMath.h"
#include "TF1.h"
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "Steps.h"

TH3D *hOriginMap;
std::vector<TH3D *> hFlux;
TFile *fout;

void geant(Long64_t maxEntry = 0)
{
   /* output file */
   time_t rawtime;
   struct tm *timeinfo;
   time(&rawtime);
   timeinfo = localtime(&rawtime);
   char output[30];
   strftime(output, 30, "%Y-%m-%d-%H-%M", timeinfo);
   TString tdateTag = TString(output);
   TString fullname = (Form("geantSim-%s-%lld.root", tdateTag.Data(), maxEntry));
   fout = new TFile(fullname, "recreate"); // DEF made to update rather than recreate so that it doesn't write over a file already mad

   // bin sizes 3D hists have same number of bins
   int nbinx = 5;
   int nbiny = 5;
   int nbinz = 5;

   hFlux.resize(3);

   hOriginMap = new TH3D("OriginMap", "origin map", nbinx, 0., 20., nbiny, -1., 1., nbinz, -TMath::Pi(), TMath::Pi());
   hOriginMap->GetXaxis()->SetTitle(" radius [cm]");
   hOriginMap->GetYaxis()->SetTitle(" cos(theta)");
   hOriginMap->GetZaxis()->SetTitle(" phi ");

   // flux maps
   for (int i = 0; i < hFlux.size(); ++i)
   {
      hFlux[i] = new TH3D(Form("FluxMapChan%i", i + 9), Form("FluxMapChan%i", i + 9), nbinx, 0., 20., nbiny, -1., 1., nbinz, -TMath::Pi(), TMath::Pi());
      hFlux[i]->GetXaxis()->SetTitle(" radius [cm]");
      hFlux[i]->GetYaxis()->SetTitle(" cos(theta)");
      hFlux[i]->GetZaxis()->SetTitle(" phi ");
   }

   // instantiate Steps class
   Steps *stp = new Steps;

   /* loop over entris */
   Long64_t nentries = stp->fChain->GetEntriesFast();

   if (maxEntry > 0)
      nentries = maxEntry;

   printf("Steps::Loop over %lld \n", nentries);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = stp->LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = stp->fChain->GetEntry(jentry);
      nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      /* cut bad entries */
      if (stp->position_phi_rad == 0)
         continue;

      /* get event orign */
      Double_t radius = stp->position_r_mm / 10.; // convert mm to cm
      Double_t cost = cos(stp->position_theta_rad);
      Double_t phi = stp->position_phi_rad;
      hOriginMap->Fill(radius, cost, phi, 1.);

      // ensure line of sight for flux fraction
      if (stp->TriggerSiPM1_LOS)
         hFlux[0]->Fill(radius, cost, phi, stp->TriggerSiPM1_fluxfraction);
      if (stp->TriggerSiPM2_LOS)
         hFlux[1]->Fill(radius, cost, phi, stp->TriggerSiPM2_fluxfraction);
      if (stp->TriggerSiPM3_LOS)
         hFlux[2]->Fill(radius, cost, phi, stp->TriggerSiPM3_fluxfraction);
   }

   double probabilityNorm = hOriginMap->GetEntries();
   printf(" norm to get probability %.0f \n", probabilityNorm);

   /* normalize */
   for (int i = 0; i < hFlux.size(); ++i) // loop over sipms
   {
      for (int xbin = 0; xbin < hOriginMap->GetNbinsX(); ++xbin)
      {
         for (int ybin = 0; ybin < hOriginMap->GetNbinsY(); ++ybin)
         {
            for (int zbin = 0; zbin < hOriginMap->GetNbinsZ(); ++zbin)
            {
               double norm = hOriginMap->GetBinContent(xbin, ybin, zbin);
               // origin norm to probability
               hOriginMap->SetBinContent(xbin, ybin, zbin, norm / probabilityNorm);
               if (norm > 0)
               {
                  printf("Origin chan %i (%i,%i,%i) nevents %E prob %E \n", i, xbin, ybin, zbin, norm, norm / probabilityNorm);
                  // flux norm
                  double val = hFlux[i]->GetBinContent(xbin, ybin, zbin);
                  hFlux[i]->SetBinContent(xbin, ybin, zbin, val / norm); // number of times bin was filled
                  printf("chan %i (%i,%i,%i) val %E norm %E ave %E \n", i, xbin, ybin, zbin, val, norm, val / norm);
               }
            }
         }
      }
   }

   /* check with projections of all bins */
   TH1D *projRadius = hOriginMap->ProjectionX("OriginRadius", 0, hOriginMap->GetNbinsY(), hOriginMap->GetNbinsZ());
   TH1D *projCos = hOriginMap->ProjectionY("OriginCos", 0, hOriginMap->GetNbinsX(), hOriginMap->GetNbinsZ());
   TH1D *projPhi = hOriginMap->ProjectionZ("OriginPhi", 0, hOriginMap->GetNbinsX(), hOriginMap->GetNbinsY());

   fout->Add(projRadius);
   fout->Add(projCos);
   fout->Add(projPhi);

   TCanvas *canx = new TCanvas("EventOrigins", "EventOrigins");
   canx->Divide(1, 3);
   canx->cd(1);
   gPad->SetLogy();
   projRadius->Draw();
   canx->cd(2);
   gPad->SetLogy();
   projCos->Draw();
   canx->cd(3);
   gPad->SetLogy();
   projPhi->Draw();
   canx->Print(".pdf");

   /* flux maps */
   std::vector<TCanvas *> fluxCan;
   fluxCan.resize(3);
   TString canName;
   for (int i = 0; i < hFlux.size(); ++i)
   { // loop over sipms
      TH1D *pRadius = hFlux[i]->ProjectionX(Form("FluxRadiusChan%i", i + 9), 0, hFlux[i]->GetNbinsY(), hFlux[i]->GetNbinsZ());
      TH1D *pCos = hFlux[i]->ProjectionY(Form("FluxCosChan%i", i + 9), 0, hFlux[i]->GetNbinsX(), hFlux[i]->GetNbinsZ());
      TH1D *pPhi = hFlux[i]->ProjectionZ(Form("FluxPhihan%i", i + 9), 0, hFlux[i]->GetNbinsX(), hFlux[i]->GetNbinsY());
      fout->Add(pRadius);
      fout->Add(pCos);
      fout->Add(pPhi);

      canName.Form("FluxChan%i", i + 9);
      fluxCan[i] = new TCanvas(canName, canName);
      fluxCan[i]->Divide(1, 3);
      fluxCan[i]->cd(1);
      gPad->SetLogy();
      pRadius->Draw();
      fluxCan[i]->cd(2);
      gPad->SetLogy();
      pCos->Draw();
      fluxCan[i]->cd(3);
      gPad->SetLogy();
      pPhi->Draw();
      fluxCan[i]->Print(".pdf");
   }

   fout->ls();
   fout->Write();
   printf("wrote file %s\n", fout->GetName());
}
