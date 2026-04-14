#include <iostream>
#include <vector>
#include <cmath>

#include "TRandom3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"

using namespace std;

namespace Masses {
  const double Kminus = 0.493677;
  const double Proton = 0.9382720813;
  const double Pion   = 0.13957039;
  const double Pi0    = 0.1349768;
  const double Gamma  = 0.0;
  const double Lambda = 1.115683;
  const double Sigma0 = 1.192642;
  const double SigmaP = 1.18937;
}

namespace Res {
  const double Sigma1670_M = 1.670;
  const double Sigma1670_W = 0.060;
  const double Lambda1405_M = 1.405;
  const double Lambda1405_W = 0.050;
}

double ThetaDeg(const TLorentzVector& p4)
{
  return p4.Vect().Theta() * TMath::RadToDeg();
}

bool SampleBW(TRandom3& rng,
              double mean,
              double width,
              double minv,
              double maxv,
              double& out)
{
  if (maxv <= minv) return false;

  for (int i = 0; i < 200000; ++i) {
    double x = rng.BreitWigner(mean, width);
    if (x > minv && x < maxv) {
      out = x;
      return true;
    }
  }
  return false;
}

// accepted only if ALL charged particles are outside blind cone
bool AllChargedDetected(const vector<TLorentzVector*>& tracks,
                        double blindMax = 10.0)
{
  for (auto tr : tracks) {
    double th = ThetaDeg(*tr);
    if (th >= 0.0 && th <= blindMax) return false;
  }
  return true;
}

void acceptance_vs_cm_sigma(Long64_t nEvent = 500000,
                            double pBeam = 1.10,
                            double blindMax = 21.0)
{
  const double mK   = Masses::Kminus;
  const double mp   = Masses::Proton;
  const double mpi  = Masses::Pion;
  const double mpi0 = Masses::Pi0;

  const double eBeam = sqrt(pBeam * pBeam + mK * mK);

  TLorentzVector beam_lab(0, 0, pBeam, eBeam);
  TLorentzVector target_lab(0, 0, 0, mp);
  TLorentzVector initial_lab = beam_lab + target_lab;

  TVector3 beta_cm = initial_lab.BoostVector();

  cout << "Beam momentum = " << pBeam << " GeV/c\n";
  cout << "sqrt(s) = " << initial_lab.M() << " GeV\n";
  cout << "Blind cone = 0 - " << blindMax << " deg\n";

  TRandom3 rng(0);

  TH1D* hGen_R1 = new TH1D("hGen_R1",
    "Reaction1 generated;cos#theta_{cm}(#Sigma(1670));Generated events", 40, -1, 1);
  TH1D* hAcc_R1 = new TH1D("hAcc_R1",
    "Reaction1 accepted;cos#theta_{cm}(#Sigma(1670));Accepted events", 40, -1, 1);
  

  TH1D* hGen_R2 = new TH1D("hGen_R2",
    "Reaction2 generated;cos#theta_{cm}(#Sigma(1670));Generated events", 40, -1, 1);
  TH1D* hAcc_R2 = new TH1D("hAcc_R2",
    "Reaction2 accepted;cos#theta_{cm}(#Sigma(1670));Accepted events", 40, -1, 1);
  TH1D* hAcc_R2_update = new TH1D("hAcc_R2_update",
			   "Reaction2 accepted;cos#theta_{cm}(#Sigma(1670));Accepted events", 40, -1, 1);
  TH2D* hSigmaVsThetaP_R1 = new TH2D(
  "hSigmaVsThetaP_R1",
  "Reaction1;#theta_{lab}(p) [deg];cos#theta_{cm}(#Sigma(1670))",
  180, 0.0, 180.0,
  80, -1.0, 1.0
				     );

  TH2D* hSigmaCos_vs_ProtonTheta_R1 = new TH2D(
  "hSigmaCos_vs_ProtonTheta_R1",
  "Reaction1;#theta_{lab}(p) [deg];cos#theta_{cm}(#Sigma(1670))",
  180, 0.0, 180.0,
  80, -1.0, 1.0
);

TH2D* hSigmaCos_vs_ProtonTheta_R2 = new TH2D(
  "hSigmaCos_vs_ProtonTheta_R2",
  "Reaction2;#theta_{lab}(p) [deg];cos#theta_{cm}(#Sigma(1670))",
  180, 0.0, 180.0,
  80, -1.0, 1.0
);
  
  TGenPhaseSpace prod, dec1, dec2, dec3;

  for (Long64_t iev = 0; iev < nEvent; ++iev) {

    // =========================================================
    // Reaction 1
    // K- p -> Sigma(1670)+ pi-
    // Sigma(1670)+ -> Sigma0 pi+ pi0
    // Sigma0 -> Lambda gamma
    // Lambda -> p pi-
    // charged final: p, pi-_prod, pi+, pi-_Lambda
    // =========================================================
    {
      double mSigmaMin = Masses::Sigma0 + mpi + mpi0;
      double mSigmaMax = initial_lab.M() - mpi;
      double mSigma = -1.0;

      if (SampleBW(rng, Res::Sigma1670_M, Res::Sigma1670_W,
                   mSigmaMin, mSigmaMax, mSigma)) {

        double massesProd[2] = {mSigma, mpi};
        TLorentzVector initCopy = initial_lab;

        if (prod.SetDecay(initCopy, 2, massesProd)) {
          prod.Generate();

          TLorentzVector* sigma_lab = prod.GetDecay(0);
          TLorentzVector* pi_prod   = prod.GetDecay(1);

          TLorentzVector sigma_cm = *sigma_lab;
          sigma_cm.Boost(-beta_cm);
          double coscm = sigma_cm.Vect().CosTheta();

          hGen_R1->Fill(coscm);

          double masses1[3] = {Masses::Sigma0, mpi, mpi0};
          TLorentzVector sigmaCopy = *sigma_lab;

          if (dec1.SetDecay(sigmaCopy, 3, masses1)) {
            dec1.Generate();

            TLorentzVector* sigma0 = dec1.GetDecay(0);
            TLorentzVector* piPlus = dec1.GetDecay(1);

            double masses2[2] = {Masses::Lambda, Masses::Gamma};
            TLorentzVector sigma0Copy = *sigma0;

            if (dec2.SetDecay(sigma0Copy, 2, masses2)) {
              dec2.Generate();

              TLorentzVector* lambda = dec2.GetDecay(0);

              double masses3[2] = {mp, mpi};
              TLorentzVector lambdaCopy = *lambda;

              if (dec3.SetDecay(lambdaCopy, 2, masses3)) {
                dec3.Generate();

                TLorentzVector* proton = dec3.GetDecay(0);
                TLorentzVector* piLam  = dec3.GetDecay(1);
		double theta_p_lab = ThetaDeg(*proton);
		hSigmaCos_vs_ProtonTheta_R1->Fill(theta_p_lab, coscm);

                vector<TLorentzVector*> tracks = {
                  proton,
                  pi_prod,
                  piPlus,
                  piLam
                };

                if (AllChargedDetected(tracks, blindMax)) {
                  hAcc_R1->Fill(coscm);
                }
		
              }
            }
          }
        }
      }
    }

    // =========================================================
    // Reaction 2
    // K- p -> Sigma(1670)+ pi-
    // Sigma(1670)+ -> Lambda(1405) pi+ pi0
    // Lambda(1405) -> Sigma+ pi-
    // Sigma+ -> p pi0
    // charged final: p, pi-_prod, pi+, pi-_L1405
    // =========================================================
    {
      double mSigmaMin = Masses::SigmaP + mpi + mpi + mpi0;
      double mSigmaMax = initial_lab.M() - mpi;
      double mSigma = -1.0;

      if (SampleBW(rng, Res::Sigma1670_M, Res::Sigma1670_W,
                   mSigmaMin, mSigmaMax, mSigma)) {

        double massesProd[2] = {mSigma, mpi};
        TLorentzVector initCopy = initial_lab;

        if (prod.SetDecay(initCopy, 2, massesProd)) {
          prod.Generate();

          TLorentzVector* sigma_lab = prod.GetDecay(0);
          TLorentzVector* pi_prod   = prod.GetDecay(1);

          TLorentzVector sigma_cm = *sigma_lab;
          sigma_cm.Boost(-beta_cm);
          double coscm = sigma_cm.Vect().CosTheta();

          hGen_R2->Fill(coscm);

          double mLMin = Masses::SigmaP + mpi;
          double mLMax = mSigma - mpi - mpi0;
          double mL = -1.0;

          if (SampleBW(rng, Res::Lambda1405_M, Res::Lambda1405_W,
                       mLMin, mLMax, mL)) {

            double masses1[3] = {mL, mpi, mpi0};
            TLorentzVector sigmaCopy = *sigma_lab;

            if (dec1.SetDecay(sigmaCopy, 3, masses1)) {
              dec1.Generate();

              TLorentzVector* lambda1405 = dec1.GetDecay(0);
              TLorentzVector* piPlus     = dec1.GetDecay(1);

              double masses2[2] = {Masses::SigmaP, mpi};
              TLorentzVector lCopy = *lambda1405;

              if (dec2.SetDecay(lCopy, 2, masses2)) {
                dec2.Generate();

                TLorentzVector* sigmaP   = dec2.GetDecay(0);
                TLorentzVector* piL1405  = dec2.GetDecay(1);

                double masses3[2] = {mp, mpi0};
                TLorentzVector sigmaPCopy = *sigmaP;

                if (dec3.SetDecay(sigmaPCopy, 2, masses3)) {
                  dec3.Generate();

                  TLorentzVector* proton = dec3.GetDecay(0);

                  vector<TLorentzVector*> tracks = {
                    proton,
                    pi_prod,
                    piPlus,
                    piL1405
                  };

                  if (AllChargedDetected(tracks, blindMax)) {
                    hAcc_R2->Fill(coscm);
                  }
		  if(AllChargedDetected(tracks, 5)){
		    hAcc_R2_update->Fill(coscm);
		  }
                }
              }
            }
          }
        }
      }
    }
  }

  hGen_R1->Sumw2();
  hAcc_R1->Sumw2();
  hGen_R2->Sumw2();
  hAcc_R2->Sumw2();
  hAcc_R2_update->Sumw2();

  TH1D* hEff_R1 = (TH1D*)hAcc_R1->Clone("hEff_R1");
  hEff_R1->SetTitle("Reaction1 acceptance;cos#theta^{cm}_{#Sigma(1670)};Acceptance");
  hEff_R1->Divide(hAcc_R1, hGen_R1, 1.0, 1.0, "B");

  TH1D* hEff_R2 = (TH1D*)hAcc_R2->Clone("hEff_R2");
  hEff_R2->SetTitle("Reaction2 acceptance;cos#theta^{cm}_{#Sigma(1670)};Acceptance");
  hEff_R2->Divide(hAcc_R2, hGen_R2, 1.0, 1.0, "B");

  TH1D* hEff_R2_update = (TH1D*)hAcc_R2_update->Clone("hEff_R2_update");
  hEff_R2_update->SetTitle("Reaction2 acceptance;cos#theta^{cm}_{#Sigma(1670)};Acceptance");
  hEff_R2_update->Divide(hAcc_R2_update, hGen_R2, 1.0, 1.0, "B");

  TCanvas* c = new TCanvas("cAcc", "acceptance", 1200, 500);
  //c->Divide(2, 1);

  c->cd(1);
  /*
  hEff_R1->SetStats(0);
  hEff_R1->SetMinimum(0.0);
  hEff_R1->SetMaximum(1.05);
  hEff_R1->SetLineWidth(2);
  hEff_R1->Draw("E1");
  */
  

  //c->cd(2);
  hEff_R2->SetStats(0);
  hEff_R2->SetMinimum(0.0);
  hEff_R2->SetMaximum(1.05);
  hEff_R2->SetLineWidth(2);
  hEff_R2->SetLineColor(kRed);
  hEff_R2->SetMarkerColor(kRed);
  hEff_R2->Draw("E1");

  hEff_R2_update->SetStats(0);
  hEff_R2_update->SetMinimum(0.0);
  hEff_R2_update->SetMaximum(1.05);
  hEff_R2_update->SetLineWidth(2);
  hEff_R2_update->Draw("E1 same");

  TLegend* leg = new TLegend(0.65, 0.7, 0.88, 0.88);
  //leg->AddEntry(hEff_R1, "#Sigma#pi#pi", "l");
  leg->AddEntry(hEff_R2, "CDS System", "l");
  leg->AddEntry(hEff_R2_update, "CDS with Forward detector systems", "l");
  leg->Draw("same");
  

  c->SaveAs("acceptance_vs_cm_sigma_allChargedOutside8deg.pdf");

  TFile* fout = new TFile("acceptance_vs_cm_sigma_allChargedOutside8deg.root", "RECREATE");
  hGen_R1->Write();
  hAcc_R1->Write();
  hEff_R1->Write();
  hGen_R2->Write();
  hAcc_R2->Write();
  hEff_R2->Write();
  fout->Close();

  auto c1 = new TCanvas("c1","c1");
  c1->cd();
  hSigmaCos_vs_ProtonTheta_R1->Draw();

  cout << "Saved:\n";
  cout << "  acceptance_vs_cm_sigma_allChargedOutside8deg.pdf\n";
  cout << "  acceptance_vs_cm_sigma_allChargedOutside8deg.root\n";
}
