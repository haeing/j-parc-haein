#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "TRandom3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;

namespace Masses {
  const double Kminus    = 0.493677;
  const double Proton    = 0.9382720813;
  const double Neutron   = 0.9395654133;
  const double Pion      = 0.13957039;
  const double Pi0       = 0.1349768;
  const double SigmaMinus= 1.197449;
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
              double& out,
              int maxTry = 200000)
{
  if (maxv <= minv) return false;

  for (int i = 0; i < maxTry; ++i) {
    double x = rng.BreitWigner(mean, width);
    if (x > minv && x < maxv) {
      out = x;
      return true;
    }
  }
  return false;
}

bool GenerateDecay(TGenPhaseSpace& event,
                   const TLorentzVector& mother,
                   const vector<double>& masses,
                   double& weight)
{
  vector<double> m = masses;
  TLorentzVector motherCopy = mother;
  if (!event.SetDecay(motherCopy, (int)m.size(), m.data())) return false;
  weight = event.Generate();
  return true;
}

bool AllChargedPionsDetected(const vector<TLorentzVector*>& chargedPions,
                             double blindMaxDeg = 8.0)
{
  for (auto tr : chargedPions) {
    double th = ThetaDeg(*tr);
    if (th >= 0.0 && th <= blindMaxDeg) return false;
  }
  return true;
}

void acceptance_sigma1670_npipipipi(Long64_t nEvent = 500000,
                                    double pBeam = 1.10,
                                    double blindMaxDeg = 10.0)
{
  gStyle->SetOptStat(0);

  const double mK   = Masses::Kminus;
  const double mp   = Masses::Proton;
  const double mpi  = Masses::Pion;
  const double mpi0 = Masses::Pi0;
  const double mn   = Masses::Neutron;
  const double mSm  = Masses::SigmaMinus;

  const double eBeam = std::sqrt(pBeam * pBeam + mK * mK);

  TLorentzVector beam_lab(0.0, 0.0, pBeam, eBeam);
  TLorentzVector target_lab(0.0, 0.0, 0.0, mp);
  TLorentzVector initial_lab = beam_lab + target_lab;

  TVector3 beta_cm = initial_lab.BoostVector();
  TRandom3 rng(0);

  cout << "=====================================\n";
  cout << " Reaction 2 neutron branch\n";
  cout << " K- p -> Sigma(1670)+ pi-_prod\n";
  cout << " Sigma(1670)+ -> Lambda(1405)^0 pi+ pi0\n";
  cout << " Lambda(1405)^0 -> Sigma- pi+\n";
  cout << " Sigma- -> n pi-\n";
  cout << "-------------------------------------\n";
  cout << " Beam momentum = " << pBeam << " GeV/c\n";
  cout << " sqrt(s)       = " << initial_lab.M() << " GeV\n";
  cout << " Blind cone    = 0 - " << blindMaxDeg << " deg\n";
  cout << " Acceptance    = all 4 charged pions outside blind cone\n";
  cout << "=====================================\n";

  // generated / accepted vs cos(theta_cm(Sigma1670))
  TH1D* hGen = new TH1D("hGen",
    "Generated;cos#theta_{cm}(#Sigma(1670));Generated events",
    50, -1.0, 1.0);

  TH1D* hAcc = new TH1D("hAcc",
    "Accepted;cos#theta_{cm}(#Sigma(1670));Accepted events",
    50, -1.0, 1.0);

  // diagnostic: pion lab-angle distributions
  TH1D* hTheta_piProd  = new TH1D("hTheta_piProd",
    "#pi^{-}_{prod};#theta_{lab} [deg];Events", 180, 0, 180);
  TH1D* hTheta_piSstar = new TH1D("hTheta_piSstar",
    "#pi^{+}_{#Sigma^{*}};#theta_{lab} [deg];Events", 180, 0, 180);
  TH1D* hTheta_piL1405 = new TH1D("hTheta_piL1405",
    "#pi^{+}_{#Lambda(1405)};#theta_{lab} [deg];Events", 180, 0, 180);
  TH1D* hTheta_piSm    = new TH1D("hTheta_piSm",
    "#pi^{-}_{#Sigma^{-}};#theta_{lab} [deg];Events", 180, 0, 180);

  // correlation: proton was replaced by neutron branch, so look at one charged pion if needed
  TH2D* hSigmaCos_vs_piProdTheta = new TH2D(
    "hSigmaCos_vs_piProdTheta",
    ";#theta_{lab}(#pi^{-}_{prod}) [deg];cos#theta_{cm}(#Sigma(1670))",
    180, 0.0, 180.0, 80, -1.0, 1.0
  );

  TH1D* hRejectReason = new TH1D(
    "hRejectReason",
    "Reject reason;particle entering 0-8 deg;counts",
    4, 0.5, 4.5
  );
  hRejectReason->GetXaxis()->SetBinLabel(1, "pi-_prod");
  hRejectReason->GetXaxis()->SetBinLabel(2, "pi+_Sigma*");
  hRejectReason->GetXaxis()->SetBinLabel(3, "pi+_L1405");
  hRejectReason->GetXaxis()->SetBinLabel(4, "pi-_Sigma-");

  TGenPhaseSpace prod, dec1, dec2, dec3;

  Long64_t nAccepted = 0;
  Long64_t nGenerated = 0;

  for (Long64_t iev = 0; iev < nEvent; ++iev) {

    // -----------------------------------------
    // K- p -> Sigma(1670)+ + pi-_prod
    // -----------------------------------------
    const double mSigmaMin = Res::Lambda1405_M + mpi + mpi0; // rough minimal allowed chain
    const double mSigmaMax = initial_lab.M() - mpi;

    double mSigma = -1.0;
    if (!SampleBW(rng, Res::Sigma1670_M, Res::Sigma1670_W,
                  mSigmaMin, mSigmaMax, mSigma)) {
      continue;
    }

    double wProd = 0.0;
    if (!GenerateDecay(prod, initial_lab, {mSigma, mpi}, wProd)) continue;

    TLorentzVector* pSigmaStar = prod.GetDecay(0);
    TLorentzVector* pPiProd    = prod.GetDecay(1);

    TLorentzVector sigma_cm = *pSigmaStar;
    sigma_cm.Boost(-beta_cm);
    double coscm = sigma_cm.Vect().CosTheta();

    hGen->Fill(coscm);
    nGenerated++;

    // -----------------------------------------
    // Sigma(1670)+ -> Lambda(1405)^0 + pi+ + pi0
    // -----------------------------------------
    const double mLMin = mSm + mpi;          // Lambda1405 -> Sigma- pi+
    const double mLMax = mSigma - mpi - mpi0;

    double mL1405 = -1.0;
    if (!SampleBW(rng, Res::Lambda1405_M, Res::Lambda1405_W,
                  mLMin, mLMax, mL1405)) {
      continue;
    }

    double w1 = 0.0;
    if (!GenerateDecay(dec1, *pSigmaStar, {mL1405, mpi, mpi0}, w1)) continue;

    TLorentzVector* pL1405   = dec1.GetDecay(0);
    TLorentzVector* pPiSstar = dec1.GetDecay(1); // pi+ from Sigma(1670)

    // -----------------------------------------
    // Lambda(1405)^0 -> Sigma- + pi+
    // -----------------------------------------
    double w2 = 0.0;
    if (!GenerateDecay(dec2, *pL1405, {mSm, mpi}, w2)) continue;

    TLorentzVector* pSigmaMinus = dec2.GetDecay(0);
    TLorentzVector* pPiL1405    = dec2.GetDecay(1); // pi+ from Lambda(1405)

    // -----------------------------------------
    // Sigma- -> n + pi-
    // -----------------------------------------
    double w3 = 0.0;
    if (!GenerateDecay(dec3, *pSigmaMinus, {mn, mpi}, w3)) continue;

    TLorentzVector* pNeutron = dec3.GetDecay(0);
    TLorentzVector* pPiSm    = dec3.GetDecay(1);   // pi- from Sigma-

    (void)pNeutron; // intentionally unused in acceptance

    // diagnostics
    double th_piProd  = ThetaDeg(*pPiProd);
    double th_piSstar = ThetaDeg(*pPiSstar);
    double th_piL1405 = ThetaDeg(*pPiL1405);
    double th_piSm    = ThetaDeg(*pPiSm);

    hTheta_piProd->Fill(th_piProd);
    hTheta_piSstar->Fill(th_piSstar);
    hTheta_piL1405->Fill(th_piL1405);
    hTheta_piSm->Fill(th_piSm);

    hSigmaCos_vs_piProdTheta->Fill(th_piProd, coscm);

    vector<TLorentzVector*> chargedPions = {
      pPiProd,
      pPiSstar,
      pPiL1405,
      pPiSm
    };

    bool accepted = AllChargedPionsDetected(chargedPions, blindMaxDeg);

    if (!accepted) {
      if (th_piProd  <= blindMaxDeg) hRejectReason->Fill(1);
      if (th_piSstar <= blindMaxDeg) hRejectReason->Fill(2);
      if (th_piL1405 <= blindMaxDeg) hRejectReason->Fill(3);
      if (th_piSm    <= blindMaxDeg) hRejectReason->Fill(4);
    }

    if (accepted) {
      hAcc->Fill(coscm);
      nAccepted++;
    }
  }

  hGen->Sumw2();
  hAcc->Sumw2();

  TH1D* hEff = (TH1D*)hAcc->Clone("hEff");
  hEff->SetTitle("Acceptance for n#pi#pi#pi#pi branch;cos#theta_{cm}(#Sigma(1670));Acceptance");
  hEff->Divide(hAcc, hGen, 1.0, 1.0, "B");

  // plotting
  TCanvas* c1 = new TCanvas("c1", "acceptance", 900, 700);
  hEff->SetStats(0);
  hEff->SetMinimum(0.0);
  hEff->SetMaximum(1.05);
  hEff->SetLineWidth(2);
  hEff->Draw("E1");
  c1->SaveAs("acceptance_sigma1670_npipipipi.pdf");

  TCanvas* c2 = new TCanvas("c2", "charged pion angles", 1200, 900);
  c2->Divide(2,2);
  c2->cd(1); hTheta_piProd->Draw();
  c2->cd(2); hTheta_piSstar->Draw();
  c2->cd(3); hTheta_piL1405->Draw();
  c2->cd(4); hTheta_piSm->Draw();
  c2->SaveAs("theta_distributions_sigma1670_npipipipi.pdf");

  TCanvas* c3 = new TCanvas("c3", "reject reason", 800, 700);
  hRejectReason->SetStats(0);
  hRejectReason->SetFillColor(kAzure-9);
  hRejectReason->Draw("hist");
  c3->SaveAs("reject_reason_sigma1670_npipipipi.pdf");

  TCanvas* c4 = new TCanvas("c4", "correlation", 900, 700);
  hSigmaCos_vs_piProdTheta->SetStats(0);
  hSigmaCos_vs_piProdTheta->Draw("colz");
  c4->SaveAs("sigmaCos_vs_piProdTheta_npipipipi.pdf");

  // save root file
  TFile* fout = new TFile("acceptance_sigma1670_npipipipi.root", "RECREATE");
  hGen->Write();
  hAcc->Write();
  hEff->Write();
  hTheta_piProd->Write();
  hTheta_piSstar->Write();
  hTheta_piL1405->Write();
  hTheta_piSm->Write();
  hRejectReason->Write();
  hSigmaCos_vs_piProdTheta->Write();
  fout->Close();

  cout << "\nGenerated events = " << nGenerated << "\n";
  cout << "Accepted events  = " << nAccepted << "\n";
  if (nGenerated > 0) {
    cout << "Overall acceptance = " << (double)nAccepted / (double)nGenerated << "\n";
  }

  cout << "\nSaved files:\n";
  cout << "  acceptance_sigma1670_npipipipi.pdf\n";
  cout << "  theta_distributions_sigma1670_npipipipi.pdf\n";
  cout << "  reject_reason_sigma1670_npipipipi.pdf\n";
  cout << "  sigmaCos_vs_piProdTheta_npipipipi.pdf\n";
  cout << "  acceptance_sigma1670_npipipipi.root\n";
}
