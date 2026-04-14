#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "TRandom3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

struct Resonance {
  string name;
  double m0;
  double gamma;
};

struct ParticleHists {
  string label;
  TH1D* hTheta = nullptr;
  TH1D* hCos   = nullptr;
  TH2D* hPTheta = nullptr;
};

struct GroupHists {
  string label;
  TH1D* hTheta = nullptr;
  TH1D* hCos   = nullptr;
  TH2D* hPTheta = nullptr;
};

struct ReactionHists {
  string tag;

  // individual charged particles
  vector<ParticleHists> charged;

  // grouped histograms
  GroupHists protonGroup;
  GroupHists pionGroup;

  // sampled masses
  TH1D* hMParent = nullptr;
  TH1D* hMSub    = nullptr;
};

namespace Masses {
  const double Kminus = 0.493677;
  const double Proton = 0.9382720813;
  const double Pion   = 0.13957039;
  const double Pi0    = 0.1349768;
  const double Lambda = 1.115683;
  const double SigmaP = 1.18937;
}

namespace Res {
  const Resonance Sigma1670  = {"Sigma1670", 1.670, 0.060};
  const Resonance Lambda1405 = {"Lambda1405", 1.405, 0.050};
}

double ThetaDeg(const TLorentzVector& p4) {
  return p4.Vect().Theta() * TMath::RadToDeg();
}

ParticleHists BookParticle(const string& tag, const string& name) {
  ParticleHists h;
  h.label = name;

  h.hTheta = new TH1D(
    Form("hTheta_%s_%s", tag.c_str(), name.c_str()),
    Form("%s;#theta_{lab} [deg];Weighted events", name.c_str()),
    180, 0.0, 180.0
  );

  h.hCos = new TH1D(
    Form("hCos_%s_%s", tag.c_str(), name.c_str()),
    Form("%s;cos#theta_{lab};Weighted events", name.c_str()),
    200, -1.0, 1.0
  );

  h.hPTheta = new TH2D(
    Form("hPTheta_%s_%s", tag.c_str(), name.c_str()),
    Form("%s;#theta_{lab} [deg];p_{lab} [GeV/c]", name.c_str()),
    180, 0.0, 180.0, 220, 0.0, 2.2
  );

  return h;
}

GroupHists BookGroup(const string& tag, const string& name) {
  GroupHists h;
  h.label = name;

  h.hTheta = new TH1D(
    Form("hTheta_%s_%s", tag.c_str(), name.c_str()),
    Form("%s;#theta_{lab} [deg];Weighted entries", name.c_str()),
    180, 0.0, 180.0
  );

  h.hCos = new TH1D(
    Form("hCos_%s_%s", tag.c_str(), name.c_str()),
    Form("%s;cos#theta_{lab};Weighted entries", name.c_str()),
    200, -1.0, 1.0
  );

  h.hPTheta = new TH2D(
    Form("hPTheta_%s_%s", tag.c_str(), name.c_str()),
    Form("%s;#theta_{lab} [deg];p_{lab} [GeV/c]", name.c_str()),
    180, 0.0, 180.0, 220, 0.0, 2.2
  );

  return h;
}

double MeanDedx_BetheBloch_MeV_per_cm(double p, double mass)
{
  // p, mass in GeV
  // return mean dE/dx in MeV/cm for plastic scintillator (polystyrene-like)

  const double K = 0.307075;        // MeV mol^-1 cm^2
  const double me = 0.00051099895;  // GeV
  const double z = 1.0;             // incident charge
  const double Z_over_A = 0.53768;  // polystyrene approx
  const double rho = 1.032;         // g/cm^3
  const double I = 68.7e-6;         // GeV (68.7 eV)
  const double delta = 0.0;         // neglect density correction in this regime

  if (p <= 0.0 || mass <= 0.0) return 0.0;

  const double E = std::sqrt(p*p + mass*mass);
  const double beta = p / E;
  const double gamma = E / mass;

  if (beta <= 0.0 || beta >= 1.0) return 0.0;

  const double bg = beta * gamma;

  // Tmax in GeV
  const double ratio = me / mass;
  const double Tmax =
    (2.0 * me * bg * bg) /
    (1.0 + 2.0 * gamma * ratio + ratio * ratio);

  if (Tmax <= 0.0) return 0.0;

  const double arg =
    (2.0 * me * bg * bg * Tmax) / (I * I);

  if (arg <= 1.0) return 0.0;

  // mass stopping power [MeV cm^2 / g]
  const double dedx_mass =
    K * z * z * Z_over_A * (1.0 / (beta * beta)) *
    (0.5 * std::log(arg) - beta * beta - 0.5 * delta);

  // convert to MeV/cm
  return dedx_mass * rho;
}

double MeanDeltaEInScint_MeV(double p,
                             double mass,
                             double thetaDeg,
                             double thickness_cm = 0.5,
                             bool useAngularPath = true)
{
  // thickness_cm = scintillator thickness in cm
  // default 0.5 cm = 5 mm

  const double dedx = MeanDedx_BetheBloch_MeV_per_cm(p, mass);
  if (dedx <= 0.0) return 0.0;

  double path = thickness_cm;

  if (useAngularPath) {
    const double th = thetaDeg * TMath::DegToRad();
    const double c = std::fabs(std::cos(th));

    // avoid divergence at large theta
    if (c > 0.2) path = thickness_cm / c;
    else         path = thickness_cm / 0.2; // cap at 5x thickness
  }

  return dedx * path;
}

TH2D* ConvertPThetaToDeltaETheta(const TH2D* hIn,
                                 double mass,
                                 const char* newName,
                                 const char* newTitle,
                                 double thickness_cm = 0.5,
                                 bool useAngularPath = true,
                                 int nY = 250,
                                 double yMax = 15.0)
{
  TH2D* hOut = new TH2D(
    newName,
    newTitle,
    hIn->GetNbinsX(),
    hIn->GetXaxis()->GetXmin(),
    hIn->GetXaxis()->GetXmax(),
    nY, 0.0, yMax
  );

  for (int ix = 1; ix <= hIn->GetNbinsX(); ++ix) {
    const double theta = hIn->GetXaxis()->GetBinCenter(ix);

    for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy) {
      const double p = hIn->GetYaxis()->GetBinCenter(iy);
      const double w = hIn->GetBinContent(ix, iy);
      if (w <= 0.0) continue;

      const double dE = MeanDeltaEInScint_MeV(
        p, mass, theta, thickness_cm, useAngularPath
      );

      hOut->Fill(theta, dE, w);
    }
  }

  return hOut;
}

void FillParticle(ParticleHists& h, const TLorentzVector& p4, double w) {
  const double th = ThetaDeg(p4);
  const double c  = p4.Vect().CosTheta();
  const double p  = p4.P();

  h.hTheta->Fill(th, w);
  h.hCos->Fill(c, w);
  h.hPTheta->Fill(th, p, w);
}

void FillGroup(GroupHists& h, const TLorentzVector& p4, double w) {
  const double th = ThetaDeg(p4);
  const double c  = p4.Vect().CosTheta();
  const double p  = p4.P();

  h.hTheta->Fill(th, w);
  h.hCos->Fill(c, w);
  h.hPTheta->Fill(th, p, w);
}

ReactionHists BookReaction2() {
  ReactionHists r;
  r.tag = "R2";

  r.charged.push_back(BookParticle(r.tag, "p"));
  r.charged.push_back(BookParticle(r.tag, "piMinus_prod"));
  r.charged.push_back(BookParticle(r.tag, "piPlus"));
  r.charged.push_back(BookParticle(r.tag, "piMinus_L1405"));

  r.protonGroup = BookGroup(r.tag, "proton_all");
  r.pionGroup   = BookGroup(r.tag, "pions_all");

  r.hMParent = new TH1D(
    "hMParent_R2",
    "Sampled M(#Sigma(1670));M [GeV/c^{2}];Weighted events",
    240, 1.55, 1.76
  );

  r.hMSub = new TH1D(
    "hMSub_R2",
    "Sampled M(#Lambda(1405));M [GeV/c^{2}];Weighted events",
    220, 1.30, 1.50
  );

  return r;
}

bool SampleMassBW(TRandom3& rng,
                  const Resonance& res,
                  double mMin,
                  double mMax,
                  double& mOut,
                  int maxTry = 200000)
{
  if (mMax <= mMin) return false;

  for (int i = 0; i < maxTry; ++i) {
    double m = rng.BreitWigner(res.m0, res.gamma);
    if (m >= mMin && m <= mMax) {
      mOut = m;
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

bool SimulateReaction2(const TLorentzVector& initial,
                       TRandom3& rng,
                       ReactionHists& h,
                       double& totalWeight)
{
  // K- p -> Sigma(1670)+ + pi-_prod
  // Sigma(1670)+ -> Lambda(1405) + pi+ + pi0
  // Lambda(1405) -> Sigma+ + pi-
  // Sigma+ -> p + pi0

  TGenPhaseSpace prod, dec1, dec2, dec3;

  const double mPi  = Masses::Pion;
  const double mPi0 = Masses::Pi0;
  const double mSp  = Masses::SigmaP;
  const double mp   = Masses::Proton;

  const double sqrtS = initial.M();

  // Sigma(1670) allowed range
  const double mSigmaMin = mSp + mPi + mPi + mPi0;
  const double mSigmaMax = sqrtS - mPi;

  double mSigma = -1.0;
  if (!SampleMassBW(rng, Res::Sigma1670, mSigmaMin, mSigmaMax, mSigma)) return false;

  double wProd = 0.0;
  if (!GenerateDecay(prod, initial, {mSigma, mPi}, wProd)) return false;

  TLorentzVector* pSigmaStar   = prod.GetDecay(0);
  TLorentzVector* pPiMinusProd = prod.GetDecay(1);

  // Lambda(1405) allowed range
  const double mLMin = mSp + mPi;
  const double mLMax = mSigma - mPi - mPi0;
  if (mLMax <= mLMin) return false;

  double mL1405 = -1.0;
  if (!SampleMassBW(rng, Res::Lambda1405, mLMin, mLMax, mL1405)) return false;

  double w1 = 0.0;
  if (!GenerateDecay(dec1, *pSigmaStar, {mL1405, mPi, mPi0}, w1)) return false;

  TLorentzVector* pL1405 = dec1.GetDecay(0);
  TLorentzVector* pPiPlus = dec1.GetDecay(1);

  double w2 = 0.0;
  if (!GenerateDecay(dec2, *pL1405, {mSp, mPi}, w2)) return false;

  TLorentzVector* pSigmaP       = dec2.GetDecay(0);
  TLorentzVector* pPiMinusL1405 = dec2.GetDecay(1);

  double w3 = 0.0;
  if (!GenerateDecay(dec3, *pSigmaP, {mp, mPi0}, w3)) return false;

  TLorentzVector* pProton = dec3.GetDecay(0);

  totalWeight = wProd * w1 * w2 * w3;

  // individual
  FillParticle(h.charged[0], *pProton,       totalWeight);
  FillParticle(h.charged[1], *pPiMinusProd,  totalWeight);
  FillParticle(h.charged[2], *pPiPlus,       totalWeight);
  FillParticle(h.charged[3], *pPiMinusL1405, totalWeight);

  // grouped
  FillGroup(h.protonGroup, *pProton, totalWeight);
  FillGroup(h.pionGroup,   *pPiMinusProd,  totalWeight);
  FillGroup(h.pionGroup,   *pPiPlus,       totalWeight);
  FillGroup(h.pionGroup,   *pPiMinusL1405, totalWeight);

  h.hMParent->Fill(mSigma, totalWeight);
  h.hMSub->Fill(mL1405, totalWeight);

  return true;
}
void Save2DOverlapPdf(const ReactionHists& r, const string& pdfName)
{
  TCanvas* c = new TCanvas("c2DOverlap", "2D overlap", 900, 700);
  c->SetMargin(0.12, 0.12, 0.12, 0.08);

  // clone so original histograms are untouched
  TH2D* hP  = (TH2D*)r.protonGroup.hPTheta->Clone("hP_overlap");
  TH2D* hPi = (TH2D*)r.pionGroup.hPTheta->Clone("hPi_overlap");

  // normalize only if you want shape comparison
  if (hP->Integral()  > 0) hP->Scale(1.0 / hP->Integral());
  if (hPi->Integral() > 0) hPi->Scale(1.0 / hPi->Integral());

  hP->SetTitle("proton vs all charged pions;#theta_{lab} [deg];p_{lab} [GeV/c]");
  hP->SetStats(0);
  hPi->SetStats(0);

  // contour levels
  const int nCont = 6;
  double pLevels[nCont]  = {0.05, 0.10, 0.20, 0.35, 0.55, 0.75};
  double piLevels[nCont] = {0.05, 0.10, 0.20, 0.35, 0.55, 0.75};

  // scale contour levels by histogram maximum
  double pMax  = hP->GetMaximum();
  double piMax = hPi->GetMaximum();

  for (int i = 0; i < nCont; ++i) {
    pLevels[i]  *= pMax;
    piLevels[i] *= piMax;
  }

  hP->SetContour(nCont, pLevels);
  hPi->SetContour(nCont, piLevels);

  hP->SetLineColor(kRed+1);
  hP->SetLineWidth(3);

  hPi->SetLineColor(kBlue+1);
  hPi->SetLineWidth(3);

  // draw empty frame first
  TH2D* frame = new TH2D("frame_overlap",
                         "proton vs all charged pions;#theta_{lab} [deg];p_{lab} [GeV/c]",
                         180, 0.0, 180.0, 220, 0.0, 2.2);
  frame->SetStats(0);
  frame->Draw();

  hP->Draw("cont3 same");
  hPi->Draw("cont3 same");

  TLegend* leg = new TLegend(0.58, 0.76, 0.88, 0.88);
  leg->AddEntry(hP,  "proton", "l");
  leg->AddEntry(hPi, "#pi^{#pm}", "l");
  leg->Draw();

  c->SaveAs(pdfName.c_str());

  //delete leg;
  //delete frame;
  //delete c;
}

void SaveDeltaEOverlapPdf(const ReactionHists& r,
                          const string& pdfName,
                          double thickness_cm = 0.5,
                          bool useAngularPath = true)
{
  const double mp  = Masses::Proton;
  const double mpi = Masses::Pion;

  TH2D* hP_dE = ConvertPThetaToDeltaETheta(
    r.protonGroup.hPTheta,
    mp,
    "hP_dE",
    "proton vs all charged pions;#theta_{lab} [deg];#DeltaE in plastic scintillator [MeV]",
    thickness_cm,
    useAngularPath,
    250, 15.0
  );

  TH2D* hPi_dE = ConvertPThetaToDeltaETheta(
    r.pionGroup.hPTheta,
    mpi,
    "hPi_dE",
    "proton vs all charged pions;#theta_{lab} [deg];#DeltaE in plastic scintillator [MeV]",
    thickness_cm,
    useAngularPath,
    250, 15.0
  );

  if (hP_dE->Integral()  > 0) hP_dE->Scale(1.0 / hP_dE->Integral());
  if (hPi_dE->Integral() > 0) hPi_dE->Scale(1.0 / hPi_dE->Integral());

  TCanvas* c = new TCanvas("cDeltaEOverlap", "DeltaE overlap", 950, 750);
  c->SetMargin(0.14, 0.12, 0.12, 0.08);

  TH2D* frame = new TH2D(
    "frame_dE",
    Form("proton vs all charged pions in %g mm plastic scintillator;#theta_{lab} [deg];#DeltaE [MeV]",
         thickness_cm * 10.0),
    180, 0.0, 180.0, 250, 0.0, 15.0
  );
  frame->SetStats(0);
  frame->Draw();

  const int nCont = 6;
  double pLevels[nCont]  = {0.05, 0.10, 0.20, 0.35, 0.55, 0.75};
  double piLevels[nCont] = {0.05, 0.10, 0.20, 0.35, 0.55, 0.75};

  const double pMax  = hP_dE->GetMaximum();
  const double piMax = hPi_dE->GetMaximum();

  for (int i = 0; i < nCont; ++i) {
    pLevels[i]  *= pMax;
    piLevels[i] *= piMax;
  }

  hP_dE->SetContour(nCont, pLevels);
  hPi_dE->SetContour(nCont, piLevels);

  hP_dE->SetLineColor(kRed+1);
  hPi_dE->SetLineColor(kBlue+1);
  hP_dE->SetLineWidth(3);
  hPi_dE->SetLineWidth(3);

  hP_dE->Draw("cont3 same");
  hPi_dE->Draw("cont3 same");

  TLegend* leg = new TLegend(0.62, 0.74, 0.88, 0.88);
  leg->AddEntry(hP_dE,  "proton", "l");
  leg->AddEntry(hPi_dE, "all charged pions", "l");
  leg->Draw();

  c->SaveAs(pdfName.c_str());

  delete leg;
  delete frame;
  delete c;
}

void WriteReaction(const ReactionHists& r) {
  if (r.hMParent) r.hMParent->Write();
  if (r.hMSub)    r.hMSub->Write();

  for (const auto& p : r.charged) {
    p.hTheta->Write();
    p.hCos->Write();
    p.hPTheta->Write();
  }

  r.protonGroup.hTheta->Write();
  r.protonGroup.hCos->Write();
  r.protonGroup.hPTheta->Write();

  r.pionGroup.hTheta->Write();
  r.pionGroup.hCos->Write();
  r.pionGroup.hPTheta->Write();
}

void SaveIndividualPdf(const ReactionHists& r, const string& pdfName) {
  TCanvas* c = new TCanvas("cIndividual", "individual", 1400, 900);

  c->Divide(2,2);
  for (int i = 0; i < 4; ++i) {
    c->cd(i+1);
    r.charged[i].hTheta->SetLineWidth(2);
    r.charged[i].hTheta->Draw("hist");
  }
  c->SaveAs((pdfName + "(").c_str());

  c->Clear();
  c->Divide(2,2);
  for (int i = 0; i < 4; ++i) {
    c->cd(i+1);
    r.charged[i].hCos->SetLineWidth(2);
    r.charged[i].hCos->Draw("hist");
  }
  c->SaveAs(pdfName.c_str());

  c->Clear();
  c->Divide(2,2);
  for (int i = 0; i < 4; ++i) {
    c->cd(i+1);
    r.charged[i].hPTheta->Draw("colz");
  }
  c->SaveAs((pdfName + ")").c_str());

  delete c;
}

void SaveGroupedPdf(const ReactionHists& r, const string& pdfName) {
  TCanvas* c = new TCanvas("cGrouped", "grouped", 1200, 900);

  c->Divide(2,2);

  c->cd(1);
  r.protonGroup.hTheta->SetLineWidth(2);
  r.protonGroup.hTheta->Draw("hist");

  c->cd(2);
  r.pionGroup.hTheta->SetLineWidth(2);
  r.pionGroup.hTheta->Draw("hist");

  c->cd(3);
  r.protonGroup.hPTheta->Draw("colz");

  c->cd(4);
  r.pionGroup.hPTheta->Draw("colz");

  c->SaveAs(pdfName.c_str());

  delete c;
}

void SaveOverlayPdf(const ReactionHists& r, const string& pdfName) {
  TCanvas* c = new TCanvas("cOverlay", "overlay", 1200, 900);
  c->Divide(2,2);

  // theta overlay
  c->cd(1);
  TH1D* hThetaP = (TH1D*)r.protonGroup.hTheta->Clone("hThetaP_tmp");
  TH1D* hThetaPi = (TH1D*)r.pionGroup.hTheta->Clone("hThetaPi_tmp");
  hThetaP->Scale(1.0 / hThetaP->Integral());
  hThetaPi->Scale(1.0 / hThetaPi->Integral());
  hThetaP->SetLineColor(kRed+1);
  hThetaPi->SetLineColor(kBlue+1);
  hThetaP->SetLineWidth(2);
  hThetaPi->SetLineWidth(2);
  hThetaP->SetTitle("Normalized #theta_{lab};#theta_{lab} [deg];Normalized entries");
  hThetaP->Draw("hist");
  hThetaPi->Draw("hist same");
  {
    TLegend* leg = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg->AddEntry(hThetaP, "proton", "l");
    leg->AddEntry(hThetaPi, "all charged pions", "l");
    leg->Draw();
  }

  // cos overlay
  c->cd(2);
  TH1D* hCosP = (TH1D*)r.protonGroup.hCos->Clone("hCosP_tmp");
  TH1D* hCosPi = (TH1D*)r.pionGroup.hCos->Clone("hCosPi_tmp");
  hCosP->Scale(1.0 / hCosP->Integral());
  hCosPi->Scale(1.0 / hCosPi->Integral());
  hCosP->SetLineColor(kRed+1);
  hCosPi->SetLineColor(kBlue+1);
  hCosP->SetLineWidth(2);
  hCosPi->SetLineWidth(2);
  hCosP->SetTitle("Normalized cos#theta_{lab};cos#theta_{lab};Normalized entries");
  hCosP->Draw("hist");
  hCosPi->Draw("hist same");
  {
    TLegend* leg = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg->AddEntry(hCosP, "proton", "l");
    leg->AddEntry(hCosPi, "all charged pions", "l");
    leg->Draw();
  }

  // p projection overlay
  c->cd(3);
  TH1D* hPP = r.protonGroup.hPTheta->ProjectionY("hPP_tmp");
  TH1D* hPPi = r.pionGroup.hPTheta->ProjectionY("hPPi_tmp");
  hPP->Scale(1.0 / hPP->Integral());
  hPPi->Scale(1.0 / hPPi->Integral());
  hPP->SetLineColor(kRed+1);
  hPPi->SetLineColor(kBlue+1);
  hPP->SetLineWidth(2);
  hPPi->SetLineWidth(2);
  hPP->SetTitle("Normalized p_{lab};p_{lab} [GeV/c];Normalized entries");
  hPP->Draw("hist");
  hPPi->Draw("hist same");
  {
    TLegend* leg = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg->AddEntry(hPP, "proton", "l");
    leg->AddEntry(hPPi, "all charged pions", "l");
    leg->Draw();
  }

  // 2D just proton
  c->cd(4);
  r.protonGroup.hPTheta->SetTitle("proton p vs #theta");
  r.protonGroup.hPTheta->Draw("colz");

  c->SaveAs(pdfName.c_str());

  //delete c;
}

void SaveMassPdf(const ReactionHists& r, const string& pdfName) {
  TCanvas* c = new TCanvas("cMass", "mass", 1000, 450);
  c->Divide(2,1);

  c->cd(1);
  r.hMParent->SetLineWidth(2);
  r.hMParent->Draw("hist");

  c->cd(2);
  r.hMSub->SetLineWidth(2);
  r.hMSub->Draw("hist");

  c->SaveAs(pdfName.c_str());
  delete c;
}

void sigma1670_bw_lab_grouped(Long64_t nEvent = 500000, double pBeam = 1.17)
{
  gStyle->SetOptStat(0);

  const double mK = Masses::Kminus;
  const double mp = Masses::Proton;

  const double eBeam = std::sqrt(pBeam*pBeam + mK*mK);

  TLorentzVector beam(0.0, 0.0, pBeam, eBeam);
  TLorentzVector target(0.0, 0.0, 0.0, mp);
  TLorentzVector initial = beam + target;

  cout << "=====================================\n";
  cout << " Beam momentum = " << pBeam << " GeV/c\n";
  cout << " sqrt(s)       = " << initial.M() << " GeV\n";
  cout << "=====================================\n";

  ReactionHists r2 = BookReaction2();
  TRandom3 rng(0);

  Long64_t acc = 0;
  double sumW = 0.0;

  for (Long64_t i = 0; i < nEvent; ++i) {
    double w = 0.0;
    if (SimulateReaction2(initial, rng, r2, w)) {
      acc++;
      sumW += w;
    }
  }

  cout << "\nAccepted events = " << acc << " / " << nEvent << endl;
  cout << "Sum of weights  = " << sumW << endl;

  TFile* fout = new TFile("sigma1670_bw_lab_grouped.root", "RECREATE");
  WriteReaction(r2);
  fout->Close();

  SaveIndividualPdf(r2, "reaction2_individual.pdf");
  SaveGroupedPdf(r2, "reaction2_grouped.pdf");
  SaveOverlayPdf(r2, "reaction2_overlay.pdf");
  SaveMassPdf(r2, "reaction2_masses.pdf");
  
  Save2DOverlapPdf(r2, "reaction2_2D_overlap.pdf");
  SaveDeltaEOverlapPdf(r2, "reaction2_deltaE_overlap_5mm.pdf", 1.0, false);

  cout << "\nOutput files:\n";
  cout << "  sigma1670_bw_lab_grouped.root\n";
  cout << "  reaction2_individual.pdf\n";
  cout << "  reaction2_grouped.pdf\n";
  cout << "  reaction2_overlay.pdf\n";
  cout << "  reaction2_masses.pdf\n";
}
