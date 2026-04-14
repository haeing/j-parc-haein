#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "TRandom3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGenPhaseSpace.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

struct Resonance {
  string name;
  double m0;     // pole mass [GeV]
  double gamma;  // width [GeV]
};

struct ParticleHists {
  string label;
  TH1D* hTheta = nullptr;      // theta_lab [deg]
  TH1D* hCos   = nullptr;      // cos(theta_lab)
  TH2D* hPTheta = nullptr;     // p_lab vs theta_lab
};

struct ReactionHists {
  string tag;
  vector<ParticleHists> charged;
  TH1D* hMParent = nullptr;
  TH1D* hMSub    = nullptr; // optional
};

namespace Masses {
  const double Kminus   = 0.493677;
  const double Proton   = 0.9382720813;
  const double Pion     = 0.13957039;
  const double Pi0      = 0.1349768;
  const double Gamma    = 0.0;
  const double Lambda   = 1.115683;
  const double Sigma0   = 1.192642;
  const double Sigmap   = 1.18937;
}

namespace Res {
  const Resonance Sigma1670 = {"Sigma1670", 1.670, 0.060};
  const Resonance Lambda1405 = {"Lambda1405", 1.405, 0.050};
}

double ThetaDeg(const TLorentzVector& p4) {
  return p4.Vect().Theta() * TMath::RadToDeg();
}

void FillParticle(ParticleHists& h, const TLorentzVector& p4, double w) {
  const double theta = ThetaDeg(p4);
  const double costh = p4.Vect().CosTheta();
  const double p     = p4.P();

  h.hTheta->Fill(theta, w);
  h.hCos->Fill(costh, w);
  h.hPTheta->Fill(theta, p, w);
}

ParticleHists BookParticle(const string& tag, const string& name) {
  ParticleHists h;
  h.label = name;

  h.hTheta = new TH1D(
    Form("hTheta_%s_%s", tag.c_str(), name.c_str()),
    Form("%s ;#theta_{lab} [deg];Weighted events", name.c_str()),
    180, 0.0, 180.0
  );

  h.hCos = new TH1D(
    Form("hCos_%s_%s", tag.c_str(), name.c_str()),
    Form("%s ;cos#theta_{lab};Weighted events", name.c_str()),
    200, -1.0, 1.0
  );

  h.hPTheta = new TH2D(
    Form("hPTheta_%s_%s", tag.c_str(), name.c_str()),
    Form("%s ;#theta_{lab} [deg];p_{lab} [GeV/c]", name.c_str()),
    180, 0.0, 180.0, 220, 0.0, 2.2
  );

  return h;
}

ReactionHists BookReaction1() {
  ReactionHists r;
  r.tag = "R1";

  r.charged.push_back(BookParticle(r.tag, "p"));
  r.charged.push_back(BookParticle(r.tag, "piMinus_prod"));
  r.charged.push_back(BookParticle(r.tag, "piPlus"));
  r.charged.push_back(BookParticle(r.tag, "piMinus_Lambda"));

  r.hMParent = new TH1D(
    "hMParent_R1",
    "Reaction1 sampled M(#Sigma(1670));M [GeV/c^{2}];Weighted events",
    240, 1.40, 1.76
  );

  r.hMSub = nullptr;
  return r;
}

ReactionHists BookReaction2() {
  ReactionHists r;
  r.tag = "R2";

  r.charged.push_back(BookParticle(r.tag, "p"));
  r.charged.push_back(BookParticle(r.tag, "piMinus_prod"));
  r.charged.push_back(BookParticle(r.tag, "piPlus"));
  r.charged.push_back(BookParticle(r.tag, "piMinus_L1405"));

  r.hMParent = new TH1D(
    "hMParent_R2",
    "Reaction2 sampled M(#Sigma(1670));M [GeV/c^{2}];Weighted events",
    240, 1.55, 1.76
  );

  r.hMSub = new TH1D(
    "hMSub_R2",
    "Reaction2 sampled M(#Lambda(1405));M [GeV/c^{2}];Weighted events",
    220, 1.30, 1.50
  );

  return r;
}

// truncated Breit-Wigner sampling
bool SampleMassBW(TRandom3& rng,
                  const Resonance& res,
                  double mMin,
                  double mMax,
                  double& mOut,
                  int maxTry = 200000)
{
  if (mMax <= mMin) return false;

  for (int i = 0; i < maxTry; ++i) {
    const double m = rng.BreitWigner(res.m0, res.gamma);
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
  TLorentzVector motherCopy = mother;   // non-const copy for ROOT
  if (!event.SetDecay(motherCopy, (int)m.size(), m.data())) return false;
  weight = event.Generate();
  return true;
}

void WriteReaction(const ReactionHists& r) {
  if (r.hMParent) r.hMParent->Write();
  if (r.hMSub)    r.hMSub->Write();

  for (const auto& p : r.charged) {
    p.hTheta->Write();
    p.hCos->Write();
    p.hPTheta->Write();
  }
}

void SaveReactionPdf(const ReactionHists& r, const string& pdfName) {
  TCanvas* c = new TCanvas(Form("c_%s", r.tag.c_str()), r.tag.c_str(), 1400, 900);

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

bool SimulateReaction1(const TLorentzVector& initial,
                       TRandom3& rng,
                       ReactionHists& h,
                       double& totalWeight)
{
  // K- p -> Sigma(1670)+ + pi-_prod
  // Sigma(1670)+ -> Sigma0 + pi+ + pi0
  // Sigma0 -> Lambda + gamma
  // Lambda -> p + pi-

  TGenPhaseSpace prod, dec1, dec2, dec3;

  const double mPi   = Masses::Pion;
  const double mPi0  = Masses::Pi0;
  const double mS0   = Masses::Sigma0;
  const double mL    = Masses::Lambda;
  const double mp    = Masses::Proton;
  const double mGam  = Masses::Gamma;

  const double sqrtS = initial.M();

  // Sigma(1670) mass allowed by production and by decay threshold
  const double mSigmaMin = mS0 + mPi + mPi0;
  const double mSigmaMax = sqrtS - mPi;

  double mSigma = -1.0;
  if (!SampleMassBW(rng, Res::Sigma1670, mSigmaMin, mSigmaMax, mSigma)) return false;

  double wProd = 0.0;
  if (!GenerateDecay(prod, initial, {mSigma, mPi}, wProd)) return false;

  TLorentzVector* pSigmaStar   = prod.GetDecay(0);
  TLorentzVector* pPiMinusProd = prod.GetDecay(1);

  double w1 = 0.0;
  if (!GenerateDecay(dec1, *pSigmaStar, {mS0, mPi, mPi0}, w1)) return false;

  TLorentzVector* pSigma0 = dec1.GetDecay(0);
  TLorentzVector* pPiPlus = dec1.GetDecay(1);

  double w2 = 0.0;
  if (!GenerateDecay(dec2, *pSigma0, {mL, mGam}, w2)) return false;

  TLorentzVector* pLambda = dec2.GetDecay(0);

  double w3 = 0.0;
  if (!GenerateDecay(dec3, *pLambda, {mp, mPi}, w3)) return false;

  TLorentzVector* pProton        = dec3.GetDecay(0);
  TLorentzVector* pPiMinusLambda = dec3.GetDecay(1);

  totalWeight = wProd * w1 * w2 * w3;

  FillParticle(h.charged[0], *pProton,        totalWeight);
  FillParticle(h.charged[1], *pPiMinusProd,   totalWeight);
  FillParticle(h.charged[2], *pPiPlus,        totalWeight);
  FillParticle(h.charged[3], *pPiMinusLambda, totalWeight);

  h.hMParent->Fill(mSigma, totalWeight);

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

  const double mPi   = Masses::Pion;
  const double mPi0  = Masses::Pi0;
  const double mSp   = Masses::Sigmap;
  const double mp    = Masses::Proton;

  const double sqrtS = initial.M();

  // parent Sigma(1670) mass must allow full chain
  // minimal chain threshold = Sigma+ + pi- + pi+ + pi0
  const double mSigmaMin = mSp + mPi + mPi + mPi0;
  const double mSigmaMax = sqrtS - mPi;

  double mSigma = -1.0;
  if (!SampleMassBW(rng, Res::Sigma1670, mSigmaMin, mSigmaMax, mSigma)) return false;

  double wProd = 0.0;
  if (!GenerateDecay(prod, initial, {mSigma, mPi}, wProd)) return false;

  TLorentzVector* pSigmaStar   = prod.GetDecay(0);
  TLorentzVector* pPiMinusProd = prod.GetDecay(1);

  // Lambda(1405) mass must fit inside Sigma* -> L1405 + pi+ + pi0
  // and must also decay to Sigma+ + pi-
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

  TLorentzVector* pSigmap         = dec2.GetDecay(0);
  TLorentzVector* pPiMinusL1405   = dec2.GetDecay(1);

  double w3 = 0.0;
  if (!GenerateDecay(dec3, *pSigmap, {mp, mPi0}, w3)) return false;

  TLorentzVector* pProton = dec3.GetDecay(0);

  totalWeight = wProd * w1 * w2 * w3;

  FillParticle(h.charged[0], *pProton,       totalWeight);
  FillParticle(h.charged[1], *pPiMinusProd,  totalWeight);
  FillParticle(h.charged[2], *pPiPlus,       totalWeight);
  FillParticle(h.charged[3], *pPiMinusL1405, totalWeight);

  h.hMParent->Fill(mSigma, totalWeight);
  h.hMSub->Fill(mL1405, totalWeight);

  return true;
}

void sigma1670_bw_lab(Long64_t nEvent = 500000, double pBeam = 1.17)
{
  gStyle->SetOptStat(0);

  const double mK = Masses::Kminus;
  const double mp = Masses::Proton;

  const double eBeam = std::sqrt(pBeam*pBeam + mK*mK);

  TLorentzVector beam(0.0, 0.0, pBeam, eBeam);
  TLorentzVector target(0.0, 0.0, 0.0, mp);
  TLorentzVector initial = beam + target;

  const double sqrtS = initial.M();

  cout << "=====================================\n";
  cout << " Beam momentum = " << pBeam << " GeV/c\n";
  cout << " sqrt(s)       = " << sqrtS << " GeV\n";
  cout << "=====================================\n";

  ReactionHists r1 = BookReaction1();
  ReactionHists r2 = BookReaction2();

  TRandom3 rng(0);

  Long64_t acc1 = 0;
  Long64_t acc2 = 0;

  double sumW1 = 0.0;
  double sumW2 = 0.0;

  for (Long64_t i = 0; i < nEvent; ++i) {
    double w = 0.0;

    if (SimulateReaction1(initial, rng, r1, w)) {
      acc1++;
      sumW1 += w;
    }

    if (SimulateReaction2(initial, rng, r2, w)) {
      acc2++;
      sumW2 += w;
    }
  }

  cout << "\nAccepted events:\n";
  cout << "  Reaction 1 : " << acc1 << " / " << nEvent << endl;
  cout << "  Reaction 2 : " << acc2 << " / " << nEvent << endl;

  cout << "\nSum of weights:\n";
  cout << "  Reaction 1 : " << sumW1 << endl;
  cout << "  Reaction 2 : " << sumW2 << endl;

  TFile* fout = new TFile("sigma1670_bw_lab.root", "RECREATE");
  WriteReaction(r1);
  WriteReaction(r2);
  fout->Close();

  SaveReactionPdf(r1, "reaction1_bw_lab.pdf");
  SaveReactionPdf(r2, "reaction2_bw_lab.pdf");

  TCanvas* cMass = new TCanvas("cMass", "Mass distributions", 1200, 500);
  cMass->Divide(3,1);

  cMass->cd(1);
  r1.hMParent->SetLineWidth(2);
  r1.hMParent->Draw("hist");

  cMass->cd(2);
  r2.hMParent->SetLineWidth(2);
  r2.hMParent->Draw("hist");

  cMass->cd(3);
  r2.hMSub->SetLineWidth(2);
  r2.hMSub->Draw("hist");

  cMass->SaveAs("sampled_bw_masses.pdf");

  cout << "\nOutput files:\n";
  cout << "  sigma1670_bw_lab.root\n";
  cout << "  reaction1_bw_lab.pdf\n";
  cout << "  reaction2_bw_lab.pdf\n";
  cout << "  sampled_bw_masses.pdf\n";
  cout << endl;
}
