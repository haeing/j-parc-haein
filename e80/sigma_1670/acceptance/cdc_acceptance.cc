// simple_cdc_mc.cc
// g++ -O2 -std=c++17 simple_cdc_mc.cc $(root-config --cflags --libs) -o simple_cdc_mc
//
// 목적:
//  K- (1.17 GeV/c) + p(at rest) -> Sigma(1670) + pi  (2-body)
//  Sigma(1670) -> K0 + p + pi-  (3-body)
//  K0S -> pi+ + pi-             (2-body, with decay length)
//  Uniform Bz=0.7 T에서 charged 트랙을 헬릭스로 전파하여 CDC 진입 여부(acceptance)를 계산

#include <iostream>
#include <vector>
#include <cmath>

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

static constexpr double kPi = 3.14159265358979323846;

// ------------------------ Geometry (mm) ------------------------
struct TargetCyl {
  double R;      // mm
  double halfZ;  // mm
};

struct CDCTube {
  double Rin;    // mm
  double Rout;   // mm
  double halfZ;  // mm
};

// ------------------------ Helix propagation ------------------------
// B field along +z. Charge q = +/-1 (in units of e).
// Momentum in GeV/c, positions in mm.
// r_L [mm] = (pT [GeV/c] / (0.3 * B[T])) * 1000
//
// "CDC에 들어갔는지"는 step으로 판단 (Geant4처럼 정밀 교차 계산 대신 간단/안정).
struct HitResult {
  bool entered = false;
  double entry_s_mm = 0.0;
  double entry_r_mm = 0.0;
  double entry_z_mm = 0.0;
};

HitResult propagateToCDC_helixBz(
    const TVector3& x0_mm,
    const TVector3& p_GeV,
    int q,
    double Bz_T,
    const CDCTube& cdc,
    double step_mm = 1.0,
    double max_s_mm = 5000.0 // 충분히 크게
) {
  HitResult res;

  const double px = p_GeV.X();
  const double py = p_GeV.Y();
  const double pz = p_GeV.Z();
  const double p  = p_GeV.Mag();

  //if (p <= 0) return res;

  const double pT = std::sqrt(px*px + py*py);

  // Neutral이면 직선 전파로 처리(여기서는 charged만 호출한다고 가정해도 됨)
  if (q == 0 || std::abs(Bz_T) < 1e-12 || pT < 1e-12) {
    // straight line: x(s) = x0 + (p_hat)*s
    TVector3 phat = p_GeV.Unit();
    for (double s = 0; s <= max_s_mm; s += step_mm) {
      TVector3 x = x0_mm + phat * s;
      double r = std::sqrt(x.X()*x.X() + x.Y()*x.Y());
      if (!res.entered &&
          r >= cdc.Rin && r <= cdc.Rout &&
          std::abs(x.Z()) <= cdc.halfZ) {
        res.entered = true;
        res.entry_s_mm = s;
        res.entry_r_mm = r;
        res.entry_z_mm = x.Z();
        return res;
      }
    }
    return res;
  }

  // curvature radius in transverse plane (mm)
  const double rL_mm = (pT / (0.3 * std::abs(Bz_T))) * 1000.0;

  // helix angular frequency in terms of path length:
  // dphi/ds = (q * 0.3 * B / p) [rad/m] => [rad/mm] = ... /1000
  // but for helix stepping, we can use: dphi = (q * step / rL) * (p/pT)
  // A more stable way:
  // In xy plane, curvature kappa = 1/rL. For arc length ds, dphi_xy = ds / rL (sign with q*B).
  const double sign = (q > 0 ? +1.0 : -1.0) * (Bz_T >= 0 ? +1.0 : -1.0);

  // initial transverse angle of momentum
  const double phi0 = std::atan2(py, px);

  // Circle center in xy:
  // For a charged particle in +Bz, the center is shifted by rL rotated by +/-90deg from momentum.
  // x_c = x0 + sign * rL * (-sin(phi0), +cos(phi0))
  const double xc = x0_mm.X() + sign * rL_mm * (-std::sin(phi0));
  const double yc = x0_mm.Y() + sign * rL_mm * ( std::cos(phi0));

  // Starting point on the circle has angle theta0 around center:
  // x0 = xc + rL*cos(theta0), y0 = yc + rL*sin(theta0)
  const double theta0 = std::atan2(x0_mm.Y() - yc, x0_mm.X() - xc);

  // pitch: dz/ds = pz/p
  const double dzds = pz / p;

  for (double s = 0; s <= max_s_mm; s += step_mm) {
    // transverse rotation angle along arc length
    const double dtheta = sign * (s / rL_mm);  // rad
    const double theta = theta0 + dtheta;

    const double x = xc + rL_mm * std::cos(theta);
    const double y = yc + rL_mm * std::sin(theta);
    const double z = x0_mm.Z() + dzds * s;

    const double r = std::sqrt(x*x + y*y);

    if (!res.entered &&
        r >= cdc.Rin && r <= cdc.Rout &&
        std::abs(z) <= cdc.halfZ) {
      res.entered = true;
      res.entry_s_mm = s;
      res.entry_r_mm = r;
      res.entry_z_mm = z;
      return res;
    }
  }

  return res;
}

// ------------------------ Helpers ------------------------
TVector3 sampleVertexInTarget(TRandom3& rng, const TargetCyl& tgt) {
  // Uniform in cylinder volume:
  // r = R*sqrt(u), phi=2pi*v, z uniform in [-halfZ, halfZ]
  const double u = rng.Uniform();
  const double v = rng.Uniform();
  const double r = tgt.R * std::sqrt(u);
  const double phi = 2.0*kPi*v;
  const double x = r * std::cos(phi);
  const double y = r * std::sin(phi);
  const double z = rng.Uniform(-tgt.halfZ, +tgt.halfZ);
  return TVector3(x, y, z);
}

// Exponential decay length (mm) with mean = beta*gamma*c*tau
double sampleDecayLength_mm(TRandom3& rng, double betaGamma, double ctau_mm) {
  const double u = rng.Uniform();
  return -betaGamma * ctau_mm * std::log(u);
}

// ------------------------ Main ------------------------
void cdc_acceptance(){
  long long N = 200000;

  // Geometry
  TargetCyl tgt{30.0, 75.0};            // R=30 mm, halfZ=75 mm
  CDCTube   cdc{150.0, 530.0, 1340.0};  // Rin=150, Rout=530, halfZ=1340 (mm)

  // B field
  const double Bz = 0.7; // Tesla

  // PDG masses (GeV)
  auto* pdg = TDatabasePDG::Instance();
  const double mKminus = pdg->GetParticle(-321)->Mass();
  const double mp      = pdg->GetParticle(2212)->Mass();
  const double mpi     = pdg->GetParticle(211)->Mass();   // pi+ mass; pi- same
  const double mK0     = pdg->GetParticle(311)->Mass();   // K0
  // Sigma(1670) mass: PDG entry may vary; we set a fixed nominal (GeV)
  const double mSigma1670 = 1.670;

  // K0S ctau (mm): c*tau ~ 2.684 cm = 26.84 mm
  // (K0S lifetime 0.895e-10 s -> 2.684 cm)
  const double ctauK0S_mm = 26.84;

  // Beam
  const double pBeam = 1.17; // GeV/c
  TLorentzVector Kbeam(0, 0, pBeam, std::sqrt(pBeam*pBeam + mKminus*mKminus));
  TLorentzVector pTarget(0, 0, 0, mp);

  TLorentzVector Pinit = Kbeam + pTarget;

  TRandom3 rng(0);

  long long nAll = 0;
  long long nAcc = 0;

  // phase space generators
  TGenPhaseSpace gen2; // K- p -> Sigma + pi
  TGenPhaseSpace gen3; // Sigma -> K0 + p + pi-
  TGenPhaseSpace genK0; // K0S -> pi+ + pi-
  
  double masses2[2] = {mSigma1670, mpi};
  double masses3[3] = {mK0, mp, mpi};      // (K0, p, pi-)
  double massesK0[2] = {mpi, mpi};         // (pi+, pi-)

  const int nb = 20;
  auto hAll = new TH1D("hAll",";cos#theta*;Events", nb, -1.0, 1.0);
  auto hAcc = new TH1D("hAcc",";cos#theta*;Accepted", nb, -1.0, 1.0);
  
  auto hCosLab_allCh = new TH1D("hCosLab_allCh", ";cos#theta_{lab};Tracks", nb, -1.0, 1.0);
  auto hCosLab_accCh = new TH1D("hCosLab_accCh", ";cos#theta_{lab};Tracks (accepted events)", nb, -1.0, 1.0);

  for (long long i = 0; i < N; i++) {
    nAll++;

    // 1) vertex in target
    TVector3 vtx = sampleVertexInTarget(rng, tgt);

    // 2) primary 2-body: K- p -> Sigma(1670) + pi
    if (!gen2.SetDecay(Pinit, 2, masses2)) continue;
    gen2.Generate();
    

    TLorentzVector* pSigma = gen2.GetDecay(0);
    TLorentzVector* pPi1   = gen2.GetDecay(1); // associated pion (charge depends on channel; here just treat as charged)


    TVector3 betaCM = Pinit.BoostVector();
    TLorentzVector p_p_cm   = pTarget;   p_p_cm.Boost(-betaCM);
    TLorentzVector p_sig_cm = *pSigma;   p_sig_cm.Boost(-betaCM);

    // cos(theta*) between initial target proton direction and Sigma direction in CM
    double cosTheta = p_p_cm.Vect().Unit().Dot(p_sig_cm.Vect().Unit());

    hAll->Fill(cosTheta);
    // 3) Sigma 3-body: Sigma -> K0 + p + pi-
    if (!gen3.SetDecay(*pSigma, 3, masses3)) continue;
    gen3.Generate();

    TLorentzVector* pK0    = gen3.GetDecay(0);
    TLorentzVector* pProt  = gen3.GetDecay(1);
    TLorentzVector* pPim   = gen3.GetDecay(2);

    // 4) K0 flight + decay to pi+ pi- (assume K0S)
    // decay point = vtx + direction*K0_decay_length
    TVector3 dirK0 = pK0->Vect().Unit();
    const double betaGammaK0 = pK0->P() / pK0->M(); // = |p|/m
    const double Ldec_mm = sampleDecayLength_mm(rng, betaGammaK0, ctauK0S_mm);
    TVector3 vtxK0 = vtx + dirK0 * Ldec_mm;

    if (!genK0.SetDecay(*pK0, 2, massesK0)) continue;
    genK0.Generate();

    TLorentzVector* pPip = genK0.GetDecay(0);
    TLorentzVector* pPim2= genK0.GetDecay(1);

    // 5) CDC entry check for charged tracks
    // charges: pi1 (use -1 as an example), proton (+1), pi- (-1), pi+ (+1), pi- (-1)
    // NOTE: 실제 채널에 맞게 pi1 charge는 바꿔줘.
    struct Track { TVector3 x0; TVector3 p; int q; };
    std::vector<Track> tracks;
    tracks.push_back({vtx,   pPi1->Vect(),  -1});
    tracks.push_back({vtx,   pProt->Vect(), +1});
    tracks.push_back({vtx,   pPim->Vect(),  -1});
    tracks.push_back({vtxK0, pPip->Vect(),  +1});
    tracks.push_back({vtxK0, pPim2->Vect(), -1});

    bool allEntered = true;
    for (const auto& tr : tracks) {
      auto hit = propagateToCDC_helixBz(tr.x0, tr.p, tr.q, Bz, cdc, 1.0, 8000.0);
      if (!hit.entered) { allEntered = false; break; }
      const double p = tr.p.Mag();
      if(p <=0) continue;
      const double cosLab = tr.p.Z() / p;
      hCosLab_allCh->Fill(cosLab);
		  
    }

    if (allEntered) {
      nAcc++;
      hAcc->Fill(cosTheta);
      for(const auto& tr : tracks){
	const double p = tr.p.Mag();
	if(p<=0)continue;
	double cosLab = tr.p.Z() / p;
	hCosLab_accCh->Fill(cosLab);
      }
    }
  }

  const double acc = (nAll > 0) ? (double)nAcc / (double)nAll : 0.0;

  std::cout << "N total = " << nAll << "\n";
  std::cout << "N accepted(all charged enter CDC) = " << nAcc << "\n";
  std::cout << "Acceptance = " << acc << "\n";

  auto c1 = new TCanvas("c1","c1");
  c1->Divide(3);
  c1->cd(1);
  hAll->Draw();
  c1->cd(2);
  hAcc->Draw();
  auto gAcc = new TGraphAsymmErrors();
  gAcc->BayesDivide(hAcc,hAll);
  c1->cd(3);
  gAcc->SetMarkerStyle(20);
  gAcc->Draw();

  auto c2 = new TCanvas("c2","c2");
  c2->Divide(2);
  c2->cd(1);
  hCosLab_allCh->Draw();
  c2->cd(2);
  hCosLab_accCh->Draw();
  
  
}
