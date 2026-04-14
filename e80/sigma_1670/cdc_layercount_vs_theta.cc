#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFile.h"
#include "TBox.h"
#include "TLatex.h"

struct LayerInfo {
  int layer;
  std::string name;
  double radius; // mm
};

struct TrackSummary {
  bool entersCDC = false;
  double thetaDeg = 0.0;
  double sIn = 0.0;
  double sOut = 0.0;
  double flightLength = 0.0; // mm
  int nLayersHit = 0;
  int maxLayerHit = 0;
  std::vector<int> hitLayers;
};

static constexpr double kCDCInnerRadius = 150.0;   // mm
static constexpr double kCDCOuterRadius = 530.0;   // mm
static constexpr double kCDCHalfLength  = 1340.0;  // mm

std::vector<LayerInfo> BuildLayers()
{
  std::vector<LayerInfo> layers = {
    { 1, "A1-1(X)" , 190.5},
    { 2, "A1-2(X')", 204.0},
    { 3, "A1-3(X)" , 217.5},
    { 4, "U1-4(U)" , 248.5},
    { 5, "U1-5(U')", 262.0},
    { 6, "V1-6(V)" , 293.0},
    { 7, "V1-7(V')", 306.5},
    { 8, "A2-8(X)" , 337.5},
    { 9, "A2-9(X')", 351.0},
    {10, "U2-10(U)", 382.0},
    {11, "U2-11(U')",395.5},
    {12, "V2-12(V)", 426.5},
    {13, "V2-13(V')",440.0},
    {14, "A3-14(X)", 471.0},
    {15, "A3-15(X')",484.5}
  };
  return layers;
}

TrackSummary AnalyzeTrack(double thetaDeg,
                          const std::vector<LayerInfo>& layers)
{
  TrackSummary out;
  out.thetaDeg = thetaDeg;

  const double th = thetaDeg * M_PI / 180.0;
  const double st = std::sin(th);
  const double ct = std::cos(th);

  if (thetaDeg <= 0.0 || st <= 0.0 || ct <= 0.0) return out;

  // enter active CDC at inner hole boundary
  const double sIn = kCDCInnerRadius / st;
  const double zIn = sIn * ct;

  // if reaches endcap before inner radius, it never enters active region
  if (zIn > kCDCHalfLength) return out;

  const double sOutR = kCDCOuterRadius / st;
  const double sOutZ = kCDCHalfLength  / ct;
  const double sOut  = std::min(sOutR, sOutZ);

  if (sOut <= sIn) return out;

  out.entersCDC = true;
  out.sIn = sIn;
  out.sOut = sOut;
  out.flightLength = sOut - sIn;

  for (const auto& L : layers) {
    const double sLayer = L.radius / st;
    const double zLayer = sLayer * ct;

    if (zLayer <= kCDCHalfLength) {
      out.hitLayers.push_back(L.layer);
      out.nLayersHit++;
      out.maxLayerHit = L.layer;
    }
  }

  return out;
}

void PrintExampleAngles(const std::vector<LayerInfo>& layers)
{
  std::vector<double> angles = {6.5, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0};

  std::cout << "\n=== Example angles ===\n";
  for (double th : angles) {
    TrackSummary s = AnalyzeTrack(th, layers);

    std::cout << "theta = " << th << " deg";
    if (!s.entersCDC) {
      std::cout << "  -> no active CDC hit\n";
      continue;
    }

    std::cout << "  flightLength = " << s.flightLength << " mm"
              << "  nLayersHit = " << s.nLayersHit
              << "  maxLayer = " << s.maxLayerHit
              << "  hit layers: ";

    for (size_t i = 0; i < s.hitLayers.size(); ++i) {
      std::cout << s.hitLayers[i];
      if (i + 1 != s.hitLayers.size()) std::cout << ",";
    }
    std::cout << "\n";
  }
}

void DrawSideViewWithLayers(const std::vector<LayerInfo>& layers, double thetaDeg)
{
  TCanvas* c = new TCanvas("cSideLayers", "CDC side view with layers", 1200, 750);

  TH2D* frame = new TH2D("frame_side", ";z [mm];r [mm]",
                         100, -1500, 1500, 100, 0, 650);
  frame->SetStats(0);
  frame->Draw();

  // CDC active region boundaries
  TLine* l1 = new TLine(-kCDCHalfLength, kCDCInnerRadius,  kCDCHalfLength, kCDCInnerRadius);
  TLine* l2 = new TLine(-kCDCHalfLength, kCDCOuterRadius,  kCDCHalfLength, kCDCOuterRadius);
  TLine* l3 = new TLine(-kCDCHalfLength, kCDCInnerRadius, -kCDCHalfLength, kCDCOuterRadius);
  TLine* l4 = new TLine( kCDCHalfLength, kCDCInnerRadius,  kCDCHalfLength, kCDCOuterRadius);

  for (auto* line : {l1,l2,l3,l4}) {
    line->SetLineColor(kBlue+1);
    line->SetLineWidth(2);
    line->Draw("same");
  }

  // layer lines
  for (const auto& L : layers) {
    TLine* ll = new TLine(-kCDCHalfLength, L.radius, kCDCHalfLength, L.radius);
    ll->SetLineColor(kGray+1);
    ll->SetLineStyle(2);
    ll->Draw("same");
  }

  // vertex
  TBox* vtx = new TBox(-12, 0, 12, 12);
  vtx->SetFillColor(kGreen+2);
  vtx->SetLineColor(kGreen+2);
  vtx->Draw("same");

  // track
  const double th = thetaDeg * M_PI / 180.0;
  const double zEnd = 1450.0;
  const double rEnd = zEnd * std::tan(th);

  TLine* tr = new TLine(0.0, 0.0, zEnd, rEnd);
  tr->SetLineColor(kBlack);
  tr->SetLineWidth(2);
  tr->Draw("same");

  TrackSummary s = AnalyzeTrack(thetaDeg, layers);
  if (s.entersCDC) {
    const double rIn  = s.sIn  * std::sin(th);
    const double zIn  = s.sIn  * std::cos(th);
    const double rOut = s.sOut * std::sin(th);
    const double zOut = s.sOut * std::cos(th);

    TLine* inside = new TLine(zIn, rIn, zOut, rOut);
    inside->SetLineColor(kRed+1);
    inside->SetLineWidth(4);
    inside->Draw("same");

    for (const auto& L : layers) {
      double sLayer = L.radius / std::sin(th);
      double zLayer = sLayer * std::cos(th);
      if (zLayer <= kCDCHalfLength) {
        TBox* mark = new TBox(zLayer - 6, L.radius - 6, zLayer + 6, L.radius + 6);
        mark->SetFillColor(kMagenta+1);
        mark->SetLineColor(kMagenta+1);
        mark->Draw("same");
      }
    }
  }

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.15, 0.92, Form("#theta = %.1f^{#circ}", thetaDeg));
  lat.DrawLatex(0.15, 0.87, Form("Flight length = %.1f mm", s.flightLength));
  lat.DrawLatex(0.15, 0.82, Form("Layers hit = %d", s.nLayersHit));

  c->SaveAs("cdc_sideview_layers.pdf");
}

void cdc_layercount_vs_theta()
{
  gStyle->SetOptStat(0);

  auto layers = BuildLayers();

  const int n = 800;
  TGraph* grLength = new TGraph();
  TGraph* grCount  = new TGraph();
  TGraph* grMaxLay = new TGraph();

  int ip = 0;
  for (int i = 0; i < n; ++i) {
    double th = 0.0 + (40.0 / (n - 1)) * i; // 0 ~ 40 deg

    TrackSummary s = AnalyzeTrack(th, layers);

    grLength->SetPoint(ip, th, s.entersCDC ? s.flightLength : 0.0);
    grCount ->SetPoint(ip, th, s.nLayersHit);
    grMaxLay->SetPoint(ip, th, s.maxLayerHit);
    ip++;
  }

  // Critical angles
  const double thEnter = std::atan(kCDCInnerRadius / kCDCHalfLength) * 180.0 / M_PI;
  const double thOuter = std::atan(kCDCOuterRadius / kCDCHalfLength) * 180.0 / M_PI;

  // Canvas 1: flight length
  TCanvas* c1 = new TCanvas("c1", "flight length", 900, 700);
  grLength->SetTitle("Forward proton flight length in CDC;#theta_{lab}(p) [deg];Flight length [mm]");
  grLength->SetLineColor(kBlue+1);
  grLength->SetLineWidth(3);
  grLength->Draw("AL");

  TLine* lEnter1 = new TLine(thEnter, 0, thEnter, 1200);
  TLine* lOuter1 = new TLine(thOuter, 0, thOuter, 1200);
  lEnter1->SetLineColor(kRed+1);
  lOuter1->SetLineColor(kGreen+2);
  lEnter1->SetLineStyle(2);
  lOuter1->SetLineStyle(2);
  lEnter1->Draw("same");
  lOuter1->Draw("same");

  TLegend* leg1 = new TLegend(0.52, 0.72, 0.88, 0.88);
  leg1->AddEntry(grLength, "CDC flight length", "l");
  leg1->AddEntry(lEnter1, Form("enter CDC: %.2f deg", thEnter), "l");
  leg1->AddEntry(lOuter1, Form("reach outer wall: %.2f deg", thOuter), "l");
  leg1->Draw();

  c1->SaveAs("cdc_flightlength_vs_theta_layers.pdf");

  // Canvas 2: number of layers hit
  TCanvas* c2 = new TCanvas("c2", "layer count", 900, 700);
  grCount->SetTitle("Number of CDC layers crossed by forward proton;#theta_{lab}(p) [deg];Number of hit layers");
  grCount->SetLineColor(kMagenta+1);
  grCount->SetLineWidth(3);
  grCount->Draw("AL");

  TLine* lEnter2 = new TLine(thEnter, 0, thEnter, 16);
  TLine* lOuter2 = new TLine(thOuter, 0, thOuter, 16);
  lEnter2->SetLineColor(kRed+1);
  lOuter2->SetLineColor(kGreen+2);
  lEnter2->SetLineStyle(2);
  lOuter2->SetLineStyle(2);
  lEnter2->Draw("same");
  lOuter2->Draw("same");

  c2->SaveAs("cdc_layercount_vs_theta.pdf");

  // Canvas 3: outermost layer reached
  TCanvas* c3 = new TCanvas("c3", "max layer", 900, 700);
  grMaxLay->SetTitle("Outermost CDC layer reached by forward proton;#theta_{lab}(p) [deg];Max layer reached");
  grMaxLay->SetLineColor(kOrange+7);
  grMaxLay->SetLineWidth(3);
  grMaxLay->Draw("AL");

  TLine* lEnter3 = new TLine(thEnter, 0, thEnter, 16);
  TLine* lOuter3 = new TLine(thOuter, 0, thOuter, 16);
  lEnter3->SetLineColor(kRed+1);
  lOuter3->SetLineColor(kGreen+2);
  lEnter3->SetLineStyle(2);
  lOuter3->SetLineStyle(2);
  lEnter3->Draw("same");
  lOuter3->Draw("same");

  c3->SaveAs("cdc_maxlayer_vs_theta.pdf");

  // Side-view example
  DrawSideViewWithLayers(layers, 15.0);

  // Save root output
  TFile* fout = new TFile("cdc_layercount_vs_theta.root", "RECREATE");
  grLength->Write("grFlightLength");
  grCount->Write("grLayerCount");
  grMaxLay->Write("grMaxLayer");
  fout->Close();

  PrintExampleAngles(layers);

  std::cout << "\nSaved files:\n";
  std::cout << "  cdc_flightlength_vs_theta_layers.pdf\n";
  std::cout << "  cdc_layercount_vs_theta.pdf\n";
  std::cout << "  cdc_maxlayer_vs_theta.pdf\n";
  std::cout << "  cdc_sideview_layers.pdf\n";
  std::cout << "  cdc_layercount_vs_theta.root\n";
}
