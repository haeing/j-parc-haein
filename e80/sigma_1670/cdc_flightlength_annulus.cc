#include <iostream>
#include <cmath>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TArc.h"
#include "TStyle.h"
#include "TLegend.h"

struct CylCDC {
  double rin;    // inner radius [mm]
  double rout;   // outer radius [mm]
  double zhalf;  // half length [mm]
};

struct TrackResult {
  bool enters = false;
  double s_in = 0.0;
  double s_out = 0.0;
  double length = 0.0;
  double r_in = 0.0, z_in = 0.0;
  double r_out = 0.0, z_out = 0.0;
};

TrackResult FlightLengthInCDC(double thetaDeg, const CylCDC& cdc)
{
  TrackResult res;

  const double th = thetaDeg * M_PI / 180.0;
  const double st = std::sin(th);
  const double ct = std::cos(th);

  if (thetaDeg <= 0.0 || st <= 0.0 || ct <= 0.0) return res;

  // enter CDC active annulus at r = rin
  const double s_in = cdc.rin / st;
  const double z_at_in = s_in * ct;

  // if it reaches endcap before inner radius, never enters CDC active region
  if (z_at_in > cdc.zhalf) return res;

  const double s_out_r = cdc.rout / st;   // outer cylindrical wall
  const double s_out_z = cdc.zhalf / ct;  // forward endcap

  const double s_out = std::min(s_out_r, s_out_z);

  if (s_out <= s_in) return res;

  res.enters = true;
  res.s_in = s_in;
  res.s_out = s_out;
  res.length = s_out - s_in;

  res.r_in = cdc.rin;
  res.z_in = s_in * ct;

  res.r_out = s_out * st;
  res.z_out = s_out * ct;

  return res;
}

void DrawSideView(const CylCDC& cdc, double thetaDeg)
{
  TCanvas* c = new TCanvas("cSide", "CDC side view", 1100, 700);

  TH2D* frame = new TH2D("frame", ";z [mm];r [mm]",
                         100, -1500, 1500, 100, 0, 700);
  frame->SetStats(0);
  frame->Draw();

  // CDC boundaries in side view
  TLine* l1 = new TLine(-cdc.zhalf, cdc.rin,  cdc.zhalf, cdc.rin);
  TLine* l2 = new TLine(-cdc.zhalf, cdc.rout, cdc.zhalf, cdc.rout);
  TLine* l3 = new TLine(-cdc.zhalf, cdc.rin, -cdc.zhalf, cdc.rout);
  TLine* l4 = new TLine( cdc.zhalf, cdc.rin,  cdc.zhalf, cdc.rout);

  l1->SetLineColor(kBlue+1); l1->SetLineWidth(3);
  l2->SetLineColor(kBlue+1); l2->SetLineWidth(3);
  l3->SetLineColor(kBlue+1); l3->SetLineWidth(3);
  l4->SetLineColor(kBlue+1); l4->SetLineWidth(3);

  l1->Draw("same"); l2->Draw("same"); l3->Draw("same"); l4->Draw("same");

  // vertex
  TLine* vx1 = new TLine(-15, 0, 15, 0);
  vx1->SetLineColor(kGreen+2);
  vx1->SetLineWidth(4);
  vx1->Draw("same");

  const double th = thetaDeg * M_PI / 180.0;
  const double zEnd = 1450.0;
  const double rEnd = zEnd * std::tan(th);

  TLine* tr = new TLine(0.0, 0.0, zEnd, rEnd);
  tr->SetLineColor(kBlack);
  tr->SetLineWidth(2);
  tr->Draw("same");

  TrackResult rr = FlightLengthInCDC(thetaDeg, cdc);
  if (rr.enters) {
    TLine* inside = new TLine(rr.z_in, rr.r_in, rr.z_out, rr.r_out);
    inside->SetLineColor(kRed+1);
    inside->SetLineWidth(4);
    inside->Draw("same");

    std::cout << "theta = " << thetaDeg << " deg"
              << "  z_in = " << rr.z_in << " mm"
              << "  z_out = " << rr.z_out << " mm"
              << "  L = " << rr.length << " mm\n";
  } else {
    std::cout << "theta = " << thetaDeg << " deg : no active CDC hit\n";
  }

  c->SaveAs("cdc_annulus_sideview.pdf");
}

void DrawCrossSection(const CylCDC& cdc)
{
  TCanvas* c = new TCanvas("cCross", "CDC cross section", 700, 700);

  TH2D* frame = new TH2D("frame2", ";x [mm];y [mm]",
                         100, -650, 650, 100, -650, 650);
  frame->SetStats(0);
  frame->Draw();

  TEllipse* outer = new TEllipse(0, 0, cdc.rout, cdc.rout);
  TEllipse* inner = new TEllipse(0, 0, cdc.rin,  cdc.rin);

  outer->SetFillStyle(0);
  inner->SetFillStyle(0);
  outer->SetLineColor(kBlue+1);
  inner->SetLineColor(kBlue+1);
  outer->SetLineWidth(3);
  inner->SetLineWidth(3);

  outer->Draw("same");
  inner->Draw("same");

  c->SaveAs("cdc_annulus_crosssection.pdf");
}

void cdc_flightlength_annulus()
{
  gStyle->SetOptStat(0);

  CylCDC cdc;
  cdc.rin   = 150.0;   // mm
  cdc.rout  = 530.0;   // mm
  cdc.zhalf = 1340.0;  // mm (2680 / 2)

  const int n = 600;
  TGraph* gr = new TGraph();

  int ip = 0;
  for (int i = 0; i < n; ++i) {
    double th = 0.0 + (40.0 / (n - 1)) * i; // 0~40 deg
    TrackResult rr = FlightLengthInCDC(th, cdc);
    double L = rr.enters ? rr.length : 0.0;
    gr->SetPoint(ip++, th, L);
  }

  TCanvas* c1 = new TCanvas("c1", "Flight length vs theta", 900, 700);
  gr->SetTitle("Forward proton flight length inside cylindrical CDC;#theta_{lab}(p) [deg];Flight length in CDC [mm]");
  gr->SetLineColor(kBlue+1);
  gr->SetLineWidth(3);
  gr->Draw("AL");

  // mark critical angles
  const double th_enter = atan(cdc.rin / cdc.zhalf) * 180.0 / M_PI;
  const double th_outer = atan(cdc.rout / cdc.zhalf) * 180.0 / M_PI;

  TLine* le1 = new TLine(th_enter, 0, th_enter, gr->GetYaxis()->GetXmax());
  TLine* le2 = new TLine(th_outer, 0, th_outer, gr->GetYaxis()->GetXmax());
  le1->SetLineStyle(2); le2->SetLineStyle(2);
  le1->SetLineColor(kRed+1); le2->SetLineColor(kGreen+2);
  le1->Draw("same");
  le2->Draw("same");

  TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
  leg->AddEntry(gr, "CDC flight length", "l");
  leg->AddEntry(le1, Form("enter active CDC: %.2f deg", th_enter), "l");
  leg->AddEntry(le2, Form("reach outer radius: %.2f deg", th_outer), "l");
  leg->Draw();

  c1->SaveAs("cdc_flightlength_vs_theta_annulus.pdf");

  DrawSideView(cdc, 12.0);
  DrawCrossSection(cdc);

  std::cout << "critical angle to enter active CDC = " << th_enter << " deg\n";
  std::cout << "critical angle to reach outer wall = " << th_outer << " deg\n";
}