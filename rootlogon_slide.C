{
  TStyle *style = new TStyle("MyStyle","My ROOT Style");

  // ===== Canvas =====
  style->SetCanvasDefW(1000);
  style->SetCanvasDefH(800);

  // ===== Pad =====
  style->SetPadTopMargin(0.05);
  style->SetPadBottomMargin(0.13);
  style->SetPadLeftMargin(0.13);
  style->SetPadRightMargin(0.16);

  // ===== Font (Times New Roman 느낌) =====
  /*
  style->SetTextFont(42);
  style->SetLabelFont(42, "XYZ");
  style->SetTitleFont(42, "XYZ");
  */

  
  // ===== Axis Title =====
  style->SetTitleSize(0.06, "XYZ");   // 🔥 크게
  style->SetTitleOffset(0.9, "X");
  style->SetTitleOffset(1.0, "Y");

  // ===== Axis Label (숫자) =====
  style->SetLabelSize(0.05, "XYZ");   // 🔥 크게

  // ===== Tick =====
  style->SetPadTickX(1);
  style->SetPadTickY(1);

  //style->SetPadGridX(true);
  //style->SetPadGridY(true);
  style->SetGridColor(kGray+1);
  // ===== Line =====
  style->SetFrameLineWidth(2);

  // ===== Legend =====
  style->SetLegendBorderSize(0);
  style->SetLegendFont(42);
  style->SetLegendTextSize(0.05);

  // ===== Stat box 제거 =====
  style->SetOptStat(0);

  // ===== Title 제거 =====
  style->SetOptTitle(0);
  style->SetCanvasColor(0);
  style->SetPadColor(0);
  style->SetFrameFillColor(0);
  style->SetStatColor(0);
  style->SetLegendFillColor(0);

  // ===== Apply =====
  gROOT->SetStyle("MyStyle");
  gROOT->ForceStyle();
}
