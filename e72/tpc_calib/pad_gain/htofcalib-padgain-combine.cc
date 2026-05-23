#include "../TPCPadHelper_260416.hh"

const int runnumber[8] = {2489,2599, 2601, 2602, 2603, 2604, 2606, 2607};

void htofcalib_padgain_combine(){
  gROOT->SetBatch(kTRUE);
  string outpdf = "htofcalib-padgain-combine-minlayer5.pdf";
  
  TH1D *hist_de[NumOfPadTPC];
  /*
  for(int i=0;i<NumOfPadTPC;i++){
    hist_de[i] = new TH1D(Form("hist_de%d",i),Form("hist_de%d",i),100,0,1000);
  }
  */
  TH2Poly *TPC_count = new TH2Poly("TPC_count","TPC_count;Z;X",MinZ,MaxZ,MinX,MaxX);
  auto TPC_gain = new TH2Poly("TPC_gain","TPC_gain;Z;X",MinZ,MaxZ,MinX,MaxX);
  TGraph *graph_gain = new TGraph();
  graph_gain->SetName("graph_gain");

  double l = (586./2.)/(1+sqrt(2.));
  Double_t px[9] = {-l*(1+sqrt(2.)),-l,l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.))};
  Double_t py[9] = {l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.)),-l,l};
  
  auto pLine = new TPolyLine(9,px,py);

  Double_t X[5];
  Double_t Y[5];
  for (Int_t l=0; l<NumOfLayersTPC; ++l) {
    Double_t pLength = tpc::padParameter[l][5];
    Double_t st      = (180.-(360./tpc::padParameter[l][3]) *
                        tpc::padParameter[l][1]/2.);
    Double_t sTheta  = (-1+st/180.)*TMath::Pi();
    Double_t dTheta  = (360./tpc::padParameter[l][3])/180.*TMath::Pi();
    Double_t cRad    = tpc::padParameter[l][2];
    Int_t    nPad    = tpc::padParameter[l][1];
    for (Int_t j=0; j<nPad; ++j) {
      X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[0] = X[4];
      Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[0] = Y[4];
      for (Int_t k=0; k<5; ++k) X[k] += ZTarget;
      TPC_count->AddBin(5, X, Y);
      TPC_gain->AddBin(5, X, Y);
    }
  }
  
  for(int i=0;i<8;i++){
    TFile *file = new TFile(Form("htofcalib-padgain-run0%d-minlayer5.root",runnumber[i]));
    TH2Poly *TPC_tr_cluster = (TH2Poly*)file->Get("TPC_tr_cluster");
    
    for(int n=0;n<TPC_tr_cluster->GetNumberOfBins();n++){
      double run_cnt = TPC_tr_cluster->GetBinContent(n+1);

      double total_cnt = TPC_count->GetBinContent(n+1);
      TPC_count->SetBinContent(n+1,total_cnt + run_cnt);

      auto hist = (TH1D*)file->Get(Form("hist_de%d",n));
      if(i==0){
        hist_de[n] = (TH1D*)hist->Clone(Form("hist_de%d",n));
	hist_de[n] ->SetDirectory(0);
      }
      else{
        hist_de[n]->Add(hist);
      }
    }
  }


  TFile *f = new TFile("htofcalib-padgain-combine-minlayer5.root","RECREATE");
  
  auto c1 = new TCanvas("c1","c1");
  gStyle->SetOptStat(0);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("htofcalib-padgain-combine.cc");
  p->AddText("TPC Pad Gain Calibration");
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());




  for(int ipad = 0;ipad <NumOfPadTPC;ipad++){
    TF1 fL("fL", "landau", 100, 700);
    fL.SetParLimits(1,10,500);
    fL.SetLineColor(kRed);
    
    if(ipad%30 == 0){
      if(ipad !=0)c1->Print(outpdf.c_str());
      c1->Clear();
      c1->Divide(6,5);
    }
    c1->cd(ipad%30+1);
    if(hist_de[ipad]->GetEntries()>80){
      hist_de[ipad]->Fit(&fL, "QR");
      auto fit_param = fL.GetParameters();
      graph_gain->AddPoint(ipad,fit_param[1]);
      TPC_gain->SetBinContent(ipad+1,fit_param[1]);
    }

    hist_de[ipad]->Draw();
    hist_de[ipad]->Write();
  }

  c1->Print(outpdf.c_str());
  
  c1->Clear();
  TPC_gain->SetMinimum(0);
  TPC_gain->SetMaximum(400);
  TPC_gain->Draw("colz");
  TPC_gain->Write();
  c1->Print(outpdf.c_str());

  c1->Clear();
  graph_gain->SetMarkerStyle(4);
  graph_gain->Draw("AP");
  graph_gain->Write();
  c1->Print(outpdf.c_str());

  
  
  c1->Clear();
  gPad->SetLogz();
  TPC_count->Draw("colz");
  TPC_count->Write();
  c1->Print((outpdf + ")").c_str());
  f->Close();
}
