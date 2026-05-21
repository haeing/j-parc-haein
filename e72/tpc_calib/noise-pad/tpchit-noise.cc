#include "../TPCPadHelper_260416.hh"

void tpchit_noise(){
  const int runnumber = 3396;
  //string dir = "/home/had/haein/data/JPARC2021May_root";
  string dir = "/hsm/had/sks/E72/JPARC2025Nov/rootfile/tpchit/v2";
  TFile *file = new TFile(Form("%s/run0%d_TPCHit.root",dir.c_str(),runnumber));
  TTree *tree = (TTree*)file->Get("tpc");

  int nhTpc;
  vector<int>* padTpc = nullptr;

  tree->SetBranchAddress("nhTpc",&nhTpc);
  tree->SetBranchAddress("padTpc",&padTpc);
  

  auto TPC_hit = new TH2Poly("TPC_hit","TPC_hit;Z;X",MinZ,MaxZ,MinX,MaxX);
  
  
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
      TPC_hit->AddBin(5, X, Y);
    }
  }

  //for(int n=0;n<50000;n++){
  for(int n=0;n<tree->GetEntries();n++){
    tree->GetEntry(n);
    cout<<nhTpc<<endl;
    for(int i=0;i<nhTpc;i++){
      double cnt = TPC_hit->GetBinContent((*padTpc)[i]+1);
      TPC_hit->SetBinContent((*padTpc)[i]+1,cnt+1);
    }
  }

  TFile *f = new TFile(Form("tpchit-noise-run0%d.root",runnumber),"RECREATE");
  
  auto c1 = new TCanvas("c1","c1");
  gPad->SetLogz();
  TPC_hit->Draw("colz");
  TPC_hit->Write();

  f->Close();

}
