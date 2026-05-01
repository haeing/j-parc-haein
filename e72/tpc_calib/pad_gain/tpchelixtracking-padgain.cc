#include "../TPCPadHelper_260416.hh"

const int runnumber = 2599;
void tpchelixtracking_padgain(){

  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/gain_calib";
  TFile *file = new TFile(Form("%s/run0%d_DstTPCHelixTracking.root",dir.c_str(),runnumber));
  TTree *tree = (TTree*)file->Get("tpc");

  int nclTpc;
  int ntTpc;

  vector<double>* mom0 = nullptr;
  vector<double>* dEdx = nullptr;

  //Track Info
  vector<int>* nhtrack = nullptr;
  vector<vector<double>>* hitlayer = nullptr;
  vector<vector<double>>* track_cluster_de = nullptr;
  vector<vector<double>>* track_cluster_size = nullptr;
  vector<vector<double>>* track_cluster_mrow = nullptr;
  vector<vector<double>>* track_cluster_row_center = nullptr;
  
  
  tree->SetBranchAddress("ntTpc",&ntTpc);
  tree->SetBranchAddress("mom0",&mom0);
  tree->SetBranchAddress("dEdx",&dEdx);
  
  tree->SetBranchAddress("nhtrack",&nhtrack);
  tree->SetBranchAddress("hitlayer",&hitlayer);
  tree->SetBranchAddress("track_cluster_de",&track_cluster_de);
  tree->SetBranchAddress("track_cluster_size",&track_cluster_size);
  tree->SetBranchAddress("track_cluster_mrow",&track_cluster_mrow);
  tree->SetBranchAddress("track_cluster_row_center",&track_cluster_row_center);
  
  
  auto TPC_gain = new TH2Poly("TPC_gain", "TPC_gain;Z;X", MinZ, MaxZ, MinX, MaxX);
  TGraph *graph_gain = new TGraph();
  auto TPC_count = new TH2Poly("TPC_count", "TPC_count;Z;X", MinZ, MaxZ, MinX, MaxX);
  auto TPC_pid = new TH2D("TPC_pid","TPC_pid;p [GeV/#it{c}];dE/dx [A.U.]",100,0,1,100,0,500);

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
      TPC_gain->AddBin(5, X, Y);
      TPC_count->AddBin(5, X, Y);
    }
  }

  for(int i=0;i<NumOfPadTPC;i++){
    TPC_count->SetBinContent(i+1,0);
  }

  double de_sum[NumOfPadTPC]={0.};
  int count[NumOfPadTPC]={0};
  TH1D *hist_de[NumOfPadTPC];
  TH1D *hist_de_cut[NumOfPadTPC];
  int total[NumOfLayersTPC]={0};
  int hit_layer[NumOfLayersTPC]={0};
  int clu_size[NumOfLayersTPC]={0};
  TH1D *hist_clu = new TH1D("hist_clu","hist_clu;Cluster;Entries",5,-0.5,4.5);
  TH1D *hist_clu_mipcut = new TH1D("hist_clu_mipcut","hist_clu_mipcut",5,-0.5,4.5);
  
  
  for(int i=0;i<NumOfPadTPC;i++){
    hist_de[i] = new TH1D(Form("hist_de%d",i),Form("hist_de%d",i),100,0,1000);
    hist_de_cut[i] = new TH1D(Form("hist_de_cut%d",i),Form("hist_de_cut%d",i),100,0,1000);
  }

  for(int n=0;n<tree->GetEntries();n++){
    //for(int n=0;n<50000;n++){
    tree->GetEntry(n);
    if(n%1000==0)std::cout<<n<<std::endl;
    if(ntTpc < 1)continue;

    //MIP cut
    for(int ntr = 0;ntr<ntTpc;ntr++){
      TPC_pid->Fill((*mom0)[ntr],(*dEdx)[ntr]);
      //Check if all pads have hit (w/ enough statistics)
      for(int nhit = 0;nhit<(*nhtrack)[ntr];nhit++){
	hist_clu->Fill((*track_cluster_size)[ntr][nhit]);
	if((*dEdx)[ntr]<50)hist_clu_mipcut->Fill((*track_cluster_size)[ntr][nhit]);
	int layerid = (*hitlayer)[ntr][nhit];
	int padid = tpc::GetPadId(layerid,(*track_cluster_row_center)[ntr][nhit]);
	double cnt = TPC_count->GetBinContent(padid+1);
	hist_de[padid]->Fill((*track_cluster_de)[ntr][nhit]);
	//if((*dEdx)[ntr]<50 && (*track_cluster_size)[ntr][nhit]==1){
	if((*dEdx)[ntr]<50){
	  TPC_count->SetBinContent(padid+1,cnt+1.);
	  hist_de_cut[padid]->Fill((*track_cluster_de)[ntr][nhit]);
	}
	
      }
    }
  }

  
  TFile *f = new TFile("padgain.root","RECREATE");
  
  gROOT->SetBatch(kTRUE);   
  TCanvas c("c","c",900,700);
  gStyle->SetOptStat(0);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("tpchelixtracking-padgain.cc");
  p->AddText("TPC Pad Gain Calibration");
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c.Print("padgain_fits.pdf(");
  c.Clear();
  gPad->SetLogz();
  TPC_pid->Draw("colz");
  c.Print("padgain_fits.pdf");

  c.Clear();
  hist_clu->Draw();
  hist_clu_mipcut->SetLineColor(kRed);
  hist_clu_mipcut->Draw("same");
  c.Print("padgain_fits.pdf");
  
  /*
  auto g_eff = new TGraph();
  for(int i=0;i<NumOfLayersTPC;i++){
    std::cout<<i<<"'s # of hit : "<<hit_layer[i]<<std::endl;
    g_eff->AddPoint(i,(double)hit_layer[i] / (double) total[i]);
  }
  c.Print("padgain_fits.pdf");
  
  for(int ipad = 0;ipad <NumOfPadTPC;ipad++){
  TF1 fL("fL", "landau", 0, 1000);
    fL.SetParLimits(1,50,300);
    if(count[ipad] == 0)continue;
    hist_de[ipad]->Fit(&fL, "QR");
    auto fit_param = fL.GetParameters();
    
    TPC_gain->SetBinContent(ipad+1,de_sum[ipad] / (double)count[ipad]);
    TPC_count->SetBinContent(ipad+1,count[ipad]);
    graph_gain->AddPoint(ipad,fit_param[1]);
    if (ipad % 100 == 0) {
      c.Clear();
      hist_de[ipad]->Draw();
      c.Print("padgain_fits.pdf");
    }
  }


  
  //TPC_gain->Draw("colz");
  graph_gain->Draw();
  c.Print("padgain_fits.pdf");
  //auto c1 = new TCanvas("c1","c1");
  */

  for(int ipad = 0;ipad <NumOfPadTPC;ipad++){
    TF1 fL("fL", "landau", 0, 1000);
    fL.SetParLimits(1,50,300);

    if(ipad%30 == 0){
      if(ipad !=0)c.Print("padgain_fits.pdf");
      c.Clear();
      c.Divide(6,5);
    }
    c.cd(ipad%30+1);
    
    hist_de_cut[ipad]->Fit(&fL, "QR");
    auto fit_param = fL.GetParameters();

    graph_gain->AddPoint(ipad,fit_param[1]);
    TPC_gain->SetBinContent(ipad+1,fit_param[1]);
    hist_de[ipad]->Draw();
    hist_de_cut[ipad]->SetLineColor(kRed);
    hist_de_cut[ipad]->Draw("same");

    
    hist_de[ipad]->Write();
    hist_de_cut[ipad]->Write();
  }
  
  f->Close();
  c.Print("padgain_fits.pdf");
  c.Clear();
  graph_gain->SetMarkerStyle(4);
  graph_gain->Draw("AP");
  c.Print("padgain_fits.pdf");
  c.Clear();
  TPC_gain->Draw("colz");
  c.Print("padgain_fits.pdf");
  c.Clear();
  TPC_count->Draw("colz");
  c.Print("padgain_fits.pdf)");   
   
}
