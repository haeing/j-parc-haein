#include "../TPCPadHelper_260416.hh"

const int runnumber = 2599;
void tpctracking_padgain(){
  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/gain_calib";
  TFile *file = new TFile(Form("%s/run0%d_DstTPCHelixTracking.root",dir.c_str(),runnumber));
  TTree *tree = (TTree*)file->Get("tpc");

  int nclTpc;
  int ntTpc;
  
  vector<double>* cluster_x = nullptr;
  vector<double>* cluster_y = nullptr;
  vector<double>* cluster_z = nullptr;
  vector<double>* cluster_de = nullptr;
  vector<int>* cluster_size = nullptr;
  vector<int>* cluster_layer = nullptr;
  vector<int>* cluster_row_center = nullptr;
  vector<double>* cluster_mrow = nullptr;

  tree->SetBranchAddress("nclTpc",&nclTpc);
  tree->SetBranchAddress("ntTpc",&ntTpc);
  tree->SetBranchAddress("cluster_x",&cluster_x);
  tree->SetBranchAddress("cluster_y",&cluster_y);
  tree->SetBranchAddress("cluster_z",&cluster_z);
  tree->SetBranchAddress("cluster_de",&cluster_de);
  tree->SetBranchAddress("cluster_size",&cluster_size);
  tree->SetBranchAddress("cluster_layer",&cluster_layer);
  tree->SetBranchAddress("cluster_row_center",&cluster_row_center);
  tree->SetBranchAddress("cluster_mrow",&cluster_mrow);
  
  auto TPC_gain = new TH2Poly("TPC_gain", "TPC_gain;Z;X", MinZ, MaxZ, MinX, MaxX);
  TGraph *graph_gain = new TGraph();
  auto TPC_count = new TH2Poly("TPC_count", "TPC_count;Z;X", MinZ, MaxZ, MinX, MaxX);

  double l = (586./2.)/(1+sqrt(2.));
  Double_t px[9] = {-l*(1+sqrt(2.)),-l,l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.))};
  Double_t py[9] = {l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.)),-l,l};
  
  auto pLine = new TPolyLine(9,px,py);

  Double_t X[5];
  Double_t Y[5];
  for (Int_t l=0; l<NumOfLayersTPC; ++l) {
    Double_t pLength = padParameter[l][5];
    Double_t st      = (180.-(360./padParameter[l][3]) *
                        padParameter[l][1]/2.);
    Double_t sTheta  = (-1+st/180.)*TMath::Pi();
    Double_t dTheta  = (360./padParameter[l][3])/180.*TMath::Pi();
    Double_t cRad    = padParameter[l][2];
    Int_t    nPad    = padParameter[l][1];
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


  double de_sum[NumOfPadTPC]={0.};
  int count[NumOfPadTPC]={0};
  TH1D *hist_de[NumOfPadTPC];
  int total[NumOfLayersTPC]={0};
  int hit_layer[NumOfLayersTPC]={0};
  int clu_size[NumOfLayersTPC]={0};
  
  
  for(int i=0;i<NumOfPadTPC;i++){
    hist_de[i] = new TH1D(Form("hist_de%d",i),Form("hist_de%d",i),100,0,1000);
  }
  for(int n=0;n<tree->GetEntries();n++){
    tree->GetEntry(n);
    if(n%1000==0)std::cout<<n<<std::endl;
    if(ntTpc < 1)continue;
    for(int i=0;i<NumOfLayersTPC;i++)
      total[i]++;
    if(nclTpc < 1)continue;
    for(int ncl=0;ncl<nclTpc;ncl++){
      int layerid = (*cluster_layer)[ncl];
      int padid = GetPadId(layerid,(*cluster_row_center)[ncl]);
      if(Noise(padid))continue; 
      hit_layer[layerid]++;
      if(padid >=0){
	count[padid]++;
	de_sum[padid]+=(*cluster_de)[ncl];
	hist_de[padid]->Fill((*cluster_de)[ncl]);
      }
      
      
    }
  }
  

  gROOT->SetBatch(kTRUE);   
  TCanvas c("c","c",900,700);

  c.Print("padgain_fits.pdf[");

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
  TPC_count->Draw("colz");
  c.Print("padgain_fits.pdf");
  c.Print("padgain_fits.pdf]");   
   
}
