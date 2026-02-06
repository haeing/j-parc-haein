const int runnumber[7] = {2570, 2581, 2585, 2587, 2589, 2590, 2592};

const double z_tgt = -143.; //mm
const double z_tpccenter_to_ff = -150.; //mm

double DCA(double xbcout, double ybcout, double ubcout, double vbcout,
           double xtpc, double ytpc, double utpc, double vtpc)
{
  const TVector3 p1(xbcout, ybcout, 0.0);
  //TPC x0, y0 is at z = TGT -> Change to z = FF;
  const TVector3 p2(xtpc + utpc*(-z_tgt -z_tpccenter_to_ff) , ytpc + vtpc*(-z_tgt - z_tpccenter_to_ff), 0.0);
  
  const TVector3 v1vec(ubcout, vbcout, 1.0);
  const TVector3 v2vec(utpc, vtpc, 1.0);

  const TVector3 w0 = p1 - p2;

  const double a = v1vec.Dot(v1vec);
  const double b = v1vec.Dot(v2vec);
  const double c = v2vec.Dot(v2vec);
  const double d = v1vec.Dot(w0);
  const double e = v2vec.Dot(w0);

  const double denom = a*c - b*b;
  const double eps = 1e-12;

  double s = 0.0;
  double t = 0.0;

  if (TMath::Abs(denom) > eps) {
    // general (not parallel)
    s = (b*e - c*d) / denom;
    t = (a*e - b*d) / denom;
  } else {
    // nearly parallel: take t=0 and minimize distance to line2 point p2
    // s = argmin |(p1 + s*v1) - p2|^2
    s = -d / a;
    t = 0.0;
  }

  const TVector3 cp1 = p1 + s * v1vec;
  const TVector3 cp2 = p2 + t * v2vec;

  /*
  if (c1) *c1 = cp1;
  if (c2) *c2 = cp2;
  if (s_out) *s_out = s;
  if (t_out) *t_out = t;
  */

  return (cp1 - cp2).Mag();
}

double cos_theta(double ubcout, double vbcout, double utpc, double vtpc)
{
    const TVector3 d1(ubcout, vbcout, 1.0);
    const TVector3 d2(utpc, vtpc, 1.0);

    const double n1 = d1.Mag();
    const double n2 = d2.Mag();
    if (n1 == 0.0 || n2 == 0.0) return -2.0; // invalid marker

    return d1.Dot(d2) / (n1 * n2);
}


void hsoff_beamthrough_padgain(){

  TH1D *hist_dca = new TH1D("hist_dca","hist_dca",1000,0,1000);
  TH1D *hist_cos = new TH1D("hist_cos","hist_cos",100,-1,1);
  TH2D *hist_dca_cos = new TH2D("hist_dca_cos","hist_dca_cos",1000,0,1000,100,-1,1);



  TH1D *hist_res_x = new TH1D("hist_res_x","hist_res_x",100,-50,50);
  TH1D *hist_res_y = new TH1D("hist_res_y","hist_res_y",100,40,80);
  TH2D *hist_res = new TH2D("hist_res","hist_res",100,-100,100,100,0,120);

  TH1D *hist_cor_track = new TH1D("hist_cor_track","hist_cor_track",6,-0.5,5.5);
  
  TH2D *hist_bcout_pass = new TH2D("hist_bcout_pass","hist_bcout_pass",300,-100,100,300,-100,100);
  TH2D *hist_bcout_fail = new TH2D("hist_bcout_fail","hist_bcout_fail",300,-100,100,300,-100,100);

  TH2D *hist_bcout_uv_pass = new TH2D("hist_bcout_uv_pass","hist_bcout_uv_pass",300,-0.5,0.5,300,-0.5,0.5);
  TH2D *hist_bcout_uv_fail = new TH2D("hist_bcout_uv_fail","hist_bcout_uv_fail",300,-0.5,0.5,300,-0.5,0.5);
  
  for (size_t i = 0; i < 1; ++i) {
    int run = runnumber[i];
    //for(int run : runnumber){
  
    string dir = "/gpfs/group/had/sks/Users/haein/JPARC2025Nov_root/hsoff_beamthrough";
    TFile *file_hodo = new TFile(Form("%s/run0%d_Hodoscope.root",dir.c_str(),run));
    TFile *file_bcout = new TFile(Form("%s/run0%d_BcOutTracking.root",dir.c_str(),run));
    TFile *file_tpc = new TFile(Form("%s/run0%d_DstTPCTracking.root",dir.c_str(),run));

    TTree *tree_hodo = (TTree*)file_hodo->Get("hodo");
    TTree *tree_bcout = (TTree*)file_bcout->Get("bcout");
    TTree *tree_tpc = (TTree*)file_tpc->Get("tpc");

    //Hodo
    double btof0;
    tree_hodo->SetBranchAddress("btof0",&btof0);

    //BcOut
    int ntrack;
    vector<double>* chisqr = nullptr;
    vector<double>* x0 = nullptr;
    vector<double>* y0 = nullptr;
    vector<double>* u0 = nullptr;
    vector<double>* v0 = nullptr;
    tree_bcout->SetBranchAddress("ntrack",&ntrack);
    tree_bcout->SetBranchAddress("chisqr",&chisqr);
    tree_bcout->SetBranchAddress("x0",&x0);
    tree_bcout->SetBranchAddress("y0",&y0);
    tree_bcout->SetBranchAddress("u0",&u0);
    tree_bcout->SetBranchAddress("v0",&v0);

    //TPCTracking
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

    vector<int>* nhtrack = nullptr;
    vector<double>* x0Tpc = nullptr;
    vector<double>* y0Tpc = nullptr;
    vector<double>* u0Tpc = nullptr;
    vector<double>* v0Tpc = nullptr;


    tree_tpc->SetBranchAddress("nclTpc",&nclTpc);
    tree_tpc->SetBranchAddress("ntTpc",&ntTpc);
    tree_tpc->SetBranchAddress("cluster_x",&cluster_x);
    tree_tpc->SetBranchAddress("cluster_y",&cluster_y);
    tree_tpc->SetBranchAddress("cluster_z",&cluster_z);
    tree_tpc->SetBranchAddress("cluster_de",&cluster_de);
    tree_tpc->SetBranchAddress("cluster_size",&cluster_size);
    tree_tpc->SetBranchAddress("cluster_layer",&cluster_layer);
    tree_tpc->SetBranchAddress("cluster_row_center",&cluster_row_center);
    tree_tpc->SetBranchAddress("cluster_mrow",&cluster_mrow);

    tree_tpc->SetBranchAddress("nhtrack",&nhtrack);
    tree_tpc->SetBranchAddress("x0Tpc",&x0Tpc);
    tree_tpc->SetBranchAddress("y0Tpc",&y0Tpc);
    tree_tpc->SetBranchAddress("u0Tpc",&u0Tpc);
    tree_tpc->SetBranchAddress("v0Tpc",&v0Tpc);


    bool IsPion = false;
    bool IsTPCTrack = false;
    bool IsBcoutTrack = false;
      
    for(int n=0;n<tree_hodo->GetEntries();n++){
    //for(int n=0;n<100;n++){
      tree_hodo->GetEntry(n);
      tree_bcout->GetEntry(n);
      tree_tpc->GetEntry(n);

      IsPion = false;
      IsTPCTrack = false;
      IsBcoutTrack = false;
      int IsPionTrack = 0;
      
      if(btof0 > -3)IsPion = true;
      //Tight cut to BcOut 
      if(ntrack == 1){
	if((*chisqr)[0] < 3.)
	  IsBcoutTrack = true;
      }
      if(ntTpc > 0)IsTPCTrack = true;

      if(IsPion && IsBcoutTrack && IsTPCTrack){
	//Find Beam Track
	IsPionTrack = 0;
	for(int nbcouttr = 0;nbcouttr<ntrack;nbcouttr++){
	  for(int ntpctr = 0;ntpctr<ntTpc;ntpctr++){
	    double min_dis = DCA((*x0)[nbcouttr],(*y0)[nbcouttr],(*u0)[nbcouttr],(*v0)[nbcouttr],(*x0Tpc)[ntpctr],(*y0Tpc)[ntpctr],(*u0Tpc)[ntpctr],(*v0Tpc)[ntpctr]);
	    hist_dca->Fill(min_dis);
	    double track_angle = cos_theta((*u0)[nbcouttr],(*v0)[nbcouttr],(*u0Tpc)[ntpctr],(*v0Tpc)[ntpctr]);
	    hist_cos->Fill(track_angle);
	    hist_dca_cos->Fill(min_dis,track_angle);


	    double xbcout_ff = (*x0)[nbcouttr];
	    double ybcout_ff = (*y0)[nbcouttr];
	    double xtpc_ff = (*x0Tpc)[ntpctr] + (*u0Tpc)[ntpctr]*(-z_tgt -z_tpccenter_to_ff);
	    double ytpc_ff = (*y0Tpc)[ntpctr] + (*v0Tpc)[ntpctr]*(-z_tgt -z_tpccenter_to_ff);

	    hist_res_x->Fill(xbcout_ff - xtpc_ff);
	    hist_res_y->Fill(ybcout_ff - ytpc_ff);
	    hist_res->Fill(xbcout_ff - xtpc_ff,ybcout_ff - ytpc_ff);

	    double res_x = xbcout_ff - xtpc_ff;
	    double res_y = ybcout_ff - ytpc_ff;

	    if(res_x > -20 && res_x < 20 && res_y > 50 && res_y < 80)IsPionTrack++;
	    
	  }
	}
	hist_cor_track->Fill(IsPionTrack);
	if(IsPionTrack ==0){
	  hist_bcout_fail->Fill((*x0)[0],(*y0)[0]);
	  hist_bcout_uv_fail->Fill((*u0)[0],(*v0)[0]);
	}
	else if(IsPionTrack >0){
	  hist_bcout_pass->Fill((*x0)[0],(*y0)[0]);
	  hist_bcout_uv_pass->Fill((*u0)[0],(*v0)[0]);
	}
      }

      
      
      
            
    }
  }

  auto c1 = new TCanvas("c1","c1");
  c1->Divide(3);
  c1->cd(1);
  hist_dca->Draw();
  c1->cd(2);
  hist_cos->Draw();
  c1->cd(3);
  hist_dca_cos->Draw("colz");

  auto c2 = new TCanvas("c2","c2");
  c2->Divide(2,2);
  c2->cd(1);
  hist_res_x->Draw();
  c2->cd(2);
  hist_res_y->Draw();
  c2->cd(3);
  hist_res->Draw("colz");
  c2->cd(4);
  hist_cor_track->Draw();

  auto c3 = new TCanvas("c3","c3");
  c3->Divide(2,2);
  c3->cd(1);
  hist_bcout_pass->Draw("colz");
  c3->cd(2);
  hist_bcout_fail->Draw("colz");
  c3->cd(3);
  hist_bcout_uv_pass->Draw("colz");
  c3->cd(4);
  hist_bcout_uv_fail->Draw("colz");
  
  
    
}
