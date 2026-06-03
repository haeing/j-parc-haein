const int NumOfSegBH2 = 15;
const int NumOfSegBHT = 63;
const double tdc_min = 708000;
const double tdc_max = 714000;

const double tdc_window_min = 710000;
const double tdc_window_max = 713000;

const double tdc_bht_window_min = 1400000;
const double tdc_bht_window_max = 1450000;

void bh2_blc(){
  int run = 2486;
  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root";
  TFile *file_hodo = new TFile(Form("%s/run0%d_Hodoscope.root",dir.c_str(),run));
  //TFile *file_bcin = new TFile(Form("%s/run0%d_BcInTracking_bcin13.root",dir.c_str(),run));
  TFile *file_bcin = new TFile("~/work/e72/ana/e72/test_bcin.root");
  TFile *file_bcout = new TFile(Form("%s/run0%d_BcOutTracking.root",dir.c_str(),run));

  TTree *tree_hodo = (TTree*)file_hodo->Get("hodo");
  TTree *tree_bcin = (TTree*)file_bcin->Get("bcin");
  TTree *tree_bcout = (TTree*)file_bcout->Get("bcout");


  //Hodo
  vector<double>* bh2_raw_seg = nullptr;
  vector<double>* bh2_adc_u = nullptr;
  vector<double>* bh2_adc_d = nullptr;
  vector<vector<double>>* bh2_tdc_u = nullptr;
  vector<vector<double>>* bh2_tdc_d = nullptr;
  vector<vector<double>>* bh2_tdc_s = nullptr;
  vector<double>* bht_raw_seg = nullptr;
  vector<vector<double>>* bht_tdc_u = nullptr;
  vector<vector<double>>* bht_tdc_d = nullptr;

  

  bh2_tdc_u = new vector<vector<double>>();
  bh2_tdc_d = new vector<vector<double>>();
  bh2_tdc_s = new vector<vector<double>>();
  bht_tdc_u = new vector<vector<double>>();
  bht_tdc_d = new vector<vector<double>>();
  
  tree_hodo->SetBranchAddress("bh2_raw_seg",&bh2_raw_seg);
  tree_hodo->SetBranchAddress("bh2_adc_u",&bh2_adc_u);
  tree_hodo->SetBranchAddress("bh2_adc_d",&bh2_adc_d);
  tree_hodo->SetBranchAddress("bh2_tdc_u",&bh2_tdc_u);
  tree_hodo->SetBranchAddress("bh2_tdc_d",&bh2_tdc_d);
  tree_hodo->SetBranchAddress("bh2_tdc_s",&bh2_tdc_s);

  tree_hodo->SetBranchAddress("bht_raw_seg",&bht_raw_seg);
  tree_hodo->SetBranchAddress("bht_tdc_u",&bht_tdc_u);
  tree_hodo->SetBranchAddress("bht_tdc_d",&bht_tdc_d);
  
  
  //BcOut
  int ntrack_out;
  vector<double>* chisqr_out = nullptr;
  vector<double>* x0_out = nullptr;
  vector<double>* y0_out = nullptr;
  vector<double>* u0_out = nullptr;
  vector<double>* v0_out = nullptr;
  tree_bcout->SetBranchAddress("ntrack",&ntrack_out);
  tree_bcout->SetBranchAddress("chisqr",&chisqr_out);
  tree_bcout->SetBranchAddress("x0",&x0_out);
  tree_bcout->SetBranchAddress("y0",&y0_out);
  tree_bcout->SetBranchAddress("u0",&u0_out);
  tree_bcout->SetBranchAddress("v0",&v0_out);

  //BcIn
  int ntrack_in;
  vector<double>* chisqr_in = nullptr;
  vector<double>* x0_in = nullptr;
  vector<double>* y0_in = nullptr;
  vector<double>* u0_in = nullptr;
  vector<double>* v0_in = nullptr;
  tree_bcin->SetBranchAddress("ntrack",&ntrack_in);
  tree_bcin->SetBranchAddress("chisqr",&chisqr_in);
  tree_bcin->SetBranchAddress("x0",&x0_in);
  tree_bcin->SetBranchAddress("y0",&y0_in);
  tree_bcin->SetBranchAddress("u0",&u0_in);
  tree_bcin->SetBranchAddress("v0",&v0_in);

  TH1D *hist_bh2[NumOfSegBH2];
  

  TH2D *hitpat_bcin_x = new TH2D("hitpat_bcin_x;BH2 Seg;X [mm]","hitpat_bcin_x",NumOfSegBH2,-0.5,NumOfSegBH2-0.5,500,-200,200);
  TH2D *hitpat_bcin_y = new TH2D("hitpat_bcin_y","hitpat_bcin_y",NumOfSegBH2,-0.5,NumOfSegBH2-0.5,500,-200,200);
  TH2D *hitpat_bcout_x = new TH2D("hitpat_bcout_x","hitpat_bcout_x",NumOfSegBH2,-0.5,NumOfSegBH2-0.5,500,-200,200);
  TH2D *hitpat_bcout_y = new TH2D("hitpat_bcout_y","hitpat_bcout_y",NumOfSegBH2,-0.5,NumOfSegBH2-0.5,500,-200,200);

  TH2D *hitpat_bht_bcin_x = new TH2D("hitpat_bht_bcin_x","hitpat_bht_bcin_x",NumOfSegBHT,-0.5,NumOfSegBHT-0.5,500,-200,200);
  TH2D *hitpat_bht_bcin_y = new TH2D("hitpat_bht_bcin_y","hitpat_bht_bcin_y",NumOfSegBHT,-0.5,NumOfSegBHT-0.5,500,-200,200);
  TH2D *hitpat_bht_bcout_x = new TH2D("hitpat_bht_bcout_x","hitpat_bht_bcout_x",NumOfSegBHT,-0.5,NumOfSegBHT-0.5,500,-200,200);
  TH2D *hitpat_bht_bcout_y = new TH2D("hitpat_bht_bcout_y","hitpat_bht_bcout_y",NumOfSegBHT,-0.5,NumOfSegBHT-0.5,500,-200,200);

  
  for(int i=0;i<NumOfSegBH2;i++){
    hist_bh2[i] = new TH1D(Form("hist_bh2%d",i),Form("hist_bh2%d",i),(tdc_max - tdc_min)/10.,tdc_min,tdc_max);
  }


  //for(int n=0;n<tree_hodo->GetEntries();n++){
  for(int n=0;n<100000;n++){
    tree_hodo->GetEntry(n);
    tree_bcin->GetEntry(n);
    tree_bcout->GetEntry(n);
    
    if(bh2_raw_seg->empty())continue;
    if(bht_raw_seg->empty())continue;
    
    for(int i=0;i<bh2_raw_seg->size();i++){
      bool bh2_pass = false;
      int bh2_seg = (*bh2_raw_seg)[i];
      for(int j=0;j<(*bh2_tdc_s)[i].size();j++){
	hist_bh2[bh2_seg]->Fill((*bh2_tdc_s)[i][j]);
	if((*bh2_tdc_s)[i][j]>tdc_window_min && (*bh2_tdc_s)[i][j]<tdc_window_max)bh2_pass = true;
      }
      if(bh2_pass){
	if(ntrack_in > 0){
	  hitpat_bcin_x->Fill(bh2_seg,(*x0_in)[0]);
	  hitpat_bcin_y->Fill(bh2_seg,(*y0_in)[0]);
	}
	if(ntrack_out > 0){
	  hitpat_bcout_x->Fill(bh2_seg,(*x0_out)[0]);
	  hitpat_bcout_y->Fill(bh2_seg,(*y0_out)[0]);
	}
      }
    }
    for(int i=0;i<bht_raw_seg->size();i++){
      int bht_seg = (*bht_raw_seg)[i];
      bool tdc_u = false;
      bool tdc_d = false;
      bool tdc_pass = false;
      for(int j=0;j<(*bht_tdc_u)[i].size();j++){
	if((*bht_tdc_u)[i][j] > tdc_bht_window_min && (*bht_tdc_u)[i][j] < tdc_bht_window_max){
	  tdc_u = true;
	  break;
	}
      }
	
      for(int j=0;j<(*bht_tdc_d)[i].size();j++){
	if((*bht_tdc_d)[i][j] > tdc_bht_window_min && (*bht_tdc_d)[i][j] < tdc_bht_window_max){
	  tdc_d = true;
	  break;
	}
      }
      if(bht_seg == 0 || bht_seg == 2 || bht_seg == 60 || bht_seg == 62)continue;
      if(tdc_u && tdc_d){
	if(ntrack_in>0){
	  hitpat_bht_bcin_x->Fill(bht_seg,(*x0_in)[0]);
	  hitpat_bht_bcin_y->Fill(bht_seg,(*y0_in)[0]);
	}
	if(ntrack_out>0){
	  hitpat_bht_bcout_x->Fill(bht_seg,(*x0_out)[0]);
	  hitpat_bht_bcout_y->Fill(bht_seg,(*y0_out)[0]);
	}
      }
    }
    
  }

  
  TCanvas* c = new TCanvas("c", "c", 900, 700);
  c->Print("bh2_blc_260602.pdf[");
  /*
  c->Divide(4,4);
  cout<<"1"<<endl;
  for(int i=0;i<NumOfSegBH2;i++){
    c->cd(i+1);
    hist_bh2[i]->Draw("HIST");
  }
  c->Print("bh2_blc.pdf");
  
  */
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.12);
  c->Divide(2);
  
  
  c->cd(1);
  hitpat_bcin_x->Draw("colz");
  c->cd(2);
  hitpat_bcin_y->Draw("colz");
  c->Print("bh2_blc_260602.pdf");
  c->Clear();
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.12);
  
  c->Divide(2);
  c->cd(1);
  hitpat_bcout_x->Draw("colz");
  c->cd(2);
  hitpat_bcout_y->Draw("colz");
  c->Print("bh2_blc_260602.pdf");

  c->Clear();
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.12);
  
  c->Divide(2);
  c->cd(1);
  hitpat_bht_bcin_x->Draw("colz");
  c->cd(2);
  hitpat_bht_bcin_y->Draw("colz");
  c->Print("bh2_blc_260602.pdf");

  c->Clear();
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetTopMargin(0.05);
  c->SetBottomMargin(0.12);
  
  c->Divide(2);
  c->cd(1);
  hitpat_bht_bcout_x->Draw("colz");
  c->cd(2);
  hitpat_bht_bcout_y->Draw("colz");
  c->Print("bh2_blc_260602.pdf");
  
  c->Print("bh2_blc_260602.pdf]");
  
    
  
}
