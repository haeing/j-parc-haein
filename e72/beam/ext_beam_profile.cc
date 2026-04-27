#include <iostream>
#include <fstream>

const double target_pos = -143 - 150;
const double start_pos = -900.; //Geant4 coordinate is different!! This position is E72 experimental coordinate 

bool pion = true;

void ext_beam_profile(int mom){

  double bac_tdc_min = 0.;
  double bac_tdc_max = 0.;
  int runnumber = 0;
  double btof_cut = 0;
  if(mom == 645){
    bac_tdc_min = 740000;
    bac_tdc_max = 765000;
    runnumber = 2978;
    btof_cut = -4;
  }

  else if(mom == 665){
    bac_tdc_min = 740000;
    bac_tdc_max = 765000;
    runnumber = 3006;
    btof_cut = -3;
  }
  
  else if(mom == 685){
    bac_tdc_min = 730000;
    bac_tdc_max = 770000;
    runnumber = 2883;
    btof_cut = -3;
  }

  else if(mom == 715){
    bac_tdc_min = 740000;
    bac_tdc_max = 765000;
    runnumber = 2629;
    btof_cut = -3;
  }

  else if(mom == 735){
    bac_tdc_min = 700000;
    bac_tdc_max = 720000;
    runnumber = 2056;
    btof_cut = -3;
  }

  else if(mom == 755){
    bac_tdc_min = 730000;
    bac_tdc_max = 765000;
    runnumber = 2723;
    btof_cut = -2;
  }

  else if(mom == 790){
    bac_tdc_min = 735000;
    bac_tdc_max = 765000;
    runnumber = 2800;
    btof_cut = -2;
  }

  else if(mom == 814){
    bac_tdc_min = 740000;
    bac_tdc_max = 765000;
    runnumber = 2912;
    btof_cut = -2;
  }
  
  else if(mom == 842){
    bac_tdc_min = 740000;
    bac_tdc_max = 760000;
    runnumber = 3031;
    btof_cut = -1;
  }

  else if(mom == 870){
    bac_tdc_min = 740000;
    bac_tdc_max = 765000;
    runnumber = 2992;
    btof_cut = -1;
  }
  
  //beam trigger
  else if(mom == 933){
    bac_tdc_min = 0.;
    bac_tdc_max = 10.;
    runnumber = 2856;
    btof_cut = -1;
  }
  
  
  TString dir = "/gpfs/group/had/sks/Users/haein/JPARC2025Nov_root";
  TFile *file_track = new TFile(Form("%s/run0%d_BcOutTracking.root",dir.Data(),runnumber)); //Beam
  TTree *tree_track = (TTree*)file_track->Get("bcout");
  
  TFile *file_hodo = new TFile(Form("%s/run0%d_Hodoscope.root",dir.Data(),runnumber)); //Beam
  TTree *tree_hodo = (TTree*)file_hodo->Get("hodo");

  TFile *file_beam_old = new TFile("/home/had/haein/Work/E72/Simul/k18geant4/hyptpc-11.0.2/param/BEAM/beam.k.run69_0130.root");
  TTree *tree_beam_old = (TTree*)file_beam_old->Get("tr");
  
  TFile *file_beam = new TFile(Form("/hsm/had/sks/E72/JPARC2025Nov/beam_simul/beam_profile_run0%d_-%d_pi.root",runnumber,mom), "RECREATE");
  //TFile *file_beam = new TFile(Form("beam_profile_run0%d_-%d.root",runnumber,mom), "RECREATE");
  TTree * tree_beam = new TTree("tr","Beam Profile");


  double px,py,pz;
  
  double pointInx,pointIny,pointInz,pInx,pIny,pInz;

  tree_beam_old->SetBranchAddress("pInx",&px);
  tree_beam_old->SetBranchAddress("pIny",&py);
  tree_beam_old->SetBranchAddress("pInz",&pz);


  tree_beam->Branch("pointInx",&pointInx,"pointInx/D");
  tree_beam->Branch("pointIny",&pointIny,"pointIny/D");
  tree_beam->Branch("pointInz",&pointInz,"pointInz/D");
  tree_beam->Branch("pInx",&pInx,"pInx/D");
  tree_beam->Branch("pIny",&pIny,"pIny/D");
  tree_beam->Branch("pInz",&pInz,"pInz/D");
  

  vector<double> *x0 = nullptr;
  vector<double> *y0 = nullptr;

  vector<double> *u0 = nullptr;
  vector<double> *v0 = nullptr;

  TBranch *b_x0 = tree_track->GetBranch("x0");
  b_x0->SetAddress(&x0);
  TBranch *b_y0 = tree_track->GetBranch("y0");
  b_y0->SetAddress(&y0);
  TBranch *b_u0 = tree_track->GetBranch("u0");
  b_u0->SetAddress(&u0);
  TBranch *b_v0 = tree_track->GetBranch("v0");
  b_v0->SetAddress(&v0);

  vector<vector<double>> *bac_tdc_u = nullptr;
  bac_tdc_u = new vector<vector<double>>();

  TBranch *b_bac_tdc_u = tree_hodo->GetBranch("bac_tdc_u");
  b_bac_tdc_u->SetAddress(&bac_tdc_u);

  double btof0;
  tree_hodo->SetBranchAddress("btof0",&btof0);

  bool btof_kaon = false;
  bool bac_pion = false;

  TH1D *hist_bac_tdc = new TH1D("hist_bac_tdc","hist_bac_tdc",1000,600000,800000);
  TH1D *hist_btof0 = new TH1D("hist_btof0","hist_btof0",400,-20,20);
  TH1D *hist_mom = new TH1D("hist_mom","hist_mom",280,mom-70,mom+70);
  


  int total_old = tree_beam_old->GetEntries();
  int total = tree_track->GetEntries();
  
  //for(int i=0;i<(total_old > total)?total : total_old;i++){
  //for(int i=0;i<total_old;i++){
  for(int i=0;i<20000;i++){
    tree_track->GetEntry(i);
    tree_hodo->GetEntry(i);
    tree_beam_old->GetEntry(i);
    
    btof_kaon = false;
    bac_pion = false;
    
    double pbeam = sqrt(px*px + py*py + pz*pz) * (double)mom / 907;

    hist_btof0->Fill(btof0);
    if(btof0<btof_cut)btof_kaon = true;
    if(!bac_tdc_u->empty()){
      if((*bac_tdc_u)[4].size()>0){
	for(int j=0;j<(*bac_tdc_u)[4].size();j++){
	  hist_bac_tdc->Fill((*bac_tdc_u)[4][j]);
	  if((*bac_tdc_u)[4][j]>bac_tdc_min && (*bac_tdc_u)[4][j]<bac_tdc_max)
	    bac_pion = true;
	}
      }
    }
    if(!x0->empty() && !y0->empty() && !u0->empty() && !v0->empty()){
      if(!pion){
	if(btof_kaon && !bac_pion){
	  pointInx = (*x0)[0]+start_pos*(*u0)[0];
	  pointIny = (*y0)[0]+start_pos*(*v0)[0];
	  pointInz = start_pos + 150.071;
	  
	  double norm = sqrt(1+(*u0)[0]*(*u0)[0] + (*v0)[0]*(*v0)[0]);
	  pInx = pbeam*(*u0)[0]/norm;
	  pIny = pbeam*(*v0)[0]/norm;
	  pInz = pbeam/norm;
	  hist_mom->Fill(sqrt(pInx*pInx+pIny*pIny+pInz*pInz)*1000);
	  
	  tree_beam->Fill();
	}
	
      }
      else if(pion){
	if(!btof_kaon){
	  pointInx = (*x0)[0]+start_pos*(*u0)[0];
	  pointIny = (*y0)[0]+start_pos*(*v0)[0];
	  pointInz = start_pos + 150.071;
	  
	  double norm = sqrt(1+(*u0)[0]*(*u0)[0] + (*v0)[0]*(*v0)[0]);
	  pInx = pbeam*(*u0)[0]/norm;
	  pIny = pbeam*(*v0)[0]/norm;
	  pInz = pbeam/norm;
	  hist_mom->Fill(sqrt(pInx*pInx+pIny*pIny+pInz*pInz)*1000);
	  
	  tree_beam->Fill();
	}
      }
    }
    
  }
  


  
  std::cout<<"write close"<<std::endl;

  std::ofstream fout(Form("E72_%d.csv",mom));


  fout << "bin,center,value,error\n";

  int nbins = hist_mom->GetNbinsX();
  for (int i = 1; i <= nbins; ++i) {
    int    bin     = i;
    double center  = hist_mom->GetBinCenter(i);
    double content = hist_mom->GetBinContent(i);
    double error   = hist_mom->GetBinError(i);
    
    fout << bin << ","
	 << center << ","
	 << content << ","
	 << error << "\n";
  }

  
  fout.close();
  hist_bac_tdc->Write();
  hist_btof0->Write();
  hist_mom->Write();
  tree_beam->Write();
  file_beam->Close();
 
}
