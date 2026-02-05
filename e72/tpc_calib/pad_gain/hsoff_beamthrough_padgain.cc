//const int runnumber[7] = {2570, 2581, 2585, 2587, 2589, 2590, 2592};
const int runnumber[6] = {2581, 2585, 2587, 2589, 2590, 2592};

void hsoff_beamthrough_padgain(){

  for(int run : runnumber){
    string dir = "/gpfs/group/had/sks/Users/haein/JPARC2025Nov_root/hsoff_beamthrough";
    TFile *file_hodo = new TFile(Form("%s/run0%d_Hodoscope.root",dir.c_str(),run));
    TFile *file_bcout = new TFile(Form("%s/run0%d_BcOutTracking.root",dir.c_str(),run));
    TFile *file_tpc = new TFile(Form("%s/run0%d_DstTPCTracking.root",dir.c_str(),run));

    TTree *tree_hodo = (TTree*)file_hodo->Get("hodo");
    TTree *tree_bcout = (TTree*)file_bcout->Get("bcout");
    TTree *tree_tpc = (TTree*)file_tpc->Get("tpc");

    for(int n=0;n<tree_hodo->GetEntries();n++){
      tree_hodo->GetEntry(n);
      tree_bcout->GetEntry(n);
      tree_tpc->GetEntry(n);

      //Hodo
      double btof0;
      tree_hodo->SetBranchAddress("btof0",&btof0);

      //BcOut
      vector<double>* x0 = nullptr;
      vector<double>* y0 = nullptr;
      vector<double>* u0 = nullptr;
      vector<double>* v0 = nullptr;
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
      
    }
  }
    
}
