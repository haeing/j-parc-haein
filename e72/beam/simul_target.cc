static const double mK = 493.677;
static const double mp = 938.27208816;
static const double mn = 939.5654205;
static const double mpi = 139.57039;
static const double mpi0 = 134.9768;
static const double me = 0.51099895;

static const double meta =547.862;
static const double mLambda = 1115.683;
static const double BAC_thre = 1/1.115;
static const double KVC_thre = 1/1.45;

double convert_roots(double M1, double P1, double M2, double P2){
  double E1 = TMath::Sqrt(M1*M1 + P1*P1);
  double E2 = TMath::Sqrt(M2*M2 + P2*P2);
  double totalE = E1 + E2;
  double totalpz = P1 + P2;

  return TMath::Sqrt(totalE*totalE - totalpz*totalpz);
}

double GetBeta(int pdg, double mom){
  double mass = -9999;
  if(TMath::Abs(pdg) == 321)mass = mK;
  else if(TMath::Abs(pdg) == 2212)mass = mp;
  else if(TMath::Abs(pdg) == 2112)mass = mn;
  else if(TMath::Abs(pdg) == 211)mass = mpi;
  else if(TMath::Abs(pdg) == 111)mass = mpi0;
  else if(TMath::Abs(pdg) == 11)mass = me;

  if(mass == -9999)return -9999;
  else{ return mom/TMath::Sqrt(mass*mass + mom*mom);}
  
}

void simul_target(int mom){
  TFile* file = new TFile(Form("/gpfs/group/had/sks/Users/haein/E72_Simul/E72_Beam_Simul_%d.root",mom));
  TTree* tree = (TTree*)file->Get("g4hyptpc");

  vector<TParticle>* PRM = nullptr;
  vector<TParticle>* BH2 = nullptr;
  vector<TParticle>* BAC = nullptr;
  vector<TParticle>* TGT = nullptr;
  vector<TParticle>* HTOF = nullptr;
  vector<TParticle>* KVC = nullptr;

  TBranch *b_PRM = tree->GetBranch("PRM");
  b_PRM->SetAddress(&PRM);
  
  TBranch *b_BH2 = tree->GetBranch("BH2");
  b_BH2->SetAddress(&BH2);

  TBranch *b_BAC = tree->GetBranch("BAC");
  b_BAC->SetAddress(&BAC);

  TBranch *b_TGT = tree->GetBranch("TGT");
  b_TGT->SetAddress(&TGT);

  bool BAC_signal = false;
  bool TGT_signal = false;
  int BH2_pdg = -9999;
  int BAC_pdg = -9999;
  int TGT_pdg = -9999;
  double BH2_p_beam = -9999;
  double TGT_p_beam = -9999;
  double BH2_roots = -9999;
  double TGT_roots = -9999;
  bool BH2_trig = false;
  double Gen_p_beam = -9999;
  double Gen_roots = -9999;

  double BAC_beta;
  double KVC_beta;

  int TGT_accept = 0;

  TH2D *hist_target = new TH2D("hist_target","hist_target",300,-60,60,30,-60,60);
  
  for(int n=0;n<tree->GetEntries();n++){
    tree->GetEntry(n);

    BAC_signal = false;
    BAC_pdg = -9999;
    TGT_pdg = -9999;
    TGT_signal = false;

    if(BAC && !BAC->empty()){
      const auto&p_BAC = (*BAC)[0];
      BAC_pdg = p_BAC.GetPdgCode();
      double BAC_mom = TMath::Sqrt(p_BAC.Px()*p_BAC.Px()+p_BAC.Py()*p_BAC.Py()+p_BAC.Pz()*p_BAC.Pz());

      BAC_beta = GetBeta(BAC_pdg,BAC_mom);
      
      if(BAC_beta > BAC_thre)BAC_signal = true;
    }
    
    if(TGT && !TGT->empty()){
      const auto&p_TGT = (*TGT)[0];
      TGT_p_beam = TMath::Sqrt(p_TGT.Px()*p_TGT.Px()+p_TGT.Py()*p_TGT.Py()+p_TGT.Pz()*p_TGT.Pz());
      TGT_roots = convert_roots(mp,0,mK,TGT_p_beam);
      TGT_pdg = p_TGT.GetPdgCode();
      if(TGT_pdg == -321){
	TGT_signal = true;
	hist_target->Fill(p_TGT.Vx(),p_TGT.Vy());
      }
      
    }
    if(!BAC_signal && TGT_signal){
      TGT_accept++;
    }
    
    
  }

  std::cout<<(double)TGT_accept / (double)tree->GetEntries()<<std::endl;
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gStyle->SetOptStat(0);
  hist_target->GetXaxis()->SetTitle("X [mm]");
  hist_target->GetYaxis()->SetTitle("Y [mm]");
  hist_target->Draw("colz");

  c1->SaveAs("beam_target.pdf");
}
