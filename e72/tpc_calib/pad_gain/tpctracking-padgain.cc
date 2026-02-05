const Int_t NumOfPadTPC       = 5768;
const Int_t NumOfLayersTPC    = 32;

const Double_t MinX = -300;
const Double_t MaxX = 300;
const Double_t MinZ = -300;
const Double_t MaxZ = 300;

const Double_t ZTarget = -143.;
//_____________________________________________________________________________
//#OfPad #division #radius padLength
static const Double_t padParameter[NumOfLayersTPC][6] =
{{0, 48,    14.75, 48, 0,  9.},
 {1, 48,    24.25, 48, 0,  9.},
 {2, 72,    33.75, 72, 0,  9.},
 {3, 96,    43.25, 96, 0,  9.},
 {4, 120,    52.75,120,0,   9.},
 {5, 144,    62.25,144,0,   9.},
 {6, 168,    71.75,168,0,   9.},
 {7, 192,    81.25,192,0,   9.},
 {8, 216,    90.75,216,0,   9.},
 {9, 240,    100.25,240,0,  9.},
 {10,208,    111.5,241, 0,  12.5},
 {11,218,    124.5,271, 0,  12.5},
 {12,230,    137.5,300, 0,  12.5},
 {13,214,    150.5,330, 0,  12.5},
 {14,212,    163.5,360, 0,  12.5},
 {15,214,    176.5,390, 0,  12.5},
 {16,220,    189.5,420, 0,  12.5},
 {17,224,    202.5,449, 0,  12.5},
 {18,232,    215.5,479, 0,  12.5},
 {19,238,    228.5,509, 0,  12.5},
 {20,244,    241.5,539, 0,  12.5},
 {21,232,    254.5,569, 0,  12.5},
 {22,218,    267.5,599, 0,  12.5},
 {23,210,    280.5,628, 0,  12.5},
 {24,206,    293.5,658, 0,  12.5},
 {25,202,    306.5,688, 0,  12.5},
 {26,200,    319.5,718, 0,  12.5},
 {27,196,    332.5,748, 0,  12.5},
 {28,178,    345.5,777, 0,  12.5},
 {29,130,    358.5,807, 0,  12.5},
 {30,108,    371.5,837, 0,  12.5},
 {31,90,     384.5,867, 0, 12.5}};

//_____________________________________________________________________________
Int_t GetPadId(Int_t layer_id, Int_t row_id)
{
  // Return -1 (invalid) for negative row_id to prevent crashes.
  // This is primarily to handle the INT_MIN(-2147483648) value resulting from (Int_t)NaN
  // when initializing TPCHit with NaN in TPCCluster as below;
  //  m_mean_hit(new TPCHit(layer, TMath::QuietNaN())) <- this
  if (row_id < 0) return -1;

  Int_t pad_id = 0;
  for (Int_t layer = 0; layer < layer_id; layer++)
    pad_id += static_cast<Int_t>(padParameter[layer][1]);
  pad_id += row_id;
  return pad_id;
}

//_____________________________________________________________________________
static const Int_t padOnSectionFrame[] =
{
//Gem section 1 frame 0
1503 ,1504, 1505 ,1506 ,1507 ,1508 ,1509 ,1510 ,1511 ,1512 ,1513 ,1514 ,1515 ,1728 ,1729 ,1730,

//Gem section 1 frame 1
72 ,134 ,135 ,222 ,332 , 408,  467 ,626 ,809 ,1016, 59 ,110 ,185 ,284,555 ,726, 921 ,1140 ,1363, 1564 ,1565 ,1778,

//Gem section 2 frame 0, frame 1 + adjacent pads
1628,1629,1630,1817,1818,1819,1820,1849,1850,1851,2035,2036,2037,2038,2067,2068,2069,2244,2245,2246,2247,2248,2276,2277,2278,2454,2455,2456,2457,2485,2486,2487,2488,2667,2668,2669,2670,2698,2699,2700,2701,2886,2887,2888,2917,2918,2919,3110,3111,3112,3113,3141,3142,3143,3341,3342,3343,3344,3373,3374,3375,3578,3579,3580,3581,3610,3611,3612,3813,3814,3815,3844,3845,3846,4034,4035,4036,4065,4066,4067,4276,4277,4278,4480,4481,4482,4680,4681,4682,4878,4879,4880,5072,5073,5074,

//Gem section 3 frame 0
4720 ,4721 ,4722 ,4723 ,4724 ,4725 ,4726 ,4727 ,4728 ,4929 ,4930 ,4931 ,4932 ,4933 ,4934 ,4935 ,4936 ,5134 ,5135 ,5136 ,5137 ,5138 ,5139 ,5140 ,5328 ,5329 ,5330 ,5331 ,5332 ,5333 ,5487 ,5488 ,5489 ,5490 ,5491 ,5492 ,5612 ,5613 ,5614 ,5615 ,5616 ,5716 ,5717 ,5718 ,5719 ,5720, 5532 ,5533, 5612 ,5613 ,5614 ,5615 ,5616 ,5654 ,5655 ,5716 ,5717 ,5718 ,5719 ,5720 ,5757 ,5758,

//Gem section 3 frame 1
3413 ,3414 ,3415 ,3416 ,3417 ,3418 ,3419 ,3660 ,3661 ,3662 ,3663 ,3664 ,3665 ,3904 ,3905 ,3906 ,3907 ,3908 ,4134 ,4135 ,4136 ,4137 ,4353 ,4354 ,4355 ,4356 ,4565 ,4566 ,4567 ,4568 ,4774 ,4775 ,4776 ,4981 ,4982 ,5183 ,5184 ,5372 ,5373 ,5374 ,5530 ,5531 ,5532 ,5533 ,5654 ,5655 ,5757 ,5758,

//Gem section 4 frame 0
4409 ,4410 ,4411 ,4412 ,4413 ,4414 ,4415 ,4416 ,4417 ,4418 ,4419 ,4420 ,4421 ,4422 ,4423 ,4424 ,4425 ,4426 ,4427 ,4428 ,4429 ,4430 ,4431 ,4432 ,4433 ,4434 ,4435 ,4436 ,4437 ,4438 ,4439 ,4440 ,4441 ,4442 ,4443 ,4444 ,4445 ,4446 ,4447 ,4448 ,4449 ,4450 ,4451 ,4452 ,4611 ,4612 ,4613 ,4614 ,4615 ,4616 ,4617 ,4618 ,4619,

//Gem section 4 frame 1
2794 ,2795 ,2796 ,2797 ,2798 ,2799 ,2800 ,2801 ,2802 ,2803 ,2804 ,2805 ,2806 ,2807 ,2808 ,2809 ,3001 ,3002 ,3003 ,3004 ,3005 ,3006 ,3007 ,3008 ,3009 ,3010 ,3011 ,3012 ,3013 ,3014 ,3015 ,3016 ,3017, 3018, 3037, 3038, 3039 ,3040 ,3041 ,3042 ,3043 ,3044 ,3045 ,3046 ,3047 ,3048 ,3049 ,3050 ,3051 ,3052 ,3053 ,3054 ,3227 ,3228 ,3229 ,3230 ,3231 ,3289 ,3290 ,3291 ,3292 ,3293 ,3294 ,3295 ,3296 ,3297 ,3540 ,3541 ,3542 ,3543 ,3544 ,3545 ,3546
};

//_____________________________________________________________________________
inline Bool_t Noise(Int_t pad_id){
  Bool_t noise = std::find(std::begin(padOnSectionFrame), std::end(padOnSectionFrame), pad_id) != std::end(padOnSectionFrame);
  if(noise) return true;
  else return false;
}


void tpctracking_padgain(){
  string dir = "/gpfs/group/had/sks/Users/haein/JPARC2025Nov_root";
  TFile *file = new TFile(Form("%s/run02489_TPCTracking.root",dir.c_str()));
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
