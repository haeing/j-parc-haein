#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include "TGeoManager.h"
#include "TGeoXtru.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"

#include "../include/TPCPadHelper.hh"

const Int_t MaxTPCHits = 50000;
const Int_t MaxHits = 50000;
const Int_t MaxTPCTracks = 5000;

Double_t Extrapolate(Double_t radius, TVector3 center, Double_t slope, Double_t dip, Double_t alpha, TVector3 &ptOnHelix);

enum CommandIdentifiers
  {
   M_FILE_OPEN,
   M_FILE_SAVE,
   M_FILE_SAVEAS,
   M_FILE_PRINT,
   M_FILE_PRINTSETUP,
   M_FILE_EXIT,

   M_ANA_RECO,
   M_ANA_FIT,
   
   M_VIEW_ENBL_DOCK,
   M_VIEW_ENBL_HIDE,
   M_VIEW_DOCK,
   M_VIEW_UNDOCK,

   M_HELP_CONTENTS,
   M_HELP_SEARCH,
   M_HELP_ABOUT,

  };

const char *filetypes[] = { "All files",     "*",
                            "ROOT files",    "*.root",
                            "ROOT macros",   "*.C",
                            "Text files",    "*.[tT][xX][tT]",
                            0,               0 };

struct anaData_t {
   const char *pixmap_name;
   const char *tip_text;
   Int_t       id;
   TGButton   *button;
};

anaData_t histo_data[] = {
   { "h1_s.xpm",        "TH1",      1001,  0 },
   { "h2_s.xpm",        "TH2",      1002,  0 },
   { "h3_s.xpm",        "TH3",      1003,  0 },
   { "profile_s.xpm",   "TProfile", 1004,  0 },
   { 0,                 0,          0,     0 }
};

anaData_t function_data[] = {
   { "f1_s.xpm",        "TF1",      2001,  0 },
   { "f2_s.xpm",        "TF2",      2002,  0 },
   { 0,                 0,          0,     0 }
};

anaData_t tree_data[] = {
   { "ntuple_s.xpm",    "TNtuple",  3001,  0 },
   { "tree_s.xpm",      "TTree",    3002,  0 },
   { "chain_s.xpm",     "TChain",   3003,  0 },
   { 0,                 0,          0,     0 }
};

struct Event {

  Int_t evnum;
  Int_t status;
  Int_t nhittpc;                 // Number of Hits                                                                              
  Int_t max_ititpc;

  Int_t nhPrm;
  Double_t xPrm[MaxTPCHits];
  Double_t yPrm[MaxTPCHits];
  Double_t zPrm[MaxTPCHits];
  Double_t pxPrm[MaxTPCHits];
  Double_t pyPrm[MaxTPCHits];
  Double_t pzPrm[MaxTPCHits];
  Double_t thetaPrm[MaxTPCHits];
  Double_t phiPrm[MaxTPCHits];
  
  Int_t ititpc[MaxTPCHits];
  Int_t nhittpc_iti[MaxTPCHits];

  Double_t xtpc[MaxTPCHits];
  Double_t ytpc[MaxTPCHits];
  Double_t ztpc[MaxTPCHits];
  Double_t x0tpc[MaxTPCHits];
  Double_t y0tpc[MaxTPCHits];
  Double_t z0tpc[MaxTPCHits];

  Double_t pxtpc[MaxTPCHits];
  Double_t pytpc[MaxTPCHits];
  Double_t pztpc[MaxTPCHits];

  Double_t edeptpc[MaxTPCHits];
  Int_t laytpc[MaxTPCHits];
  Int_t rowtpc[MaxTPCHits];
  
  Int_t nhTpc;
  Int_t nh_cluster_Tpc;

  std::vector<Double_t> *raw_hitpos_x;
  std::vector<Double_t> *raw_hitpos_y;
  std::vector<Double_t> *raw_hitpos_z;
  std::vector<Double_t> *raw_de;
  std::vector<Int_t> *raw_layerid;
  std::vector<Int_t> *raw_rowid;
  std::vector<Int_t> *raw_padid;
  std::vector<Double_t> *raw_hitpatpos_x;
  std::vector<Double_t> *raw_hitpatpos_y;
  std::vector<Double_t> *raw_hitpatpos_z;


  Int_t nGFTracks;

  std::vector<Int_t> *GFStatus;
  std::vector<Int_t> *GenFitID;

  std::vector<Double_t>  *nhit; 

  std::vector<Double_t> *radius;
  std::vector<Double_t> *xcenter;
  std::vector<Double_t> *zcenter;
  std::vector<Double_t> *yoffset;
  std::vector<Double_t> *slope;
  std::vector<Double_t> *dipangle;

  std::vector<Double_t> *alphaHead;
  std::vector<Double_t> *alphaTail;

  std::vector<Double_t> *helixmom;
  std::vector<Double_t> *tracklength;

  std::vector<Double_t> *rmsW;
  std::vector<Double_t> *rmsH;

  std::vector<Double_t> *charge;

  std::vector<Int_t> *isBeam;

  std::vector<Int_t> *helixnclusters;
  std::vector< std::vector<Double_t> > *helixxcluster;
  std::vector< std::vector<Double_t> > *helixycluster;
  std::vector< std::vector<Double_t> > *helixzcluster;
  
  std::vector<Double_t> *GFChisqr;
  std::vector<Int_t> *GFCharge;
  std::vector<Double_t> *GFpValue;
  std::vector<Double_t> *gfPx;
  std::vector<Double_t> *gfPy;
  std::vector<Double_t> *gfPz;

  std::vector<Double_t> *vtxPx;
  std::vector<Double_t> *vtxPy;
  std::vector<Double_t> *vtxPz;

  std::vector<Double_t> *vtxPosx;
  std::vector<Double_t> *vtxPosy;
  std::vector<Double_t> *vtxPosz;
  
  std::vector<Double_t> *htofPx;
  std::vector<Double_t> *htofPy;
  std::vector<Double_t> *htofPz;

  std::vector<Double_t> *htofPosx;
  std::vector<Double_t> *htofPosy;
  std::vector<Double_t> *htofPosz;

  std::vector<Double_t> *tgtPx;
  std::vector<Double_t> *tgtPy;
  std::vector<Double_t> *tgtPz;

  std::vector<Double_t> *tgtPosx;
  std::vector<Double_t> *tgtPosy;
  std::vector<Double_t> *tgtPosz;

  std::vector<Double_t > *gfnclusters;
  std::vector<std::vector<Double_t>> *gfxcluster;
  std::vector<std::vector<Double_t>> *gfycluster;
  std::vector<std::vector<Double_t>> *gfzcluster;

  std::vector<Double_t> GFTrackLength;
  std::vector<Double_t> GFTrackLengthHTOF;
  std::vector<Double_t> GFTrackTOF;
  std::vector<Double_t> GFTrackTOFHTOF;

  std::vector<Double_t> vtxWeight;

  Int_t numVertices;
  std::vector<Int_t> VtxNumTrack;
  std::vector<Int_t> RaveVtxID;
  std::vector<Double_t> RaveChisqr;
  std::vector<Double_t> RaveVtxPosx;
  std::vector<Double_t> RaveVtxPosy;
  std::vector<Double_t> RaveVtxPosz;
  std::vector<Double_t> RaveVtxCovx;
  std::vector<Double_t> RaveVtxCovy;
  std::vector<Double_t> RaveVtxCovz;

  Int_t nKstar;

  std::vector<Double_t>* mKstar;
  std::vector<Double_t>* mKstar_init;
  std::vector<Double_t>* chisqrKstar;
  std::vector<Double_t>* fitprobKstar;

  std::vector<Double_t>* xKstar;
  std::vector<Double_t>* yKstar;
  std::vector<Double_t>* zKstar;

  std::vector<Double_t>* pxKstar;
  std::vector<Double_t>* pyKstar;
  std::vector<Double_t>* pzKstar;

  Int_t nKs0;

  std::vector<Double_t>* mKs0;
  std::vector<Double_t>* mKs0_init;
  std::vector<Double_t>* chisqrKs0;
  std::vector<Double_t>* fitprobKs0;

  std::vector<Double_t>* xKs0;
  std::vector<Double_t>* yKs0;
  std::vector<Double_t>* zKs0;

  std::vector<Double_t>* pxKs0;
  std::vector<Double_t>* pyKs0;
  std::vector<Double_t>* pzKs0;

};

class Analysis {

  RQ_OBJECT("Analysis")

  private:
  
  TGTransientFrame *fMain;
  TGCompositeFrame *fFrame1, *fFrame2, *fFrame3;
  TGButton *fOkButton, *fCancelButton, *fStartButton, *fStopButton;
  TGLayoutHints *fL1, *fL2, *fL3, *fL4;
  TRootEmbeddedCanvas *fECanvas1, *fECanvas2;
  Bool_t fFillHistos;
  TH1D  *fHpx1;
  TH1D  *fHpx2;

public:

  Analysis(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h);
  virtual ~Analysis();

  void DoOK();
  void DoCancel();
  void DoClose();
  void CloseWindow();
};

class MyMainFrame {
   RQ_OBJECT("MyMainFrame")

public:
  MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
  virtual ~MyMainFrame();

  void HandleMenu(Int_t id);

  Bool_t InitializeTPC();
  TGeoVolume *TPCGeometry();
  Bool_t loadFile();
  Bool_t loadEvent();

  void GoForward();
  void GoBackward();
  void GoTo();
  
  void CloseWindow();
  
private:
  TGMainFrame         *fMain;
  TGDockableFrame     *fMenuDock;
  TRootEmbeddedCanvas *fEcanvas;
  TCanvas *fCanvas;
  TGeoVolume *fTPC;
  TH2Poly             *fTPC2dPoly;
  TH2Poly             *fHTOF2dPoly;
  TH2D             *fTPCxz;
  TH2D             *fTPCxz_raw;
  TH2D             *fTPCxz_helix;
  TH2D             *fTPCxz_genfit;
  TH2D             *fTPCKstar;
  TH2D             *fTPCKs0;
  TGTextEntry      *fStatusBar;
  TGTextEntry      *fEvtHandler;

  TGMenuBar          *fMenuBar;
  TGPopupMenu        *fMenuFile, *fMenuAnalysis, *fMenuView, *fMenuHelp;
  TGLayoutHints      *fMenuBarLayout, *fMenuBarItemLayout, *fMenuBarHelpLayout;

};

MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h) {
   // Create a main frame
   fMain = new TGMainFrame(p,w,h);

   fMain->SetCleanup(kDeepCleanup);
   fMain->Connect("CloseWindow()","MyMainFrame",this,"CloseWindow()");
   // Create canvas widget

   fMenuDock = new TGDockableFrame(fMain);
   fMain->AddFrame(fMenuDock, new TGLayoutHints(kLHintsExpandX, 0, 0, 1, 0));
   fMenuDock->SetWindowName("Menu");

   fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsExpandX);
   fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
   fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

   fMenuFile = new TGPopupMenu(gClient->GetRoot());
   fMenuFile->AddEntry("&Open...", M_FILE_OPEN);
   fMenuFile->AddEntry("&Save", M_FILE_SAVE);
   fMenuFile->AddEntry("S&ave as...", M_FILE_SAVEAS);
   fMenuFile->AddEntry("&Close", -1);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("E&xit", M_FILE_EXIT);

   fMenuAnalysis = new TGPopupMenu(gClient->GetRoot());
   fMenuAnalysis->AddLabel("Menu for Analysis");
   fMenuAnalysis->AddSeparator();
   fMenuAnalysis->AddEntry("&Reconstruction",M_ANA_RECO);
   fMenuAnalysis->AddEntry("&Fit",M_ANA_FIT);
      
   fMenuView = new TGPopupMenu(gClient->GetRoot());
   fMenuView->AddEntry("&Dock", M_VIEW_DOCK);
   fMenuView->AddEntry("&Undock", M_VIEW_UNDOCK);
   fMenuView->AddSeparator();
   fMenuView->AddEntry("Enable U&ndock", M_VIEW_ENBL_DOCK);
   fMenuView->AddEntry("Enable &Hide", M_VIEW_ENBL_HIDE);
   fMenuView->DisableEntry(M_VIEW_DOCK);
      
   fMenuHelp = new TGPopupMenu(gClient->GetRoot());
   fMenuHelp->AddEntry("&Contents", M_HELP_CONTENTS);
   fMenuHelp->AddEntry("&Search...", M_HELP_SEARCH);
   fMenuHelp->AddSeparator();
   fMenuHelp->AddEntry("&About", M_HELP_ABOUT);

   fMenuFile->Connect("Activated(Int_t)", "MyMainFrame", this,
                      "HandleMenu(Int_t)");
   fMenuAnalysis->Connect("Activated(Int_t)", "MyMainFrame", this,
                      "HandleMenu(Int_t)");
   fMenuView->Connect("Activated(Int_t)", "MyMainFrame", this,
                      "HandleMenu(Int_t)");
   fMenuHelp->Connect("Activated(Int_t)", "MyMainFrame", this,
                      "HandleMenu(Int_t)");
   
   fMenuBar = new TGMenuBar(fMenuDock, 1, 1, kHorizontalFrame);
   fMenuBar->AddPopup("&File", fMenuFile, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Analysis", fMenuAnalysis, fMenuBarItemLayout);
   fMenuBar->AddPopup("&View", fMenuView, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Help", fMenuHelp, fMenuBarHelpLayout);
   
   fMenuDock->AddFrame(fMenuBar, fMenuBarLayout);
   

   fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,1600,800);

   fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX |
                   kLHintsExpandY, 10,10,10,1));
   
   fTPC2dPoly = new TH2Poly("fTPC2dPoly","fTPC2dPoly",-400,400,-400,400);
   fHTOF2dPoly = new TH2Poly("fHTOF2dPoly","fHTOF2dPoly",-400,400,-400,400);
   fTPCxz = new TH2D("fTPCxz","fTPCxz",300,-300,300,300,-300,300);
   fTPCxz_helix = new TH2D("fTPCxz_helix","fTPCxz_helix",300,-300,300,300,-300,300);
   fTPCxz_genfit = new TH2D("fTPCxz_genfit","fTPCxz_genfit",300,-300,300,300,-300,300);
   fTPCxz_raw = new TH2D("fTPCxz_raw","fTPCxz_raw",300,-300,300,300,-300,300);

   fTPCKstar = new TH2D("fTPCKstar","fTPCKstar",300,-300,300,300,-300,300);
   fTPCKs0 = new TH2D("fTPCKs0","fTPCKs0",300,-300,300,300,-300,300);

   fTPC = TPCGeometry();

   // Create a horizontal frame widget with buttons
   TGHorizontalFrame *hframe1 = new TGHorizontalFrame(fMain,1600,40);
   TString iconDir(TString::Format("%s/icons",gSystem->Getenv("ROOTSYS")));
   TGPictureButton *bwd = new TGPictureButton(hframe1, gClient->GetPicture(iconDir +"GoBack.gif"));
   bwd->Connect("Clicked()","MyMainFrame",this, "GoBackward()");
   hframe1->AddFrame(bwd,new TGLayoutHints(kLHintsCenterX | kLHintsExpandX,
					  2,2,2,2));

   TGPictureButton *fwd = new TGPictureButton(hframe1, gClient->GetPicture(iconDir +"GoFoward.gif"));
   fwd->Connect("Clicked()","MyMainFrame",this, "GoForward()");
   hframe1->AddFrame(fwd, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX,
					  2,2,2,2));

   TGTextButton *shortcut = new TGTextButton(hframe1);
   shortcut->Connect("Clicked()","MyMainFrame",this, "GoTo()");
   hframe1->AddFrame(shortcut, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
                                            2,2,2,2));

   fEvtHandler = new TGTextEntry(hframe1);
   fEvtHandler->SetEnabled(true);
   fEvtHandler->SetInsertMode();
   hframe1->AddFrame(fEvtHandler,new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsExpandY,
						     2,2,2,2));
   
   TGTextButton *exit = new TGTextButton(hframe1,"&Exit",
                                "gApplication->Terminate(0)");
   hframe1->AddFrame(exit, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
                                            2,2,2,2));
   TGHorizontalFrame *hframe2 = new TGHorizontalFrame(fMain,1600,40);

   fStatusBar = new TGTextEntry(hframe2);
   fStatusBar->SetEnabled(kFALSE);
   hframe2->AddFrame(fStatusBar,new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2)); 

   fMain->AddFrame(hframe1, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX,
                                             2,2,2,2));
   fMain->AddFrame(hframe2, new TGLayoutHints(kLHintsCenterX| kLHintsExpandX ,
                                             2,2,2,2));


   // Set a name to the main frame
   fMain->SetWindowName("Simple Example");

   // Map all subwindows of main frame
   fMain->MapSubwindows();

   // Initialize the layout algorithm
   fMain->Resize(fMain->GetDefaultSize());

   // Map main frame
   fMain->MapWindow();

   fCanvas = fEcanvas->GetCanvas();
   fCanvas->Divide(2);

   if (loadFile())
     loadEvent();
}

Bool_t MyMainFrame::InitializeTPC() {

   Double_t X[5]; Double_t Y[5];
   
   for (int l = 0 ; l < 32; l++) {
     Double_t pLength = tpc::padParameter[l][5];
     Double_t st = (180. - ( 360./tpc::padParameter[l][3]) * tpc::padParameter[l][1]/2.);
     Double_t sTheta = (-1 + st/180.) * TMath::Pi();
     Double_t dTheta = (360./tpc::padParameter[l][3])/180.*TMath::Pi();
     Double_t cRad = tpc::padParameter[l][2];
     Double_t nPad = tpc::padParameter[l][1];

     for (int j = 0 ; j < nPad ; j++) {

       Bool_t dead = tpc::Dead(l,j);
       if (!dead)  {
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

	 for (int k = 0 ; k < 5 ; k++)
	   X[k] += tpc::ZTarget;
       }

       fTPC2dPoly->AddBin(5,X,Y);
     }
   }
   
   fTPC2dPoly->SetStats(0);
   {
     const Double_t L = 337;
     const Double_t t = 10;
     const Double_t w = 68;
     Double_t theta[8];
     Double_t X[5];
     Double_t Y[5];
     Double_t seg_X[5];
     Double_t seg_Y[5];
     for( Int_t i=0; i<8; i++ ){
       theta[i] = (-180+45*i)*acos(-1)/180.;
       for( Int_t j=0; j<4; j++ ){
	 seg_X[1] = L-t/2.;
	 seg_X[2] = L+t/2.;
	 seg_X[3] = L+t/2.;
	 seg_X[4] = L-t/2.;
	 seg_X[0] = seg_X[4];
	 seg_Y[1] = w*j-2*w;
	 seg_Y[2] = w*j-2*w;
	 seg_Y[3] = w*j-1*w;
	 seg_Y[4] = w*j-1*w;
	 seg_Y[0] = seg_Y[4];
	 for( Int_t k=0; k<5; k++ ){
	   X[k] = cos(theta[i])*seg_X[k]-sin(theta[i])*seg_Y[k];
	   Y[k] = sin(theta[i])*seg_X[k]+cos(theta[i])*seg_Y[k];
	 }
	 fHTOF2dPoly->AddBin(5, X, Y);
       }
     }
   }
   

   return true;
}

MyMainFrame::~MyMainFrame() {
   // Clean up used widgets: frames, buttons, layout hints  
  fMain->Cleanup();
  delete fMain;
}

namespace {
  Event event;
  Int_t eventID;
  TFile *datafile;
  TTree *datatree;
  //
  const char *filename = "../../macro/kstar.root";
  //  const char *filename = "./geant4/proton_300MeV_EM_GenfitTPCTrackingRiemannGeant4.root";
  //  const char *filename = "./rootfile/run05764_DstTPCRiemannTracking.root";
  //  const char *filename = "~/Bureau/JPARC/data/E42/TPC/run05641_DstTPCKKAnaHS.root";
  //  const char *filename = "../../data/run05641_DstTPCKKAnaHS.root";
  //  const char *filename = "./rootfile/run05641_DstTPCRiemannKKAna.root";
  //
}
Bool_t MyMainFrame::loadFile() {

  printf("################### Opening File ####################\n");
  
  datafile = TFile::Open(filename,"READ");
  //  datatree = (TTree*)datafile->Get("tpc");
  datatree = (TTree*)datafile->Get("ktpc_g");
  
  if (!datafile) {
    Warning("TPCEventDisplay","No such file or directory");
    return false;
  }
  
  {
    datatree->SetBranchAddress("evnum",&event.evnum);
    /*    
    datatree->SetBranchAddress("nhTpc", &event.nhTpc);
    datatree->SetBranchAddress("nh_cluster_Tpc", &event.nh_cluster_Tpc);
    datatree->SetBranchAddress("raw_hitpos_x", &event.raw_hitpos_x);
    datatree->SetBranchAddress("raw_hitpos_y", &event.raw_hitpos_y);
    datatree->SetBranchAddress("raw_hitpos_z", &event.raw_hitpos_z);
    datatree->SetBranchAddress("raw_de", &event.raw_de);
    datatree->SetBranchAddress("raw_layerid", &event.raw_layerid);
    datatree->SetBranchAddress("raw_rowid", &event.raw_rowid);
    datatree->SetBranchAddress("raw_padid", &event.raw_padid);
    datatree->SetBranchAddress("raw_hitpatpos_x", &event.raw_hitpatpos_x);
    datatree->SetBranchAddress("raw_hitpatpos_y", &event.raw_hitpatpos_y);
    datatree->SetBranchAddress("raw_hitpatpos_z", &event.raw_hitpatpos_z);    
    */
    datatree->SetBranchAddress("nGFTracks",&event.nGFTracks);
    datatree->SetBranchAddress("GFStatus",&event.GFStatus);
    
    datatree->SetBranchAddress("radius",&event.radius);
    datatree->SetBranchAddress("xcenter",&event.xcenter);
    datatree->SetBranchAddress("zcenter",&event.zcenter);
    datatree->SetBranchAddress("yoffset",&event.yoffset);
    datatree->SetBranchAddress("slope",&event.slope);
    datatree->SetBranchAddress("dipangle",&event.dipangle);
    datatree->SetBranchAddress("alphaHead",&event.alphaHead);
    datatree->SetBranchAddress("alphaTail",&event.alphaTail);
    datatree->SetBranchAddress("helixmom",&event.helixmom);
    datatree->SetBranchAddress("tracklength",&event.tracklength);

    datatree->SetBranchAddress("rmsW",&event.rmsW);
    datatree->SetBranchAddress("rmsH",&event.rmsH);
  
    datatree->SetBranchAddress("charge",&event.charge);

    datatree->SetBranchAddress("isBeam",&event.isBeam);

    datatree->SetBranchAddress("helixnclusters",&event.helixnclusters);
    datatree->SetBranchAddress("helixxcluster",&event.helixxcluster);
    datatree->SetBranchAddress("helixycluster",&event.helixycluster);
    datatree->SetBranchAddress("helixzcluster",&event.helixzcluster);
   
    datatree->SetBranchAddress("gfnclusters",&event.gfnclusters);
    datatree->SetBranchAddress("gfxcluster",&event.gfxcluster);
    datatree->SetBranchAddress("gfycluster",&event.gfycluster);
    datatree->SetBranchAddress("gfzcluster",&event.gfzcluster);
    
    datatree->SetBranchAddress("nKstar",&event.nKstar);

    datatree->SetBranchAddress("chisqrKstar",&event.chisqrKstar);
    datatree->SetBranchAddress("fitprobKstar",&event.fitprobKstar);
    datatree->SetBranchAddress("mKstar",&event.mKstar);
    datatree->SetBranchAddress("mKstar_init",&event.mKstar_init);
  
    datatree->SetBranchAddress("xKstar",&event.xKstar);
    datatree->SetBranchAddress("yKstar",&event.yKstar);
    datatree->SetBranchAddress("zKstar",&event.zKstar);

    datatree->SetBranchAddress("pxKstar",&event.pxKstar);
    datatree->SetBranchAddress("pyKstar",&event.pyKstar);
    datatree->SetBranchAddress("pzKstar",&event.pzKstar);

    datatree->SetBranchAddress("nKs0",&event.nKs0);

    datatree->SetBranchAddress("chisqrKs0",&event.chisqrKs0);
    datatree->SetBranchAddress("fitprobKs0",&event.fitprobKs0);
    datatree->SetBranchAddress("mKs0",&event.mKs0);
    datatree->SetBranchAddress("mKs0_init",&event.mKs0_init);
  
    datatree->SetBranchAddress("xKs0",&event.xKs0);
    datatree->SetBranchAddress("yKs0",&event.yKs0);
    datatree->SetBranchAddress("zKs0",&event.zKs0);

    datatree->SetBranchAddress("pxKs0",&event.pxKs0);
    datatree->SetBranchAddress("pyKs0",&event.pyKs0);
    datatree->SetBranchAddress("pzKs0",&event.pzKs0);



  }

  InitializeTPC();  
  datatree->GetEntry(eventID);
  fStatusBar->SetText(Form("Data file %s Loaded sucessfully" , filename));

  return true;
}


Bool_t MyMainFrame::loadEvent() {

  printf("Loading Event %d.\n", eventID);

  fStatusBar->SetTextColor(0xff0000);
  fStatusBar->SetText(Form("Loading Event %d ......", eventID));
  //  gSystem->ProcessEvents();

  fTPC2dPoly->Reset("");
  fHTOF2dPoly->Reset("");
  fTPCxz->Reset("");
  fTPCxz_raw->Reset("");
  fTPCxz_helix->Reset("");
  fTPCxz_genfit->Reset("");

  fTPCKstar->Reset("");
  fTPCKs0->Reset("");
  
  datatree->GetEntry(eventID);
  
  Int_t nhit = event.nhittpc;
  //  Int_t ncluster = event.nhTpc;
  TPolyMarker3D *tpcRawHit3d = new TPolyMarker3D(nhit,2);

  for (int hit = 0 ; hit < nhit ; hit++) {
    Double_t x = event.xtpc[hit];
    Double_t y = event.ytpc[hit];
    Double_t z = event.ztpc[hit];
    Double_t charge = event.edeptpc[hit];
    Int_t pad = tpc::findPadID(z,x);
    Bool_t dead = tpc::Dead(pad);

    if (!dead)
      fTPC2dPoly->SetBinContent(pad+1,charge);
    else 
      fTPC2dPoly->SetBinContent(pad+1,0);

    tpcRawHit3d->SetPoint(hit,z,x,y);
    fTPCxz_raw->Fill(z,x);
    //    fTPCyz_raw->Fill(z,y);
  }
  
  Int_t numtrack = event.nGFTracks;
  
  if (numtrack == 0 ) return true;

  TArc *track_xz[numtrack];
  TLine *track_yz[numtrack];

  TPolyMarker3D *tpcCluster3d[numtrack];
  TPolyLine3D *tpcTrack3d[numtrack];
  
  for (int ntrack = 0 ; ntrack < numtrack ; ntrack++) {

    auto status = event.GFStatus->at(ntrack);
    //    if (status != 1) continue;
    
    Double_t radius = event.radius->at(ntrack);
    Double_t xcenter = event.xcenter->at(ntrack);
    Double_t zcenter = event.zcenter->at(ntrack);

    Double_t yoffset = event.yoffset->at(ntrack);
    Double_t slope = event.slope->at(ntrack);
    Double_t dipangle = event.dipangle->at(ntrack);
    
    Double_t alphaHead = event.alphaHead->at(ntrack);
    Double_t alphaTail = event.alphaTail->at(ntrack);
      
    Double_t HeadAngle = alphaHead*180./TMath::Pi();
    Double_t TailAngle = alphaTail*180./TMath::Pi();

    Double_t xhead = radius*TMath::Sin(alphaHead) + xcenter;
    Double_t zhead = radius*TMath::Cos(alphaHead) + zcenter;

    Double_t xtail = radius*TMath::Sin(alphaTail) + xcenter;
    Double_t ztail = radius*TMath::Cos(alphaTail) + zcenter;

    Double_t helixmom = event.helixmom->at(ntrack);

    Double_t rmsW = event.rmsW->at(ntrack);
    Double_t rmsH = event.rmsH->at(ntrack);

    Double_t charge = event.charge->at(ntrack);

    Double_t dip = event.dipangle->at(ntrack);
    
    //    track_xz->at(ntrack)= new TEllipse(zcenter,xcenter,radius);
    track_xz[ntrack]= new TArc(zcenter,xcenter,radius, HeadAngle, TailAngle);
    track_yz[ntrack] = new TLine(ztail,yoffset + slope*alphaTail,zhead,yoffset+slope*alphaHead);
    
    //    Int_t nPoints = 1./ratio;
    tpcTrack3d[ntrack]= new TPolyLine3D();

    TVector3 center(xcenter, yoffset, zcenter);
    for (auto i = 0. ; i <= 1. ; i += 0.1) {
      TVector3 ptOnHelix;
      Double_t alpha = i*alphaHead + (1-i)*alphaTail;
      Extrapolate(radius, center, slope, dip, alpha, ptOnHelix);
      //      std::cout << i << "\t" << ptOnHelix.X() << "\t" << ptOnHelix.Y() << "\t" << ptOnHelix.Z() << std::endl;
      tpcTrack3d[ntrack]->SetNextPoint(ptOnHelix.Z(), ptOnHelix.X(), ptOnHelix.Y());
      
    }
    
    Int_t nhelixcluster = event.helixnclusters->at(ntrack);
    if (nhelixcluster != 0) 
      tpcCluster3d[ntrack] = new TPolyMarker3D(nhelixcluster,8);

    for (int icl = 0 ; icl < nhelixcluster ; icl++) {

      Double_t xpos = event.helixxcluster->at(ntrack).at(icl);
      Double_t ypos = event.helixycluster->at(ntrack).at(icl);
      Double_t zpos = event.helixzcluster->at(ntrack).at(icl);
      
      //	  tpcTrack->at(ntrack)->SetPoint(numStableCluster,ztrcluster,xtrcluster,ytrcluster);
      fTPCxz_helix->Fill(zpos,xpos);
      fTPCxz_helix->SetMarkerColor(kSpring);
      tpcCluster3d[ntrack]->SetPoint(icl,zpos, xpos, ypos);
      tpcCluster3d[ntrack]->SetMarkerColor(kSpring);
    }

    Int_t nCl = event.gfnclusters->at(ntrack);
    for (int iCl = 0 ; iCl < nCl ; iCl++) {

      auto xpos = event.gfxcluster->at(ntrack).at(iCl);
      auto ypos = event.gfycluster->at(ntrack).at(iCl);
      auto zpos = event.gfzcluster->at(ntrack).at(iCl);

      fTPCxz_genfit->Fill(zpos,xpos);
      fTPCxz_genfit->SetMarkerColor(kOrange);
    }
  }
  
  Int_t nKstar = event.nKstar;
  Int_t nKs0 = event.nKs0;

  TPolyMarker3D *tpcKstarvtx = new TPolyMarker3D(1,2);
  TPolyMarker3D *tpcKs0vtx = new TPolyMarker3D(1,2);

  TVector3 KstarVtx(-999,-999,-999);
  TVector3 Ks0Vtx(-999,-999,-999);
  
  if (nKstar == 1) {
    for (int n = 0 ; n < nKstar ; n++) {
      Double_t fitprob = event.fitprobKstar->at(n);
      //      hkinprobKstar->Fill(fitprob);
      if (fitprob > 0.01) {
	
	Double_t xKstar = event.xKstar->at(n);
	Double_t yKstar = event.yKstar->at(n);
	Double_t zKstar = event.zKstar->at(n);
	KstarVtx.SetXYZ(xKstar,yKstar,zKstar);
	tpcKstarvtx->SetPoint(n,zKstar,xKstar,yKstar);
	fTPCKstar->Fill(zKstar,xKstar);
      }
    }
  }
  if (nKs0 > 0) {
    for (int n = 0 ; n < 1 ; n++) {
      Double_t fitprob = event.fitprobKs0->at(n);
      //      hkinprobKs0_all->Fill(fitprob);
      if (fitprob < 0.01) continue;
      Double_t xKs0 = event.xKs0->at(n);
      Double_t yKs0 = event.yKs0->at(n);
      Double_t zKs0 = event.zKs0->at(n);
      Ks0Vtx.SetXYZ(xKs0,yKs0,zKs0);
      if (n == 0) tpcKs0vtx->SetPoint(n,zKs0,xKs0,yKs0);	
      fTPCKs0->Fill(zKs0,xKs0);    
    }
  }
  //  tree->Print();
  
  fCanvas->cd(1);
  
  fTPC2dPoly->Draw("col");
  fHTOF2dPoly->Draw("colsame");

  fTPCxz_raw->Draw("SAME");
  fTPCxz_raw->SetMarkerSize(1);
  fTPCxz_raw->SetMarkerStyle(2);
  fTPCxz_raw->SetMarkerColor(kOrange-2);

  fTPCxz->Draw("SAME");
  fTPCxz->SetMarkerStyle(8);
  fTPCxz->SetMarkerColor(kRed);

  fTPCxz_helix->Draw("SAME");
  fTPCxz_helix->SetMarkerStyle(8);
  fTPCxz_genfit->Draw("SAME");
  fTPCxz_genfit->SetMarkerStyle(8);
  //  fTPCxz_cluster->SetMarkerColor(kSpring);

  fTPCKstar->Draw("SAME");
  fTPCKstar->SetMarkerSize(1);
  fTPCKstar->SetMarkerStyle(30);
  fTPCKstar->SetMarkerColor(kRed);

  fTPCKs0->Draw("SAME");
  fTPCKs0->SetMarkerSize(1);
  fTPCKs0->SetMarkerStyle(33);
  fTPCKs0->SetMarkerColor(kMagenta);

  gPad->SetLogz(1);
  
  fCanvas->Modified();
  fCanvas->Update();

  for (int nt = 0 ; nt < numtrack ; nt++) {
    //    if (event.isBeam->at(nt) == -1) continue;
    track_xz[nt]->Draw("SAMEONLY");
    track_xz[nt]->SetFillStyle(0);
    track_xz[nt]->SetFillColorAlpha(10, 0);
    if (event.charge->at(nt) > 0)
      track_xz[nt]->SetLineColor(kAzure);
    else 
      track_xz[nt]->SetLineColor(kRed);
    track_xz[nt]->SetLineWidth(2);
  }

  fCanvas->Update();
  fCanvas->Modified();

  fCanvas->cd(2);
  TView3D *view = (TView3D*) TView::CreateView(1);
  fTPC->Draw();
  view->ZoomView(gPad,1.5);
  
  for (int nt = 0 ; nt < numtrack ; nt++){
    //    tpcRawHit3d->Draw("SAME");
    //    tpcRawHit3d->SetMarkerSize(0.7);
    //    tpcRawHit3d->SetMarkerColor(kOrange-2);

    //    if (event.isBeam->at(nt) == -1) continue;
    //    tpcHit3d[nt]->Draw("SAME");
    //    tpcHit3d[nt]->SetMarkerSize(0.7);
    //    tpcHit3d[nt]->SetMarkerColor(kRed);
    //    tpcHitTrack3d[nt]->Draw("SAME");
    tpcCluster3d[nt]->Draw("SAME");
    tpcCluster3d[nt]->SetMarkerSize(0.7);

    //    tpcTrack[nt]->Draw("SAME");
    //    tpcTrack[nt]->SetLineWidth(2);
    //    tpcTrack[nt]->SetLineColor(nt + 2);

    tpcTrack3d[nt]->Draw("SAME");
    tpcTrack3d[nt]->SetLineWidth(2);
    tpcTrack3d[nt]->SetLineColor(nt + 2);
//    tpcClusterTrack3d[nt]->Draw("SAME");
    //    if (event.charge[nt] > 0)
    //      tpcClusterTrack3d[nt]->SetLineColor(kAzure);
    //    else
    //      tpcClusterTrack3d[nt]->SetLineColor(kRed);
   
  }
  tpcKstarvtx->Draw("SAMEONLY");
  tpcKstarvtx->SetMarkerSize(1.2);
  tpcKstarvtx->SetMarkerStyle(30);
  tpcKstarvtx->SetMarkerColor(kRed);

  tpcKs0vtx->Draw("SAMEONLY");
  tpcKs0vtx->SetMarkerSize(1.2);
  tpcKs0vtx->SetMarkerStyle(33);
  tpcKs0vtx->SetMarkerColor(kMagenta);
  fCanvas->Update();
  fCanvas->Modified();
    
  fStatusBar->SetText(Form("Event %d Loaded sucessfully" , eventID));

  return true;
  
}

void MyMainFrame::GoForward() {

  fCanvas->Clear("D");
  
  if (eventID < datatree -> GetEntries() -1) {
    eventID++;
    loadEvent();
  }
  else {
    fStatusBar->SetTextColor(0xff0000);
    fStatusBar->SetText("Already at the last event");
    printf("Already at the last event.\n");
  }
}

void MyMainFrame::GoBackward() {

  fCanvas->Clear("D");

  //  fTPC->CleanAll();
  
  if (eventID > 0 ) {
    eventID--;
    loadEvent();
  }
  else {
    fStatusBar->SetTextColor(0xff0000);
    fStatusBar->SetText("Already at the first event");
    printf("Already at the first event.\n");
  }
}

void MyMainFrame::GoTo() {

  fCanvas->Clear("D");

  //  fTPC->CleanAll();
  
  TString evnum = fEvtHandler->GetMarkedText();
  eventID = evnum.Atoi();
  if (eventID >=0 && eventID < datatree->GetEntries()) {
    loadEvent();
  }
  else {
    fStatusBar->SetTextColor(0xff0000);
    fStatusBar->SetText("Out of range");
    printf("Out of range\n");
  }
}

void MyMainFrame::CloseWindow()
{
  gApplication->Terminate();
}

void MyMainFrame::HandleMenu(Int_t id)
{
   // Handle menu items.

   switch (id) {

      case M_FILE_OPEN:
         {
            static TString dir(".");
            TGFileInfo fi;
            fi.fFileTypes = filetypes;
            fi.SetIniDir(dir);
            printf("fIniDir = %s\n", fi.fIniDir);
            new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
            printf("Open file: %s (dir: %s)\n", fi.fFilename, fi.fIniDir);
            dir = fi.fIniDir;
         }
         break;

      case M_FILE_SAVE:
         printf("M_FILE_SAVE\n");
         break;

      case M_FILE_EXIT:
         CloseWindow();   // terminate theApp no need to use SendCloseMessage()
         break;

   case M_ANA_RECO:
         new Analysis(gClient->GetRoot(), fMain, 400, 200);
         break;

	 /*
      case M_ANA_FIT
         new TestShutter(gClient->GetRoot(), fMain, 400, 200);
         break;
	 */
      case M_VIEW_ENBL_DOCK:
         fMenuDock->EnableUndock(!fMenuDock->EnableUndock());
         if (fMenuDock->EnableUndock()) {
            fMenuView->CheckEntry(M_VIEW_ENBL_DOCK);
            fMenuView->EnableEntry(M_VIEW_UNDOCK);
         } else {
            fMenuView->UnCheckEntry(M_VIEW_ENBL_DOCK);
            fMenuView->DisableEntry(M_VIEW_UNDOCK);
         }
         break;

      case M_VIEW_ENBL_HIDE:
         fMenuDock->EnableHide(!fMenuDock->EnableHide());
         if (fMenuDock->EnableHide()) {
            fMenuView->CheckEntry(M_VIEW_ENBL_HIDE);
         } else {
            fMenuView->UnCheckEntry(M_VIEW_ENBL_HIDE);
         }
         break;

       case M_VIEW_DOCK:
         fMenuDock->DockContainer();
         fMenuView->EnableEntry(M_VIEW_UNDOCK);
         fMenuView->DisableEntry(M_VIEW_DOCK);
         break;

       case M_VIEW_UNDOCK:
         fMenuDock->UndockContainer();
         fMenuView->EnableEntry(M_VIEW_DOCK);
         fMenuView->DisableEntry(M_VIEW_UNDOCK);
         break;

      default:
         printf("Menu item %d selected\n", id);
         break;
   }
}
TGeoVolume* MyMainFrame::TPCGeometry() {

    //TPC Drawing
  Double_t flength = 586;
  Double_t fheight = 550;
  Double_t edge = flength/(1+sqrt(2));
  Double_t tan = TMath::Tan(22.5*TMath::Pi()/180.);
  
  TGeoManager *geom = new TGeoManager("K1.8HS","K1.8HS");
  TGeoMaterial *mat_world = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium *med_world = new TGeoMedium("med_world",1,mat_world);
  TGeoVolume *world = geom->MakeBox("world",med_world,flength*2,flength*2,flength*2);
  geom->SetTopVolume(world);
  world->SetVisibility(kFALSE);
  
  TGeoMaterial *mat_tpc = new TGeoMaterial("Al",26.98,13,2.7);
  TGeoMedium *med_tpc = new TGeoMedium("med_tpc",1,mat_tpc);
  TGeoXtru *xtru_tpc_inside = new TGeoXtru(2);
  Double_t zin[] = { -tan*flength/2, -flength/2, -flength/2, -tan*flength/2,
		     tan*flength/2, flength/2, flength/2, tan*flength/2};
  Double_t xin[] = { -flength/2, -edge/2, edge/2, flength/2,
		   flength/2, edge/2, -edge/2, -flength/2};  
  Double_t yin[] = { -fheight/2, fheight/2};
  Double_t scale[] = {1., 1.};
  Double_t x0[] = {0., 0.};
  Double_t z0[] = {0., 0.};

  Int_t nxz = sizeof(xin)/sizeof(Double_t);
  xtru_tpc_inside->DefinePolygon(nxz,xin,zin);
  
  Int_t n;
  Int_t ny = sizeof(yin)/sizeof(Double_t);
  for (n = 0 ; n < ny ; n++)
    xtru_tpc_inside->DefineSection(n,yin[n],x0[n],z0[n],scale[n]);  

  TGeoVolume *tpc_inside = new TGeoVolume("tpc_inside",xtru_tpc_inside,med_tpc);

  ////////////////////////////////////////////////////////////////////////////
  flength = 586 - 10;
  edge = flength/(1+sqrt(2));

  TGeoXtru *xtru_tpc_outside = new TGeoXtru(2);
  Double_t z[] = { -tan*flength/2, -flength/2, -flength/2, -tan*flength/2,
		     tan*flength/2, flength/2, flength/2, tan*flength/2};
  Double_t x[] = { -flength/2, -edge/2, edge/2, flength/2,
		   flength/2, edge/2, -edge/2, -flength/2};  
  Double_t y[] = { -fheight/2, fheight/2};
  
  xtru_tpc_outside->DefinePolygon(nxz,x,z);  
  for (n = 0 ; n < ny ; n++)
    xtru_tpc_outside->DefineSection(n,y[n],x0[n],z0[n],scale[n]);  
  
  TGeoVolume *tpc_outside = new TGeoVolume("tpc_outside",xtru_tpc_outside,med_tpc);

  //  new TGeoManager("Target","Tgt inside TPC");
  TGeoMaterial *mat_tgt = new TGeoMaterial("C",24,12,1.3);
  TGeoMedium *med_tgt = new TGeoMedium("med_tgt",1,mat_tgt);

  TGeoBBox *box_tgt = new TGeoBBox("box_tgt",30/2,20/2,20/2);
  TGeoVolume *tgt = new TGeoVolume("tgt",box_tgt,med_tgt);
  //  TGeoVolume *tgt = geom->MakeBox("tgt",med_tgt,20/2,30/2,20/2);

  //  tpc->SetLineColor(kBlack);
  //  tpc->SetLineWidth(1);
  TGeoVolumeAssembly *tpc = new TGeoVolumeAssembly("tpc");  
  tpc->AddNode(tgt,1,new TGeoTranslation(-143.,0.,0.));
  tpc->AddNode(tpc_outside,1);
  tpc->AddNode(tpc_inside,1);
  
  world->AddNode(tpc,1);

  geom->CloseGeometry();
  geom->SetVisLevel(4);

  return world;
}


Analysis::Analysis(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h) {

  fMain = new TGTransientFrame(p, main, w, h);
  fMain->Connect("CloseWindow()","Analysis",this,"DoClose()");
  fMain->DontCallClose();
  fMain->SetCleanup(kDeepCleanup);

  fFrame1 = new TGHorizontalFrame(fMain, 60, 20, kFixedWidth);
  fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);
  fOkButton = new TGTextButton(fFrame1, "&OK", 1);
  fOkButton->Connect("Clicked()", "Analysis", this, "DoOK()");
  fCancelButton = new TGTextButton(fFrame1, "&Cancel", 2);
  fCancelButton->Connect("Clicked()", "Analysis", this, "DoCancel()");
  fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 2, 2, 2, 2);
  fL2 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1);
  fFrame1->AddFrame(fOkButton, fL1);
  fFrame1->AddFrame(fCancelButton, fL1);
  fFrame1->Resize(150, fOkButton->GetDefaultHeight());
  fMain->AddFrame(fFrame1, fL2);
  
  fFillHistos = kFALSE;
  fHpx1 = 0;
  fHpx2 = 0;

  fFrame2 = new TGCompositeFrame(fMain, 60, 20, kHorizontalFrame);
  fStartButton = new TGTextButton(fFrame2, "Start Filling Histogram");
  fStopButton = new TGTextButton(fFrame2, "Stop Filling Histogram");
  fStartButton->Connect("Clicked()","Analysis",this,"HandleButtons()");
  fStopButton->Connect("Clicked()","Analysis",this,"HandleButtons()");
  fFrame2->AddFrame(fStartButton,fL3);
  fFrame2->AddFrame(fStopButton,fL3);

  fFrame3 = new TGCompositeFrame(fMain, 60,60, kHorizontalFrame);
  fL4 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5);

  fECanvas1 = new TRootEmbeddedCanvas("eCanvas1", fFrame3, 100, 100);
  fFrame3->AddFrame(fECanvas1, fL4);
  fECanvas2 = new TRootEmbeddedCanvas("eCanvas2", fFrame3, 100, 100);
  fFrame3->AddFrame(fECanvas2, fL4);

  fECanvas1->GetCanvas()->SetBorderMode(0);
  fECanvas2->GetCanvas()->SetBorderMode(0);
  fECanvas1->SetBit(kNoContextMenu);

  TGLayoutHints *fL5 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX | kLHintsExpandY, 2, 2, 5, 1);

  fMain->AddFrame(fFrame2, fL3);
  fMain->AddFrame(fFrame3, fL4);
  //  fMain->AddFrame(fFrame1, fL5);

  fMain->MapSubwindows();
  fMain->Resize();

  fMain->CenterOnParent();
  fMain->SetWindowName("Analysis");
  fMain->MapWindow();
}

Analysis::~Analysis(){

  fMain->DeleteWindow();
}
void Analysis::DoOK() {
  
  fFillHistos = kFALSE;
  printf("\nTerminating dialog: OK pressed\n");
  // Add protection against double-clicks
  fOkButton->SetState(kButtonDisabled);
  fCancelButton->SetState(kButtonDisabled);
  
  // Send a close message to the main frame. This will trigger the
  // emission of a CloseWindow() signal, which will then call
  // TestDialog::CloseWindow(). Calling directly CloseWindow() will cause
  // a segv since the OK button is still accessed after the DoOK() method.
  // This works since the close message is handled synchronous (via
  // message going to/from X server).
  fMain->SendCloseMessage();
  
  // The same effect can be obtained by using a singleshot timer:
//  TTimer::SingleShot(150, "Analysis", this, "CloseWindow()");
  
  // Close the Ged editor if it was activated.
  if (TVirtualPadEditor::GetPadEditor(kFALSE) != 0)
    TVirtualPadEditor::Terminate();
}

void Analysis::DoCancel()
{
   fFillHistos = kFALSE;
   printf("\nTerminating dialog: Cancel pressed\n");
   // Add protection against double-clicks
   fOkButton->SetState(kButtonDisabled);
   fCancelButton->SetState(kButtonDisabled);
   TTimer::SingleShot(150, "Analysis", this, "CloseWindow()");
   // Close the Ged editor if it was activated.
   if (TVirtualPadEditor::GetPadEditor(kFALSE) != 0)
      TVirtualPadEditor::Terminate();
}

void Analysis::DoClose()
{
   printf("\nTerminating dialog: via window manager\n");
   if (fFillHistos) {
      fFillHistos = kFALSE;
      TTimer::SingleShot(150, "TestDialog", this, "CloseWindow()");
   } else
      CloseWindow();

   // Close the Ged editor if it was activated.
   if (TVirtualPadEditor::GetPadEditor(kFALSE) != 0)
      TVirtualPadEditor::Terminate();
}

void Analysis::CloseWindow() {

  delete this;
}


void GFTPCEventDisplay() {

  gSystem->Load("../../macro/AutoDict_std__vector_std__vector_std__vector_Double_t____cxx.so");
  
  //  gStyle->SetPalette(kViridis);
  gStyle->SetPalette(kCividis);
  gStyle->SetNumberContours(255);
   // Popup the GUI...
   new MyMainFrame(gClient->GetRoot(),1600,800);

}


Double_t Extrapolate(Double_t radius, TVector3 center, Double_t slope, Double_t dip, Double_t alpha, TVector3 &ptOnHelix) {

  ptOnHelix.SetXYZ(radius * TMath::Sin(alpha) + center.X(),
		   alpha * slope + center.Y(),
		   radius * TMath::Cos(alpha) + center.Z());
  Double_t length = alpha * radius / TMath::Cos(dip);

  return length;
}
