#include <TRootEmbeddedCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TQObject.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TPaveText.h>
#include <TEllipse.h>
#include <TPolyMarker3D.h>
#include <TArc.h>
#include <THelix.h>
#include <TPolyLine3D.h>
#include <RQ_OBJECT.h>
#include <TVector3.h>


#include "TPC3DFrameBuilder.hh"
#include "TPC2DFrameBuilder.hh"
#include "TPCPadHelper.hh"
#include "EventLoader.hh"
#include "EventData.hh"
#include "BcOutEventLoader.hh"
#include "BcOutEventData.hh"
//#include "G4EventLoader.hh"
//#include "G4EventData.hh"

const string pdg[10] = {"e^{-}","e^{+}","#pi^{+}","#pi^{-}","K^{+}","K^{-}","p","#mu^{+}","#mu^{-}","other"};
const int pdg_code[10] = {11,-11,211,-211,321,-321,2212,-13,13,-9999};
const int pdg_color[10] = {1,44,2,3,29,6,7,8,9,14};


class TGWindow;
class TGMainFrame;

class DstEventDisplay_E42 {
  RQ_OBJECT("DstEventDisplay_E42")

public:
  DstEventDisplay_E42(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~DstEventDisplay_E42();
  void GoForward();
  void GoBackward();
  void GoTo();
  void DoCanvasDraw();
  void DrawTrackG4();
  void DrawTrackHit();
  void DrawVertex();
  
private:
  TGMainFrame *fMain;
  TRootEmbeddedCanvas *fEcanvas;

  TGButton *pre_Button, *next_Button, *go_Button;
  TGTextEntry      *fEvtHandler;
  
  TFile *fFile;
  TTree *fTree;
  TFile *fBcOutFile;
  TTree *fBcOutTree;

  //TFile *fFile_g4;
  //TTree *fTree_g4;

  TCanvas *fCanvas;
  TGraph *Track_g4[MaxTrack + 1];
  
  TGraph *Vertex_PRM_g4;
  TGraph *Vertex_SEC_g4;
  
  TGraph *Track_hit2D[MaxTrack];
  TPolyMarker3D *Track_hit3D[MaxTrack];
  TArc *Track_xz[MaxTrack];
  //THelix *Track_3D[MaxTrack];
  TPolyLine3D *Track_3D[MaxTrack];
  TPolyLine3D *BcOutTrack_3D[MaxTrack];
  EventData fEvent;
  BcOutEventData fBcOutEvent;
  
  TH2Poly *TPC_2d_xz;
  TGeoVolume *fTPC;
  TEllipse *target_zx;
  
  TPC2DFrameBuilder builder_2d;
  TPC3DFrameBuilder builder_3d;
  int fNtrack_g4;

  TPaveText *fDef;

  int fCurrentEvent;
  
};
