#ifdef __CINT__
#pragma link C++ enum EMessageTypes;
#pragma link C++ enum EMainCommandIdentifiers;
#endif


#include <vector>  
#include <iostream>

#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TPad.h>
#include <TSystem.h>
#include <TGeoManager.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TButton.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TGClient.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <TF1.h>
#include <TRandom.h>
#include <TPolyLine3D.h>

#include "TPC3DFrameBuilder.hh"
#include "TPC2DFrameBuilder.hh"
#include "TPCPadHelper.hh"
#include "EventLoader.hh"
#include "EventData.hh"
#include "BcOutEventLoader.hh"
#include "BcOutEventData.hh"
#include "DstEventDisplay_E42.h"

namespace{
using namespace std;

}

void DstEventDisplay_E42::DrawTrackHit(){

  double x_bcout = 0;
  double y_bcout = 0;
  //BcOut
  for(int i=0;i<fBcOutEvent.ntrack;i++){
    double x_ff = (*fBcOutEvent.x0)[i];
    double y_ff = (*fBcOutEvent.y0)[i];

    double ubcout = (*fBcOutEvent.u0)[i];
    double vbcout = (*fBcOutEvent.v0)[i];

    double x_start = x_ff + ubcout*(-400);
    double y_start = y_ff + vbcout*(-400);

    double x_end = x_ff;
    double y_end = y_ff;

    BcOutTrack_3D[i] = new TPolyLine3D(2);
    BcOutTrack_3D[i]->SetPoint(0,-x_start,-400,y_start);
    BcOutTrack_3D[i]->SetPoint(1,-x_end,0,y_end);
    x_bcout = x_ff;
    y_bcout = y_ff;
  }
  
  
  int nhitpoint = 0;
  
  for(int i=0;i<fEvent.ntTpc;i++){
    Track_hit2D[i] = new TGraph();
    Track_hit3D[i] = new TPolyMarker3D(MaxPad,8);

    Track_hit2D[i]->SetMarkerStyle(20);
    Track_hit2D[i]->SetMarkerColor(pdg_color[i]);
    Track_hit2D[i]->SetMarkerSize(0.7);

    Track_hit3D[i]->SetMarkerStyle(20);
    Track_hit3D[i]->SetMarkerSize(0.5);
    Track_hit3D[i]->SetMarkerColor(pdg_color[i]);
    
    nhitpoint = 0;
    int nhittrack = (*fEvent.nhtrack)[i];

    /*
    double mint = (*fEvent.helix_t)[i][0];
    double maxt = (*fEvent.helix_t)[i][nhittrack-1];
    double min_t = mint*TMath::RadToDeg() + 270.;
    double max_t = maxt*TMath::RadToDeg() + 270.;

    Track_xz[i] = new TArc((*fEvent.helix_cy)[i]-143.,-1*(*fEvent.helix_cx)[i],(*fEvent.helix_r)[i],min_t,max_t);

    double vx0 = -1*(*fEvent.helix_r)[i]*TMath::Sin(mint);
    double vy0 = (*fEvent.helix_r)[i]*TMath::Cos(mint);
    double vz0 = (*fEvent.helix_r)[i]*(*fEvent.helix_dz)[i];

    double x0 = (*fEvent.helix_cx)[i] + (*fEvent.helix_r)[i]*TMath::Cos(mint);
    double y0 = (*fEvent.helix_cy)[i] + (*fEvent.helix_r)[i]*TMath::Sin(mint) - 143.;
    double z0 = (*fEvent.helix_z0)[i] + vz0*mint;


    double initvel[3] = {vx0,vy0,vz0};
    double initpos[3] = {x0,y0,z0};
    
    double range[2] = {z0,z0 + vz0*(maxt-mint)};
    

    Track_3D[i] = new THelix(initpos,initvel,-1.,range);
    */
    double x_target = (*fEvent.x0Tpc)[i];
    double y_target = (*fEvent.y0Tpc)[i];

    double utpc = (*fEvent.u0Tpc)[i];
    double vtpc = (*fEvent.v0Tpc)[i];

    double x_start = x_target + utpc*(-200);
    double y_start = y_target + vtpc*(-200);

    double x_end = x_target + utpc*(400);
    double y_end = y_target + vtpc*(400);
    
    double xtpc = (x_target) + (143+150)*utpc;
    double ytpc = (y_target) + (143+150)*vtpc;

    double x_res = x_bcout - xtpc;
    double y_res = y_bcout - ytpc;

    
    Track_3D[i] = new TPolyLine3D(2);
    Track_3D[i]->SetPoint(0,-x_start,-200-143,y_start+64.45);
    Track_3D[i]->SetPoint(1,-x_end,400-143,y_end+64.45);
    if(x_res > -10 && x_res < 10 && y_res > 55 && y_res < 80){
      Track_3D[i]->SetLineColor(kBlue);
    }
    else{Track_3D[i]->SetLineColor(kGray);}

    
    for(int j=0;j<nhittrack;j++){
      if((*fEvent.hitpos_x)[i][j]==-9999 || (*fEvent.hitpos_z)[i][j]==-9999)continue;

      
      Track_hit2D[i]->SetPoint(Track_hit2D[i]->GetN(),(*fEvent.hitpos_z)[i][j],(*fEvent.hitpos_x)[i][j]);
      Track_hit3D[i]->SetPoint(nhitpoint,-(*fEvent.hitpos_x)[i][j],(*fEvent.hitpos_z)[i][j],(*fEvent.hitpos_y)[i][j]+64.45);

      nhitpoint++;
    }
    
  }

  
}

void DstEventDisplay_E42::DoCanvasDraw(){

  fCanvas->cd(2);
  fDef->Draw();
  
  fCanvas->cd(3);
  
  TPC_2d_xz->Draw();
  for(int i=0;i<fEvent.ntTpc;i++){

    Track_hit2D[i]->Draw("same P");
    /*
    Track_xz[i]->SetLineColor(pdg_color[i]);
    Track_xz[i]->SetLineWidth(2);
    Track_xz[i]->SetFillStyle(0);
    Track_xz[i]->SetFillColorAlpha(10,0);
    Track_xz[i]->Draw("sameonly");
    */
  }
  
  target_zx->Draw("same");
  
  fCanvas->cd(4);
  
  fTPC->Draw();
  
  for(int i=0;i<fEvent.ntTpc;i++){
    Track_hit3D[i]->Draw("SAME");
    //Track_3D[i]->SetLineColor(pdg_color[i]);
    Track_3D[i]->SetLineWidth(2);
    Track_3D[i]->Draw("same P");
  }
  for(int i=0;i<fBcOutEvent.ntrack;i++){
    BcOutTrack_3D[i]->SetLineColor(kRed);
    BcOutTrack_3D[i]->SetLineWidth(2);
    BcOutTrack_3D[i]->Draw("same P");
  }

  fCanvas->Update();
  
}

DstEventDisplay_E42::DstEventDisplay_E42(const TGWindow* p, UInt_t w, UInt_t h) : fCurrentEvent(0) {

  fMain = new TGMainFrame(p,w,h);
  
  fFile = new TFile("../../data/run02570_DstTPCTracking.root");
  
  fTree = (TTree*)fFile->Get("tpc");

  fBcOutFile = new TFile("../../data/run02570_BcOutTracking.root");
  
  fBcOutTree = (TTree*)fBcOutFile->Get("bcout");
  
  
  LoadEventData(fTree, fCurrentEvent, fEvent);
  LoadBcOutEventData(fBcOutTree,fCurrentEvent,fBcOutEvent);

  DrawTrackHit();
  
  TPC_2d_xz = builder_2d.TPC2DGeometry();
  fTPC = builder_3d.TPC3DGeometry();

  target_zx = new TEllipse(-143.,0.,40.);
  target_zx->SetLineColor(kRed); 
  target_zx->SetLineWidth(2);   
  target_zx->SetFillStyle(0);
  
  fDef = new TPaveText(0.4214918,0.0538522,0.5775999,0.9867662,"blNDC");
  fDef->SetName("title");
  fDef->SetBorderSize(0);
  fDef->SetFillColor(0);
  fDef->SetTextFont(42);
  TText *Def_pid[10];
  for(int i=0;i<10;i++){
    Def_pid[i] = fDef->AddText(Form("%s",pdg[i].c_str()));
    Def_pid[i]->SetTextColor(pdg_color[i]);
  }
  
  fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,1000,800);
 

  fCanvas = fEcanvas->GetCanvas();
  fCanvas->Divide(2,2);
  DoCanvasDraw();
  
  
  fMain->AddFrame(fEcanvas,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));
  TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,200,40);
  
  pre_Button = new TGTextButton(hframe,"&Previous",1);
  pre_Button->Connect("Clicked()","DstEventDisplay_E42",this,"GoBackward()");
  hframe->AddFrame(pre_Button, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

  next_Button = new TGTextButton(hframe, "&Next",1);
  next_Button->Connect("Clicked()","DstEventDisplay_E42",this,"GoForward()");
  hframe->AddFrame(next_Button, new TGLayoutHints(kLHintsCenterX ,5,5,3,4));


  fEvtHandler = new TGTextEntry(hframe);
  fEvtHandler ->SetEnabled(true);
  fEvtHandler->SetInsertMode();
  hframe->AddFrame(fEvtHandler, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));
  
  go_Button = new TGTextButton(hframe, "&Go");
  go_Button->Connect("Clicked()", "DstEventDisplay_E42", this, "GoTo()");
  hframe->AddFrame(go_Button, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));

  

  TGTextButton *exit = new TGTextButton(hframe,"&Exit","gApplication->Terminate(0)");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
  fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

   
  fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
  fMain->SetWindowName("HypTPC Event Display");
  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->MapWindow();
}

void DstEventDisplay_E42::GoBackward() {

  fCanvas->Clear("D");
  
  if (fCurrentEvent == 0) {
    fCurrentEvent = 0;
    std::cout << "First Event" << std::endl;
  } else {
    fCurrentEvent--;
  }
  
  LoadEventData(fTree, fCurrentEvent, fEvent);
  LoadBcOutEventData(fBcOutTree,fCurrentEvent,fBcOutEvent);
  DrawTrackHit();
  DoCanvasDraw();

  
}

void DstEventDisplay_E42::GoForward(){
  if(fCurrentEvent == fTree->GetEntries()-1){
    std::cout<< "Last Event" <<std::endl;
  }
  else{
    fCurrentEvent++;
  }

  LoadEventData(fTree, fCurrentEvent, fEvent);
  LoadBcOutEventData(fBcOutTree,fCurrentEvent,fBcOutEvent);
  DrawTrackHit();
  DoCanvasDraw();

}

void DstEventDisplay_E42::GoTo(){
  
  TString eventNumberStr = fEvtHandler ->GetMarkedText();
  Int_t evnum = eventNumberStr.Atoi();

  if(evnum < 0 || evnum >= fTree->GetEntries()){
    cout<<"No Event" <<endl;
  }

  else{
    fCurrentEvent = evnum;
    LoadEventData(fTree, fCurrentEvent, fEvent);
    LoadBcOutEventData(fBcOutTree,fCurrentEvent,fBcOutEvent);
    DrawTrackHit();
    DoCanvasDraw();
  }

}


DstEventDisplay_E42::~DstEventDisplay_E42(){
  
  fMain->Cleanup();
  delete fMain;
}
    

int main(int argc, char **argv){

  gStyle->SetOptStat(0);
  
  TApplication theApp("App", &argc, argv);
  new DstEventDisplay_E42(gClient->GetRoot(),1000,800);
  theApp.Run();  
  
}
