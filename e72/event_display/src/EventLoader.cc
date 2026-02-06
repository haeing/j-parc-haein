#include <TTree.h>
#include <TFile.h>
#include <iostream>

#include "EventLoader.hh"
#include "TPCPadHelper.hh"
#include "EventData.hh"

using namespace std;

void LoadEventData(TTree *tree, int evnum, EventData &data) {
  
  tree->SetBranchAddress("nhittpc",&data.nhittpc);
  tree->SetBranchAddress("ntTpc",&data.ntTpc);

  tree->SetBranchAddress("nhtrack",&data.nhtrack);

  tree->SetBranchAddress("x0Tpc",&data.x0Tpc);
  tree->SetBranchAddress("y0Tpc",&data.y0Tpc);
  tree->SetBranchAddress("u0Tpc",&data.u0Tpc);
  tree->SetBranchAddress("v0Tpc",&data.v0Tpc);
  
  /*
  tree->SetBranchAddress("trackid",&data.trackid);

  tree->SetBranchAddress("helix_cx",&data.helix_cx);
  tree->SetBranchAddress("helix_cy",&data.helix_cy);
  tree->SetBranchAddress("helix_z0",&data.helix_z0);
  tree->SetBranchAddress("helix_r",&data.helix_r);
  tree->SetBranchAddress("helix_dz",&data.helix_dz);
  tree->SetBranchAddress("dz_factor",&data.dz_factor);

  tree->SetBranchAddress("helix_t",&data.helix_t);
  */
  tree->SetBranchAddress("hitpos_x",&data.hitpos_x);
  tree->SetBranchAddress("hitpos_y",&data.hitpos_y);
  tree->SetBranchAddress("hitpos_z",&data.hitpos_z);

  tree->GetEvent(evnum);
}
