#include <TTree.h>
#include <TFile.h>
#include <iostream>

#include "BcOutEventLoader.hh"
#include "TPCPadHelper.hh"
#include "EventData.hh"

using namespace std;

void LoadBcOutEventData(TTree *tree, int evnum, BcOutEventData &data) {
  
  tree->SetBranchAddress("ntrack",&data.ntrack);
  tree->SetBranchAddress("x0",&data.x0);
  tree->SetBranchAddress("y0",&data.y0);
  tree->SetBranchAddress("u0",&data.u0);
  tree->SetBranchAddress("v0",&data.v0);

  tree->GetEvent(evnum);
}
