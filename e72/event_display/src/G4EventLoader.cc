#include <TTree.h>
#include <TFile.h>
#include <iostream>

#include "G4EventLoader.hh"
#include "TPCPadHelper.hh"
#include "EventData.hh"

using namespace std;

void LoadEventDataG4(TTree *tree, int evnum, G4EventData &data) {

  tree->SetBranchAddress("BEAM",&data.BEAM);
  tree->SetBranchAddress("PRM",&data.PRM);
  tree->SetBranchAddress("SEC",&data.SEC);
  
  tree->GetEvent(evnum);
}
