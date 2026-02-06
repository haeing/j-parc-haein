#ifndef G4EVENTLOADER_HH
#define G4EVENTLOADER_HH

#include "G4EventData.hh"
#include <TTree.h>

void LoadEventDataG4(TTree *tree, int evnum, G4EventData &data);

#endif
