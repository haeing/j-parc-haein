#ifndef BCOUTEVENTLOADER_HH
#define BCOUTEVENTLOADER_HH

#include "BcOutEventData.hh"
#include <TTree.h>

void LoadBcOutEventData(TTree *tree, int evnum, BcOutEventData &data);

#endif
