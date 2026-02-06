#ifndef EVENTLOADER_HH
#define EVENTLOADER_HH

#include "EventData.hh"
#include <TTree.h>

void LoadEventData(TTree *tree, int evnum, EventData &data);

#endif
