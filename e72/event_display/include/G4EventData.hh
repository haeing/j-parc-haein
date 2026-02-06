#ifndef G4EVENTDATA_HH
#define G4EVENTDATA_HH

#include "TPCPadHelper.hh"
#include <vector>
#include <TParticle.h>

using namespace std;

struct G4EventData {


  //Geant4
  vector<TParticle> *BEAM = nullptr;
  vector<TParticle> *PRM = nullptr;
  vector<TParticle> *SEC = nullptr;

};

#endif
