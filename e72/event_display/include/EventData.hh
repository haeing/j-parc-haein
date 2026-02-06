#ifndef EVENTDATA_HH
#define EVENTDATA_HH

#include "TPCPadHelper.hh"
#include <vector>

using namespace std;

struct EventData {

  int nhittpc;
  int ntTpc;
  vector<int> *nhtrack = nullptr;
  vector<int> *trackid = nullptr;

  vector<double>* x0Tpc = nullptr;
  vector<double>* y0Tpc = nullptr;
  vector<double>* u0Tpc = nullptr;
  vector<double>* v0Tpc = nullptr;

  vector<double> *helix_cx = nullptr;
  vector<double> *helix_cy = nullptr;
  vector<double> *helix_z0 = nullptr;
  vector<double> *helix_r = nullptr;
  vector<double> *helix_dz = nullptr;
  vector<double> *dz_factor = nullptr;

  vector<vector<double>> *helix_t = nullptr;

  vector<vector<double>> *hitpos_x = nullptr;
  vector<vector<double>> *hitpos_y = nullptr;
  vector<vector<double>> *hitpos_z = nullptr;

  
};

#endif
