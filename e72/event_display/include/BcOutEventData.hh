#ifndef BCOUTEVENTDATA_HH
#define BCOUTEVENTDATA_HH

#include <vector>

using namespace std;

struct BcOutEventData {

  int ntrack;
  
  vector<double> *x0 = nullptr;
  vector<double> *y0 = nullptr;
  vector<double> *u0 = nullptr;
  vector<double> *v0 = nullptr;
  

};

#endif
