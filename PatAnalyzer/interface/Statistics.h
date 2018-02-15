#ifndef STATS_H
#define STATS_H
#include "TMath.h"
using namespace std;

namespace stats{
  double LowBound(int n , int N);
  double LowBound(int n);
  double HighBound(int n, int N);
  double HighBound(int n );

}

#endif
