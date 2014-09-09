/* xrf.hpp */
#ifndef XRF_HPP
#define XRF_HPP
#include <vector>

using namespace std;

#define K_LINES -29
#define	L_LINES -113
#define M_LINES -219
#define N_LINES -324
#define O_LINES -374
#define P_LINES -383

vector<int> mac_xrf(double ev0, int Z, vector<double> & ev_vec, vector<double> & y_vec, double weight = 1.0, int line_end = P_LINES);

#endif


