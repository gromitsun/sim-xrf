/* xrf.cpp */
#include <iostream>
#include <vector>
#include "xraylib.h"

using namespace std;

#define K_LINES -29
#define	L_LINES -113
#define M_LINES -219
#define N_LINES -324
#define O_LINES -374
#define P_LINES -383

vector<int> mac_xrf(double ev0, int Z, vector<double> & ev_vec, vector<double> & y_vec, double weight = 1.0, int line_end = P_LINES)
{
	vector<int> lines;
	int line_start;
	if (ev0 >= 1e3*EdgeEnergy(Z, K_SHELL))
		line_start = -1;
	else if (line_end <= L_LINES && ev0 >= 1e3*EdgeEnergy(Z, L3_SHELL))
		line_start = K_LINES-1;
	else if (line_end <= M_LINES && ev0 >= 1e3*EdgeEnergy(Z, M5_SHELL))
		line_start = L_LINES-1;
	else if (line_end <= N_LINES && ev0 >= 1e3*EdgeEnergy(Z, N7_SHELL))
		line_start = M_LINES-1;
	else if (line_end <= O_LINES && ev0 >= 1e3*EdgeEnergy(Z, O7_SHELL))
		line_start = N_LINES-1;	
	else if (line_end <= P_LINES && ev0 >= 1e3*EdgeEnergy(Z, P5_SHELL))
		line_start = N_LINES-1;
	else 
	{
		cout << "XRF lines not found for Z = " << Z << " @ " << ev0 << " eV!" << endl;
		return lines;
	}
	
	for (int line = line_start; line >= line_end; line--)
	{
		double ev = LineEnergy(Z, line)*1e3;
		if (ev > 1e-6)
		{
			double y = CS_FluorLine_Kissel_Cascade(Z, line, ev0/1e3);
			if (y > 1e-50)
			{
				ev_vec.push_back(ev);
				y_vec.push_back(y*weight);
				lines.push_back(line);
			}
		}
	}
	return lines;
}




