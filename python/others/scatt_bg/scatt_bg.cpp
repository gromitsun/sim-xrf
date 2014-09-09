#include <iostream>
#include <cmath>
#include "xraylib.h"

using namespace std;

#ifndef PI
#define PI  3.1415926535897932384626433832795
#endif

extern "C"
{
	double subtend(double theta_rad_start, double theta_rad_end, double beta_rad_start, double beta_rad_end)
	{
		return (cos(theta_rad_start)-cos(theta_rad_end))*(beta_rad_end-beta_rad_start);
	}
	
	double radians(double deg)
	{
		return deg/180.*PI;
	}

	void scatt_bg(double kev, double *out, int Z_max, int theta_max)
	{
		double x_rad, omega, domega, temp, prev,foo;
		for (int Z = 1; Z <= Z_max; Z++)
		{
			cout << "Z = " << Z <<endl;
			for (int x = 1; x <= theta_max; x++)
			{
				x_rad = radians(x);
				omega = subtend(PI/2-x_rad, PI/2+x_rad, -x_rad, x_rad);
				domega = sin(PI/2-x_rad)*PI*PI/180/180;
				if (x-1 != 0)
					out[Z*theta_max+x-1] = prev;
				temp = 0;
				for (int i = -x; i <= x; i++)
				{
					if (PI/2-x_rad > 0)
					{
						temp += DCSP_Rayl(Z,kev,PI/2-x_rad,i*PI/180);
						temp += DCSP_Rayl(Z,kev,PI/2+x_rad,i*PI/180);
						temp += DCSP_Compt(Z,kev,PI/2-x_rad,i*PI/180);
						temp += DCSP_Compt(Z,kev,PI/2+x_rad,i*PI/180);
					}
				}
				out[(Z-1)*theta_max+x-1] += temp*domega;
				
				for (int i = -x+1; i <= x-1; i++)
				{
					temp = 0;
					domega = sin(PI/2+i*PI/180)*PI*PI/180/180;
					temp += DCSP_Rayl(Z,kev,PI/2+i*PI/180,-x_rad);
					temp += DCSP_Rayl(Z,kev,PI/2+i*PI/180,x_rad);
					temp += DCSP_Compt(Z,kev,PI/2+i*PI/180,-x_rad);
					temp += DCSP_Compt(Z,kev,PI/2+i*PI/180,x_rad);
					out[(Z-1)*theta_max+x-1] += temp*domega;
				}
				prev = out[(Z-1)*theta_max+x-1];
				out[(Z-1)*theta_max+x-1] /= omega*omega;
			}
		}
	}
}