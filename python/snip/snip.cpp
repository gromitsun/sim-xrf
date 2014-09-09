#include <iostream>
#include <cmath>

using namespace std;

double FWHM(double x, double noise = 100, double fano = 0.114)
{
	return 2.3548*sqrt((noise/2.3548)*(noise/2.3548)+3.58*fano*x);
}

double energy_to_channel(double energy, double offset = 0, double gain = 10)
{
	return (energy - offset)/gain;
	
}


// Low statistics digital filter
double reduce(int x, int length_start, double *y,
				double f = 1.5,
				double A = 75,
				double M = 10,
				double r = 1.3)
{
	for (int i=0; i < length_start; i++)
	{
		int length = length_start - i;
		if (x < length)
			return -1;
		
		double L = 0, R = 0;
		for (int i=0; i < length; i++)
		{
			L += y[x-length+i];
			R += y[x+1+i];
		}
		double S = y[x] + L + R;
		double slope = (R + 1.)/(L + 1.);
		if ((S < M) || (S < A*sqrt(y[x])) || ((slope >= 1./r) && (slope <= r)))
			return S/(2.*length+1);
	}
	cout << "Not found for x = " << x << '!' << endl;
	return y[x];
}


void lsdf(double *E, double *y, int lenE, 
			double (*FWHM)(double),
			double f = 1.5,
			double A = 75,
			double M = 10,
			double r = 1.3)
{
	for (int x = 0; x < lenE; x++)
	{
		int len_0 = (int) energy_to_channel(f*FWHM(E[x]), E[0], E[1]-E[0]);
		double temp = reduce(x, len_0, y, f, A, M, r);
		if (temp != -1)
			y[x] = temp;
	}
}

// Peak-clipping

void snip(double *E, double *y, int lenE,
		double (*FWHM)(double), 
		double offset = 0, double gain = 10,
		int loops = 24, int end_loops = 8, double factor = 2, double reduce_factor = sqrt(2))
{	
	double z[lenE],z_out[lenE],*p=z,*q=z_out;
	for (int i = 0; i < lenE; i++)
	{
		z[i] = log(log(y[i]+1)+1);
		z_out[i] = z[i];
	}
	int w[lenE];
	for (int i = 0; i < loops - end_loops; i++)
	{
		for (int x=0; x<lenE; x++)
		{
			#ifndef width
			w[x] = (int) energy_to_channel(factor*FWHM(E[x]), offset, gain);
			#endif
			if ((w[x] <= x) && ((w[x]+x)<lenE))
				q[x] = min(p[x],(p[x+w[x]]+p[x-w[x]])/2.);
		}
		#define width
		swap(p,q);
	}
	
	for (int i = 0; i < end_loops; i++)
	{
		factor /= 1.*reduce_factor;
		for (int x=0; x<lenE; x++)
		{
			w[x] = (int) energy_to_channel(factor*FWHM(E[x]), offset, gain);
			if ((w[x] <= x) && ((w[x]+x)<lenE))
				q[x] = min(p[x],(p[x+w[x]]+p[x-w[x]])/2.);
		}		
		swap(p,q);
	}
	
	//inverse G
	for (int i = 0; i < lenE; i++)
		y[i] = exp(exp(p[i])-1)-1;
		
}