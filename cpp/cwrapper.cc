#include <math.h>
#include "snip.hpp"
extern "C"
{
	double FWHM_c(double x, double noise = 100, double fano = 0.114)
	{
		return FWHM(x, noise, fano);
	}
	
	
	double energy_to_channel_c(double energy, double offset = 0, double gain = 10)
	{
		return energy_to_channel(energy, offset, gain);
	}
	
	void lsdf_c(double *E, double *y, int lenE, 
					double (*FWHM)(double),
					double f = 1.5,
					double A = 75,
					double M = 10,
					double r = 1.3)
	{
		return lsdf(E, y, lenE, FWHM, f, A, M, r);
	}
		
	void snip_c(double *E, double *y, int lenE,
			double (*FWHM)(double), 
			double offset = 0, double gain = 10,
			int loops = 24, int end_loops = 8, double factor = 2, double reduce_factor = sqrt(2))
	{
		return snip(E, y, lenE, FWHM, offset, gain, loops, end_loops, factor, reduce_factor);
	}
}
