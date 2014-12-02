#ifndef SNIP_HPP
#define SNIP_HPP

double FWHM(double x, double noise, double fano);

double energy_to_channel(double energy, double offset, double gain);

void lsdf(double *E, double *y, int lenE, 
			double (*FWHM)(double),
			double f,
			double A,
			double M,
			double r);
			
void snip(double *E, double *y, int lenE,
		double (*FWHM)(double), 
		double offset, double gain,
		int loops, int end_loops, double factor, double reduce_factor);


#endif
