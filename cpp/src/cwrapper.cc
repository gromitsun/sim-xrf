#include <iostream>
#include <fstream>
#include <string>

#include "spectrum.hpp"
#include "input.hpp"
extern "C"
{
	__declspec(dllexport) void sim(char * input_file, 
	char * output_file, 
	// Total
	double * y_vec,
	// XRF
	int * Z_vec,
	int * row,
	int * lines,
	double * xrf_ev,
	double * xrf_y,
	// Compton
	double * comp_ev,
	double * comp_y,
	// Rayleigh
	double * ray_y,
	// Detector
	double * ev_offset,
	double * ev_gain,
	int * n_channels);
}


void sim(char * input_file, 
	char * output_file, 
	// Total
	double * y_vec,
	// XRF
	int * Z_vec,
	int * row,
	int * lines,
	double * xrf_ev,
	double * xrf_y,
	// Compton
	double * comp_ev,
	double * comp_y,
	// Rayleigh
	double * ray_y,
	// Detector
	double * ev_offset,
	double * ev_gain,
	int * n_channels)
{
	Sample sp;
	Illumination il;
	solid_angle omega;
	Detector det;
	
	// Read input file
	readfile(input_file, sp, il, omega, det);
	
	// Calculate spectrum
	Spectrum spec(sp, il, omega, det);
	
	// Save results to file
	std::ofstream fout;
	fout.open(output_file);
	spec.out(fout);
	
	// Return results
	
	// Total
	for (auto i : spec.y_vec)
		*(y_vec++) = i;
	*y_vec = -1;
	
	// XRF
	for (auto i : spec.xrf.Z_vec)
		*(Z_vec++) = i;
	*Z_vec = -1;
	
	for (auto i : spec.xrf.row)
		*(row++) = i;
	*row = -1;
	
	for (auto i : spec.xrf.lines)
		*(lines++) = i;
	*lines = -1;
	
	for (auto i : spec.xrf.ev_vec)
		*(xrf_ev++) = i;
	*xrf_ev = -1;
	
	for (auto i : spec.xrf.y_vec)
		*(xrf_y++) = i;
	*xrf_y = -1;

	// Compton
	for (auto i : spec.comp.ev_vec)
		*(comp_ev++) = i;
	*comp_ev = -1;
	
	for (auto i : spec.comp.y_vec)
		*(comp_y++) = i;
	*comp_y = -1;

	// Rayleigh
	*ray_y = spec.ray.y;
	
	// Detector
	*ev_offset = det.channel.ev_offset;
	*ev_gain = det.channel.ev_gain;
	*n_channels = det.channel.n_channels;
}

