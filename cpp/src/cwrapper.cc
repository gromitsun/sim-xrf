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
	// Separate
	double * y_sep,
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
	double * det_,
	int * n_channels,
	int * win_mat,
	// Illumination
	double * il_,
	// Solid angle
	double * sa,
	int * nout);
}


void sim(char * input_file, 
	char * output_file, 
	// Total
	double * y_vec,
	// Separate
	double * y_sep,
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
	double * det_,
	int * n_channels,
	int * win_mat,
	// Illumination
	double * il_,
	// Solid angle
	double * sa,
	int * nout)
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
	fout.close();
	
	// Return results
	
	// Total
	for (auto i : spec.y_vec)
		*(y_vec++) = i;
	// *y_vec = -1;
	*(nout++) = spec.y_vec.size();
	
	// Separate
	for (auto i : spec.y_sep)
		*(y_sep++) = i;
	// *y_sep = -1;
	
	
	// XRF
	for (auto i : spec.xrf.Z_vec)
		*(Z_vec++) = i;
	// *Z_vec = -1;
	*(nout++) = spec.xrf.Z_vec.size();
	
	for (auto i : spec.xrf.row)
		*(row++) = i;
	*row = -1;
	
	for (auto i : spec.xrf.lines)
		*(lines++) = i;
	// *lines = -1;
	*(nout++) = spec.xrf.lines.size();
	
	for (auto i : spec.xrf.ev_vec)
		*(xrf_ev++) = i;
	// *xrf_ev = -1;
	
	for (auto i : spec.xrf.y_vec)
		*(xrf_y++) = i;
	// *xrf_y = -1;

	// Compton
	for (auto i : spec.comp.ev_vec)
		*(comp_ev++) = i;
	// *comp_ev = -1;
	*(nout++) = spec.comp.ev_vec.size();
	
	for (auto i : spec.comp.y_vec)
		*(comp_y++) = i;
	// *comp_y = -1;

	// Rayleigh
	*ray_y = spec.ray.y;
	
	// Detector
	*(det_++) = det.channel.ev_offset;
	*(det_++) = det.channel.ev_gain;
	*n_channels = det.channel.n_channels;
	*(det_++) = det.response.noise;
	*(det_++) = det.response.fano;
	*(det_++) = det.response.gamma;
	*(det_++) = det.response.fs;
	*(det_++) = det.response.ft;
	*(det_++) = det.window.thickness;
	*det_ = det.window.density;
	for (int i = 0; i < det.window.material.length(); i++)
		*(win_mat++) = int(det.window.material.at(i));
	*(win_mat) = int('\n');
	
	// Illumination
	*(il_++) = il.ev0;
	*(il_++) = il.psi;
	*(il_) = il.alpha;
	
	// Solid angle
	for (int i = 0; i < 4; i++)
		*(sa++) = omega.angle_range[i];
	*(sa++) = omega.theta_inc;
	*(sa) = omega.beta_inc;
	
}

