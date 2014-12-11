#include <iostream>
#include <fstream>
#include <string>

#include "spectrum.hpp"
#include "input.hpp"
#include "dose.hpp"
extern "C"
{
	// __declspec(dllexport) 
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
	char * win_mat,
	// Illumination
	double * il_,
	// Solid angle
	double * sa,
	double * dose,
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
	char * win_mat,
	// Illumination
	double * il_,
	// Solid angle
	double * sa,
	// Radiation dose
	double * dose,
	int * nout)
{
	Sample sp;
	Illumination il;
	solid_angle omega;
	Detector det;
	
	// Read input file
	readfile(input_file, sp, il, omega, det);
	
	// Calculate spectrum
	Spectrum spec(sp, il, omega, det, true, true, true);
	Dose ds(sp, il);
	
	// Save results to file
	std::ofstream fout;
	fout.open(output_file);
	spec.out(fout);
	fout << std::endl;
	ds.out(fout);
	il.out(fout);
	det.out(fout);
	fout.close();

   	// Return results
	
	// Total
	if (spec.y_vec.size() > nout[0])
		std::cerr << "Error: Output for y_vec out of range! ("
		<< nout[0] << " given, needs "
		<< spec.y_vec.size() << ".)" << std::endl;
	else
		for (auto i : spec.y_vec)
			*(y_vec++) = i;
		// *y_vec = -1;

	// Separate
	if (spec.y_sep.size() > (nout[1]+2)*nout[0])
		std::cerr << "Error: Output for y_sep out of range! ("
		<< (nout[1]+2)*nout[0] << " given, needs "
		<< spec.y_sep.size() << ".)" << std::endl;
	else
		for (auto i : spec.y_sep)
			*(y_sep++) = i;
		// *y_sep = -1;


	// XRF
	if (spec.xrf.Z_vec.size() > nout[1])
		std::cerr << "Error: Output for xrf.Z_vec out of range! ("
		<< nout[1] << " given, needs "
		<< spec.xrf.Z_vec.size() << ".)" << std::endl;
	else
	{
		for (auto i : spec.xrf.Z_vec)
			*(Z_vec++) = i;
		// *Z_vec = -1;

		for (auto i : spec.xrf.row)
			*(row++) = i;
		// *row = -1;
	}

	if (spec.xrf.lines.size() > nout[2])
		std::cerr << "Error: Output for xrf.lines out of range! ("
		<< nout[2] << " given, needs "
		<< spec.xrf.lines.size() << ".)" << std::endl;
	else
	{
		for (auto i : spec.xrf.lines)
			*(lines++) = i;
		// *lines = -1;

		for (auto i : spec.xrf.ev_vec)
			*(xrf_ev++) = i;
		// *xrf_ev = -1;

		for (auto i : spec.xrf.y_vec)
			*(xrf_y++) = i;
		// *xrf_y = -1;
	}

	// Compton
	if (spec.comp.ev_vec.size() > nout[3])
		std::cerr << "Error: Output for comp.ev_vec out of range! ("
		<< nout[3] << " given, needs "
		<< spec.comp.ev_vec.size() << ".)" << std::endl;
	else
	{
		for (auto i : spec.comp.ev_vec)
			*(comp_ev++) = i;
		// *comp_ev = -1;

		for (auto i : spec.comp.y_vec)
			*(comp_y++) = i;
		// *comp_y = -1;
	}

	// nout
	nout[0] = spec.y_vec.size();
	nout[1] = spec.xrf.Z_vec.size();
	nout[2] = spec.xrf.lines.size();
	nout[3] = spec.comp.ev_vec.size();
	nout[4] = ds.dose_vec.size();


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
//	for (int i = 0; i < det.window.material.length(); i++)
//		*(win_mat++) = int(det.window.material.at(i));
//	*(win_mat) = int('\n');

    for (auto i: det.window.material)
        *(win_mat++) = i;

	// Illumination
	*(il_++) = il.ev0;
	*(il_++) = il.psi;
	*(il_++) = il.alpha;
	*(il_++) = il.n_photons;
	*(il_) = il.beam_cross_section;

	// Solid angle
	for (int i = 0; i < 4; i++)
		*(sa++) = omega.angle_range[i];
	*(sa++) = omega.theta_inc;
	*(sa) = omega.beta_inc;

	// Radiation dosage
	for (int i=0; i<ds.dose_vec.size(); i++)
		*(dose+i) = ds.dose_vec[i];
}

