//dose.hpp
#ifndef DOSE_HPP
#define DOSE_HPP

#include <vector>
#include <iostream>

#include "geometry.hpp"
#include "sample.hpp"

class Dose
{
private:
	double _n_photons;
	double _beam_cross_section;
	std::vector<double> _dose_vec;
public:
	Dose(Illumination il, Sample sp, double n_photons_=1, double beam_cross_section_=1);
	~Dose();
	Dose & operator=(const Dose & d);
	const double & n_photons;
	const double & beam_cross_section;
	const std::vector<double> & dose_vec;
	void show() const;
	void out(std::ostream & ost) const;
};


#endif