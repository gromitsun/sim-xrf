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
	double _total_dose;
	std::vector<double> _dose_vec;
public:
	Dose(const Sample & sp, const Illumination & il);
	~Dose();
	Dose & operator=(const Dose & d);
	const double & n_photons;
	const double & beam_cross_section;
	const double & total_dose;
	const std::vector<double> & dose_vec;
	void show() const;
	void out(std::ostream & ost) const;
};


#endif