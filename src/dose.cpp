//dose.cpp
#include <cmath>

#include "dose.hpp"
#include "constants.hpp"
#include "math.hpp"

Dose::Dose(Illumination il, Sample sp, double n_photons_, double beam_cross_section_) : n_photons(_n_photons), beam_cross_section(_beam_cross_section), dose_vec(_dose_vec)
{
    _n_photons = n_photons_;
    _beam_cross_section = beam_cross_section_;
    double joules_absorbed, matrix_mass;
    for (auto ml : sp.layer_vec)
    {
        joules_absorbed = n_photons*(1 - std::exp(-ml.mac_tot(il.ev0)*ml.thickness))*il.ev0*eV_in_Joules;
        matrix_mass = ml.thickness*Pi*sq(beam_cross_section)*ml.density*1000/std::sin(il.psi);
        _dose_vec.push_back(joules_absorbed/matrix_mass);
    }
}


Dose::~Dose()
{

}

Dose & Dose::operator=(const Dose & d)
{
	if (this != &d)
	{
		_n_photons = d.n_photons;
		_beam_cross_section = d.beam_cross_section;
		_dose_vec = d.dose_vec;
	}
	return *this;
}

void Dose::show() const
{
	out(std::cout);
}

void Dose::out(std::ostream & ost) const
{
	ost << "# # # # # Radiation dosage: # # # # #*" << std::endl;
	ost << "For " << n_photons << " photons incident, "
	    << "beam cross-sectional area = " << beam_cross_section << "cm^2." << std::endl;
	for (int i=0; i<dose_vec.size(); i++)
	    ost << "Layer #" << i+1 << ":\t" << dose_vec[i] << " Gy" << std::endl;
	ost << "# # # # # End of radiation dosage # # # # #*" << std::endl;
}
