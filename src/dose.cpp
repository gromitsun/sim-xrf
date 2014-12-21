//dose.cpp
#include <cmath>

#include "dose.hpp"
#include "constants.hpp"

Dose::Dose(const Sample & sp, const Illumination & il) : n_photons(_n_photons), beam_cross_section(_beam_cross_section), dose_vec(_dose_vec), total_dose(_total_dose)
{
    _n_photons = il.n_photons;
	_beam_cross_section = il.beam_cross_section;
    _total_dose = 0;

	double _n_photons_calc, _beam_cross_section_calc;
    if ( n_photons > 0 && beam_cross_section > 0)
    {
		_n_photons_calc = n_photons;
		_beam_cross_section_calc = beam_cross_section;
    }
    else
    {
    	_n_photons_calc = 1;
		_beam_cross_section_calc = 1;
    }

    double joules_absorbed, matrix_mass, n_photons_in=_n_photons_calc, n_photons_out;
    for (auto ml : sp.layer_vec)
    {
        if (ml.thickness > 0)
        {
        	n_photons_out = n_photons_in*std::exp(-ml.mac_tot(il.ev0)*ml.thickness);
			joules_absorbed = (n_photons_in - n_photons_out)*il.ev0*eV_in_Joules;
			matrix_mass = ml.thickness*_beam_cross_section_calc*ml.density*1000/std::sin(il.psi);
			_dose_vec.push_back(joules_absorbed/matrix_mass);
			_total_dose += dose_vec.back();
			n_photons_in = n_photons_out;
        }
        else
        	_dose_vec.push_back(-1);
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
		_total_dose = d.total_dose;
	}
	return *this;
}

void Dose::show() const
{
	out(std::cout);
}

void Dose::out(std::ostream & ost) const
{
	ost << "# ======================================== #" << std::endl;
	ost << "# # # # # Radiation dosage: # # # # #" << std::endl;
	if ( n_photons > 0 && beam_cross_section > 0)
    {
		ost << "# For " << n_photons << " photon(s) incident, "
			<< "beam cross-sectional area = " << beam_cross_section << " cm^2 "
			<< "(d = " << std::sqrt(beam_cross_section/Pi)*2e7 << " nm):" << std::endl;
	}
	else
		ost << "# Dose per incident photon number density (1 photon/cm^2)" << std::endl;
	for (int i=0; i<dose_vec.size(); i++)
	    ost << "# \tLayer #" << i+1 << ":\t" << dose_vec[i] << " Gy" << std::endl;
	ost << "# \tTotal:\t" << total_dose << " Gy" << std::endl;
	ost << "# # # # # End of radiation dosage # # # # #" << std::endl;
	ost << "# ======================================== #" << std::endl;
}
