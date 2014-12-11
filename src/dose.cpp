//dose.cpp
#include <cmath>

#include "dose.hpp"
#include "constants.hpp"

Dose::Dose(const Sample & sp, const Illumination & il) : n_photons(_n_photons), beam_cross_section(_beam_cross_section), dose_vec(_dose_vec), total_dose(_total_dose)
{
    _n_photons = il.n_photons;
    _beam_cross_section = il.beam_cross_section;
    _total_dose = 0;
    double joules_absorbed, matrix_mass;
    for (auto ml : sp.layer_vec)
    {
        if (ml.thickness > 0)
        {
			joules_absorbed = n_photons*(1 - std::exp(-ml.mac_tot(il.ev0)*ml.thickness))*il.ev0*eV_in_Joules;
			matrix_mass = ml.thickness*beam_cross_section*ml.density*1000/std::sin(il.psi);
			_dose_vec.push_back(joules_absorbed/matrix_mass);
			_total_dose += dose_vec.back();
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
	ost << "# # # # # Radiation dosage: # # # # #" << std::endl;
	ost << "For " << n_photons << " photon(s) incident, "
	    << "beam cross-sectional area = " << beam_cross_section << " cm^2 "
	    << "(d = " << std::sqrt(beam_cross_section/Pi)*2e7 << " nm):" << std::endl;
	for (int i=0; i<dose_vec.size(); i++)
	    ost << "\tLayer #" << i+1 << ":\t" << dose_vec[i] << " Gy" << std::endl;
	ost << "\tTotal:\t" << total_dose << " Gy" << std::endl;
	ost << "# # # # # End of radiation dosage # # # # #" << std::endl;
}
