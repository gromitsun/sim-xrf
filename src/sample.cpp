//sample.cpp
#include "xraylib.h"
#include "cs.hpp"
#include "constants.hpp"
#include "sample.hpp"


//Compound
Compound::Compound() : Z_vec(_Z_vec), p_vec(_p_vec), molecular_weight(_molecular_weight)
{

}

Compound::Compound(const std::vector<int> & Z_vec_, const std::vector<double> & p_vec_) : Z_vec(_Z_vec), p_vec(_p_vec), molecular_weight(_molecular_weight)
{
	_Z_vec = Z_vec_;
	_p_vec = p_vec_;
	set_mw();
}

Compound::~Compound()
{

}

Compound & Compound::operator=(const Compound & c)
{
	if (this == &c)
		return *this;
	_Z_vec = c.Z_vec;
	_p_vec = c.p_vec;
	set_mw();
	return *this;
}

void Compound::set_mw()
{
	_molecular_weight = 0;
	for (int i=0; i<Z_vec.size(); i++)
		_molecular_weight += p_vec[i]*AtomicWeight(Z_vec[i]);
}

double Compound::mac_tot(const double & ev) const
{
	double mac = 0;
	for (int i=0; i<Z_vec.size(); i++)
		mac += CS_Total(Z_vec[i], ev/1e3)*p_vec[i]*AtomicWeight(Z_vec[i])/molecular_weight;
	return mac;
}

double Compound::dmac_rayleigh_pol(const double & ev, const double & theta, const double & beta) const
{
	double dmac = 0;
	for (int i=0; i<Z_vec.size(); i++)
		dmac += dcs_rayleigh_pol(theta, beta, ev, Z_vec[i])*p_vec[i]*N_A/molecular_weight;
	return dmac;
}

double Compound::dmac_compton_pol(const double & ev, const double & theta, const double & beta) const
{
	double dmac = 0;
	for (int i=0; i<Z_vec.size(); i++)
		dmac += dcs_compton_pol(theta, beta, ev, Z_vec[i])*p_vec[i]*N_A/molecular_weight;
	return dmac;
}

void Compound::show() const
{
	out(std::cout);
}

void Compound::out(std::ostream & ost) const
{
	ost << "# \tZ = ";
	for (auto Z : Z_vec)
		ost << Z << ", ";
	ost << std::endl;
	
	ost << "# \tp = ";
	for (auto p : p_vec)
		ost << p << ", ";
	ost << std::endl;
	
	ost << "# \tMolecular weight = ";
	ost << molecular_weight;
	ost << std::endl;
}

//Monolayer
Monolayer::Monolayer() : Compound::Compound(), density(_density), thickness(_thickness), layer(_layer)
{

}

// Copy constructor: used when passing through function calls
Monolayer::Monolayer(const Monolayer & ml) : Compound::Compound(), density(_density), thickness(_thickness), layer(_layer)
{
	operator=(ml);
}

Monolayer::Monolayer(const std::vector<int> & _Z_vec, const std::vector<double> & _p_vec, const double & density_, const double & thickness_, const double & layer_) : Compound::Compound(_Z_vec, _p_vec), density(_density), thickness(_thickness), layer(_layer)
{
	_density = density_;
	_thickness = thickness_;
	_layer = layer_;
}

Monolayer::~Monolayer()
{
	
}

// overloaded assignment operator
Monolayer & Monolayer::operator=(const Monolayer & ml)
{
	if (this == &ml)
		return *this;
	Compound::operator=(ml);
	_density = ml.density;
	_thickness = ml.thickness;
	_layer = ml.layer;
	return *this;
}


void Monolayer::set_layer(int i)
{
	_layer = i;
}

void Monolayer::show() const
{
	out(std::cout);
}

void Monolayer::out(std::ostream & ost) const
{
	ost << "# Parameters for layer #" << layer << ":" << std::endl;
	Compound::out(ost);
	ost << "# \tDensity = " << density << std::endl;
	ost << "# \tThickness = " << thickness << std::endl;
}

//Sample
Sample::Sample() : layer_vec(_layer_vec), nlayers(_nlayers)
{
	update();
}

Sample::~Sample()
{

}

Sample & Sample::operator=(const Sample & s)
{
	if (this == & s)
		return *this;
	_layer_vec = s.layer_vec;
	_nlayers = s.nlayers;
	return *this;
}

void Sample::update()
{
	for (int i = 0; i < _layer_vec.size(); i++)
		_layer_vec[i].set_layer(i);	
	_nlayers = _layer_vec.size();
}

void Sample::add_layer(const Monolayer & monolayer_)
{
	_layer_vec.push_back(monolayer_);
	if (monolayer_.layer == 0)
		_layer_vec.back().set_layer(nlayers);
	_nlayers++;
}

void Sample::show() const
{
	out(std::cout);
}

void Sample::out(std::ostream & ost) const
{
	ost << "# ======================================== #" << std::endl;
	ost << "# # # # # Sample parameters: # # # # #" << std::endl;
	ost << "# Total number of layers = " << nlayers << std::endl;
	for (auto m : layer_vec)
		m.out(ost);
	ost << "# # # # # End of sample parameters # # # # #" << std::endl;
	ost << "# ======================================== #" << std::endl;
}