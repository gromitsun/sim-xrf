//geometry.cpp
#include <iostream>
#include <cmath>

#include "geometry.hpp"

using std::vector;


//Illumination
Illumination::Illumination(double ev0_, double psi_, double alpha_, double n_photons_, double beam_cross_section_) : ev0(_ev0), psi(_psi), alpha(_alpha), n_photons(_n_photons), beam_cross_section(_beam_cross_section)
{
	_ev0 = ev0_;
	_psi = psi_;
	_alpha = alpha_;
	_n_photons = n_photons_;
	_beam_cross_section = beam_cross_section_;
}
Illumination::~Illumination()
{

}

Illumination & Illumination::operator=(const Illumination & il)
{
	if (this != &il)
	{
		_ev0 = il.ev0;
		_psi = il.psi;
		_alpha = il.alpha;
		_n_photons = il.n_photons;
		_beam_cross_section = il.beam_cross_section;
	}
	return *this;
}

void Illumination::show() const
{
	out(std::cout);
}

void Illumination::out(std::ostream & ost) const
{
	ost << "# ======================================== #" << std::endl;
	ost << "# # # # # Illumination: # # # # #" << std::endl;
	ost << "# ev0 = " << ev0 << std::endl;
	ost << "# psi = " << psi << std::endl;
	ost << "# alpha = " << alpha << std::endl;
	ost << "# Number of incident photons = " << n_photons << std::endl;
	ost << "# Beam cross sectional area = " << beam_cross_section << " cm^2" << std::endl;
	ost << "# # # # # End of illumination # # # # #" << std::endl;
	ost << "# ======================================== #" << std::endl;
	
}


void spherical2cartesian(const double & r, const double & theta, const double & beta, double & x, double & y, double & z) //Convert spherical coordinates (r,theta,beta) to Cartesian coordinates (x,y,z).
{
	using std::sin;
	using std::cos;
	x = r*sin(theta)*cos(beta);
	y = r*sin(theta)*sin(beta);
	z = r*cos(theta);
}

double arctan2(const double & x) //like arctan, but the returned value is in [0,np.pi]
{	
	using std::atan;
	double rad = atan(x);
	if (rad<0)
		rad += Pi;
	return rad;
}

void rotationy(const double & rad, double & x, double & y, double & z)
//rotate the axes about original y axis, angle = rad.
{
	if (rad)
	{
		using std::cos;
		using std::sin;
		double x1;
		x1 = cos(rad)*x - sin(rad)*z;
		z = sin(rad)*x + cos(rad)*z;
		x = x1;
	}
}

void rotationx(const double & rad, double & x, double & y, double & z)
//rotate the axes about original x axis, angle = rad.
{
	if (rad)
	{
		using std::cos;
		using std::sin;
		double y1;
		y1 = cos(rad)*y + sin(rad)*z;
		z = -sin(rad)*y + cos(rad)*z;
		y = y1;
	}
}
double Illumination::psi_prime(const double & theta, const double & beta) const
{
	using std::sqrt;
	double theta0 = Pi/2.+psi;
	double x, y, z;
	spherical2cartesian(1, theta, beta, x, y, z);
	rotationy(theta0, x, y, z);
	rotationx(alpha, x, y, z);
	double theta_prime = arctan2(sqrt(x*x+y*y)/(1.*z));
	return Pi/2.-theta_prime;
}
	
// solid_angle
solid_angle::solid_angle():theta(_theta), beta(_beta), subtend(_subtend), theta_inc(_theta_inc), beta_inc(_beta_inc)//, angle_range(_angle_range)
{
	_angle_range[0] = 0;
	_angle_range[1] = Pi;
	_angle_range[2] = 0;
	_angle_range[3] = Pi/2;
	_theta_inc = Pi/180;
	_beta_inc = Pi/180;
	update();
}
solid_angle::solid_angle(const double *ar, double th_inc, double be_inc):theta(_theta), beta(_beta), subtend(_subtend), theta_inc(_theta_inc), beta_inc(_beta_inc)//, angle_range(_angle_range)
{
	for (int i=0; i<4; i++)
		_angle_range[i] = ar[i];
	_theta_inc = th_inc;
	_beta_inc = be_inc;
	update();
}
solid_angle::~solid_angle()
{

}

solid_angle & solid_angle::operator=(const solid_angle & sa)
{
	if (this == &sa)
		return *this;
	for (int i=0; i<4; i++)
		_angle_range[i] = sa.angle_range[i];
	_theta_inc = sa._theta_inc;
	_beta_inc = sa._beta_inc;
	update();
	return *this;
}

double solid_angle::domega(const double & theta) const
{
	return std::sin(theta)*_theta_inc*_beta_inc;
}

void solid_angle::update()
{
	_theta.clear();
	for (double x = angle_range[0]; x <= angle_range[1]; x += _theta_inc)
		_theta.push_back(x);
	_beta.clear();
	for (double x = angle_range[2]; x <= angle_range[3]; x += _beta_inc)
		_beta.push_back(x);
	_subtend = (std::cos(angle_range[0])-std::cos(angle_range[1]))*(angle_range[3]-angle_range[2]);
}

// const vector<double> & solid_angle::get_theta(double th_inc)
// {
	// if ((th_inc > 0) && (th_inc != theta_inc))
	// {
		// theta_inc = th_inc;
		// _theta.clear();
		// for (double x = angle_range[0]; x <= angle_range[1]; x += theta_inc)
			// _theta.push_back(x);
	// }
	// else if (theta.empty())
		// for (double x = angle_range[2]; x <= angle_range[3]; x += theta_inc)
			// _theta.push_back(x);
	// return theta;
// }



// const vector<double> & solid_angle::get_beta(double be_inc)
// {
	// if ((be_inc > 0) && (be_inc != beta_inc))
	// {
		// beta_inc = be_inc;
		// _beta.clear();
		// for (double x = angle_range[2]; x <= angle_range[3]; x += beta_inc)
			// _beta.push_back(x);
	// }
	// else if (beta.empty())
		// for (double x = angle_range[2]; x <= angle_range[3]; x += beta_inc)
			// _beta.push_back(x);
	// return beta;
// } 

void solid_angle::show() const
{
	out(std::cout);
}

void solid_angle::out(std::ostream & ost) const
{
	ost << "# ======================================== #" << std::endl;
	ost << "# # # # # Solid angle definition # # # # #" << std::endl;
	ost << "# angle range: ";
	for (int i = 0; i < 4; i++)
		ost << angle_range[i] << ' ';
	ost << std::endl;
	ost << "# theta_inc = " << _theta_inc << std::endl;
	ost << "# beta_inc = " << _beta_inc << std::endl;
	ost << "# theta: ";
	for (vector<double>::const_iterator i = theta.begin(); i < theta.end(); i++)
		ost << *i << ' ';
	ost << std::endl;
	ost << "# beta: ";
	for (vector<double>::const_iterator i = beta.begin(); i < beta.end(); i++)
		ost << *i << ' ';
	ost << std::endl;
	ost << "# Total subtended solid angle: " << subtend << std::endl;
	ost << "# # # # # End of solid angle definition # # # # #" << std::endl;
	ost << "# ======================================== #" << std::endl;
}

double atten_mono(const double & ev0, 
		const double & ev, 
		const double & psi, 
		const double & psiprime, 
		const Monolayer & ml)
{
	if ((psiprime < 1e-6) && (psiprime > -1e-6)) // Check if psiprime == 0
		return 0;
	double t = ml.thickness;
	double rho = ml.density;
	if (t == 0 || rho == 0)
		return 1;
	double mac0 = ml.mac_tot(ev0);
	double mac1 = ml.mac_tot(ev);
	double sp0 = std::sin(psi);
	double sp1 = std::sin(psiprime);
	double temp = (mac0/sp0+mac1/sp1)*rho;
	if (mac0+mac1*sp0/sp1 < 1e-50 && mac0+mac1*sp0/sp1 > -1e-50) // To avoid division by zero
		return 0;
	if (psiprime > 0)
		return (1-std::exp(-(mac0/sp0+mac1/sp1)*rho*t))/(mac0+mac1*sp0/sp1);
	else
		return (std::exp(mac1*rho*t/sp1)-std::exp(-mac0*rho*t/sp0))/(mac0+mac1*sp0/sp1);
}

double atten_refl(const double & ev0, 
		const double & ev, 
		const double & psi, 
		const double & psiprime, 
		const Monolayer & ml)
{
	double t = ml.thickness;
	double rho = ml.density;
	if (t == 0 || rho == 0)
		return 1;
	else if (psi == 0 || psiprime == 0)
	    return 0;
	return std::exp(-(ml.mac_tot(ev0)/std::sin(psi)+ml.mac_tot(ev)/std::sin(psiprime))*ml.density*ml.thickness);
}		
		
double atten_trans_in(const double & ev0, 
		const double & psi, 
		const Monolayer & ml)
{
	if (ml.density == 0 || ml.thickness == 0)
		return 1;
	else if (psi == 0)
	    return 0;
	return std::exp(-ml.mac_tot(ev0)/std::sin(psi)*ml.density*ml.thickness);
}		
		
double atten_trans_out(const double & ev, 
		const double & psiprime, 
		const Monolayer & ml)
{
	if (ml.density == 0 || ml.thickness == 0)
		return 1;
	else if (psiprime == 0)
	    return 0;
	return std::exp(ml.mac_tot(ev)/std::sin(psiprime)*ml.density*ml.thickness);
}