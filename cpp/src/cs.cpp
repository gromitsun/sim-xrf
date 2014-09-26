#include <cmath>

#include "xraylib.h"

#include "constants.hpp"
#include "math.hpp"

inline double ev2nm(double ev) //Convert photon energy (in ev) to wavelength (in nm).
{
	return 1239.84187/ev;
}

double theta2x(double theta_rad, double ev) //theta in rad -> x in cm^-1
{
	return sin(theta_rad/2.)/(1.e-7*ev2nm(ev));
}


/*differential cross section for single electron*/

double thomson_pol(double theta_rad, double beta_rad)
{
	return sq(r_e)*(sq(cos(beta_rad)*cos(theta_rad))+sq(sin(beta_rad)));
}

double thomson_unpol(double theta_rad)
{
	return sq(r_e)*(1+sq(cos(theta_rad)))/2.;
}

double dcs_rayleigh_pol(double theta_rad, double beta_rad, double ev, int Z)
{
	return sq(FF_Rayl(Z, 1e-8*theta2x(theta_rad, ev)))*thomson_pol(theta_rad, beta_rad);
}
	
	
double dcs_rayleigh_unpol(double theta_rad, double ev, int Z)
{
	return sq(FF_Rayl(Z, 1e-8*theta2x(theta_rad, ev)))*thomson_unpol(theta_rad);
}	

	
/* Compton scattering cross sections */
	
/* compton dcs of single electron */
//dcs

double ev_scattered(double theta_rad, double ev) //return the energy in ev of the Compton scattered photons.
{
	return ev/(1+ev/m_e*(1-cos(theta_rad)));
}

double klein_nishina_unpol(double theta_rad, double ev) //Differential Klein-Nishina cross section in cm^2 for unpolarized radiation.

{
	double ev1 = ev_scattered(theta_rad, ev);
	return sq(r_e*(ev1/ev))/2.*(ev1/ev+ev/ev1-sq(sin(theta_rad)));
}

double klein_nishina_pol(double theta_rad, double beta_rad, double ev) //Differential Klein-Nishina cross section in cm^2 for polarized radiation.

{
	double ev1 = ev_scattered(theta_rad, ev);
	return sq(r_e*(ev1/ev))/2.*(ev1/ev+ev/ev1-2*sq(sin(theta_rad)*cos(beta_rad)));
}

//total cs
double klein_nishina_total_col(double ev) // Total Klein-Nishina collision cross section in cm^2.
{
	double a = ev/m_e;
	return 2*Pi*sq(r_e)*((1+a)/(a*a*a)*(2*a*(1+a)/(1+2*a)-log(1+2*a))+log(1+2*a)/(2*a)-(1+3*a)/sq(1+2*a));
}


double klein_nishina_total_sca(double ev) //Total Klein-Nishina scattering cross section in cm^2.
{
	double a = ev/m_e;
	return Pi*sq(r_e)*(log(1+2*a)/(a*a*a)+2*(1+a)*(2*sq(a)-2*a-1)/sq(a*(1+2*a))+8*sq(a)/(3*sq(1+2*a)*(1+2*a)));
}
//compton dcs of atom
double dcs_compton_pol(double theta_rad, double beta_rad, double ev, int Z) //Differential Compton cross section in cm^2 for polarized radiation.
{
	return klein_nishina_pol(theta_rad, beta_rad, ev)*SF_Compt(Z, theta2x(theta_rad,ev)*1e-8);
}
	

double dcs_compton_unpol(double theta_rad, double ev, int Z) //Differential Compton cross section in cm^2 for unpolarized radiation.
{
	return klein_nishina_unpol(theta_rad, ev)*SF_Compt(Z, theta2x(theta_rad,ev)*1e-8);
}