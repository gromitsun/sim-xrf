#ifndef CS_HPP
#define CS_HPP

double sq(double x);

double ev2nm(double ev); //Convert photon energy (in ev) to wavelength (in nm).

double theta2x(double theta_rad, double ev); //theta in rad -> x in cm^-1


/*differential cross section for single electron*/

double thomson_pol(double theta_rad, double beta_rad);

double thomson_unpol(double theta_rad);

double dcs_rayleigh_pol(double theta_rad, double beta_rad, double ev, int Z);
	
double dcs_rayleigh_unpol(double theta_rad, double ev, int Z);

	
/* Compton scattering cross sections */
	
/* compton dcs of single electron */
//dcs

double ev_scattered(double theta_rad, double ev); //return the energy in ev of the Compton scattered photons.

double klein_nishina_unpol(double theta_rad, double ev); //Differential Klein-Nishina cross section in cm^2 for unpolarized radiation.

double klein_nishina_pol(double theta_rad, double beta_rad, double ev); //Differential Klein-Nishina cross section in cm^2 for polarized radiation.

//total cs
double klein_nishina_total_col(double ev); // Total Klein-Nishina collision cross section in cm^2.

double klein_nishina_total_sca(double ev); //Total Klein-Nishina scattering cross section in cm^2.


//compton dcs of atom
double dcs_compton_pol(double theta_rad, double beta_rad, double ev, int Z); //Differential Compton cross section in cm^2 for polarized radiation.

double dcs_compton_unpol(double theta_rad, double ev, int Z); //Differential Compton cross section in cm^2 for unpolarized radiation.

#endif