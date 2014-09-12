//geometry.hpp
#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>
#include <ostream>
#include "constants.hpp"
#include "sample.hpp"

void spherical2cartesian(const double & r, const double & theta, const double & beta, double & x, double & y, double & z); //Convert spherical coordinates (r,theta_rad,beta_rad) to Cartesian coordinates (x,y,z).
double arctan2(const double & x);
void rotationx(const double & rad, double & x, double & y, double & z);
void rotationy(const double & rad, double & x, double & y, double & z);

class Illumination
{
private:
	double _ev0;
	double _psi;
	double _alpha;
public:
	Illumination(double ev0_ = 1e4, double psi_ = Pi/4, double alpha_ = 0);
	~Illumination();
	Illumination & operator=(const Illumination & il);
	const double & ev0;
	const double & psi;
	const double & alpha;
	double psi_prime(const double & theta, const double & beta) const;
	void show() const;
	void out(std::ostream & ost) const;
};

class solid_angle
{
private:
	double angle_range[4];
	std::vector<double> _theta;
	std::vector<double> _beta;
	double _theta_inc;
	double _beta_inc;
	double _subtend;
	void update();
public:
	solid_angle();
	solid_angle(const double *ar, double th_inc = Pi/180, double be_inc = Pi/180);
	~solid_angle();
	solid_angle & operator=(const solid_angle & sa);
	// const std::vector<double> & get_theta(double th_inc = -1);
	// const std::vector<double> & get_beta(double be_inc = -1);
	const std::vector<double> & theta;
	const std::vector<double> & beta;
	const double & subtend;
	double domega(const double & theta, const double & beta) const;
	void show() const;
	void out(std::ostream & ost) const;
};

double atten_mono(const double & ev0, 
		const double & ev, 
		const double & psi, 
		const double & psiprime, 
		const Monolayer & c);
		
double atten_refl(const double & ev0, 
		const double & ev, 
		const double & psi, 
		const double & psiprime, 
		const Monolayer & ml);	
		
double atten_trans_in(const double & ev0, 
		const double & psi, 
		const Monolayer & ml);	
		
double atten_trans_out(const double & ev, 
		const double & psiprime, 
		const Monolayer & ml);

#endif