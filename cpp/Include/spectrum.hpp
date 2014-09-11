//spectrum.hpp
#ifndef SPECTRUM_HPP
#define SPECTRUM_HPP

#include <vector>

#include "sample.hpp"
#include "geometry.hpp"
#include "detector.hpp"

class Xrf
{
private:
	std::vector<double> _y_mat;
	solid_angle _omega;
	std::vector<double> _ev_vec;
	std::vector<double> _y_vec;
	std::vector<int> _lines;
	std::vector<int> _Z_vec;
	std::vector<int> _row;
public:
	Xrf(); //default constructor
	Xrf(double ev0, const Compound & c);
	Xrf(double ev0, const Compound & c, const solid_angle omega);
	Xrf(const Sample & s, const Illumination & il, const solid_angle & omega_);
	// ~Xrf();
	Xrf & operator=(const Xrf & x);
	const std::vector<double> & y_mat;
	const solid_angle & omega;
	const std::vector<double> & ev_vec;
	const std::vector<double> & y_vec;
	const std::vector<int> & lines;
	const std::vector<int> & Z_vec;
	const std::vector<int> & row;
	// Xrf & operator+(const Xrf & x) const;
	void add(const Xrf & x, bool mat_only = false);
	void sum();
	void show() const;
	void out(std::ostream & ost) const;
};

class Rayleigh
{
private:
	std::vector<double> _y_mat;
	solid_angle _omega;
	double _y;
	double _ev0;
public:
	Rayleigh(); //default constructor
	Rayleigh(double ev0_, const Compound & c, const solid_angle & omega_, bool calc_sum = true);
	Rayleigh(const Sample & s, const Illumination & il, const solid_angle & omega_);
	// ~Rayleigh();
	Rayleigh & operator=(const Rayleigh & r);
	const std::vector<double> & y_mat;
	const solid_angle & omega;
	const double & y;
	const double & ev0;
	// Rayleigh & operator+(const Rayleigh & r) const;
	void add(const Rayleigh & x, bool mat_only = false);
	void sum();
	void show() const;
	void out(std::ostream & ost) const;
};

class Compton
{
private:
	double _ev0;
	std::vector<double> _y_mat;
	solid_angle _omega;
	std::vector<double> _ev_vec;
	std::vector<double> _y_vec;
public:
	Compton(); //default constructor
	Compton(double ev0_, const Compound & c, const solid_angle & omega_, bool calc_ev = true, bool calc_sum = true);
	Compton(const Sample & s, const Illumination & il, const solid_angle & omega_);
	// ~Compton();
	Compton & operator=(const Compton & c);
	const std::vector<double> & y_mat;
	const solid_angle & omega;
	const std::vector<double> & ev_vec;
	const std::vector<double> & y_vec;
	// Compton & operator+(const Compton & c) const;
	void add(const Compton & x, bool mat_only = false);
	const double & ev0;
	void sum();
	void ev();
	void sum_ev();
	void show() const;
	void out(std::ostream & ost) const;
};

class Spectrum
{
private:
	std::vector<double> _y_vec;
	// std::vector<int> _Z_vec;
	// std::vector<int> _lines;
	// std::vector<int> _row;
	Xrf _xrf;
	Rayleigh _ray;
	Compton _comp;
	solid_angle _omega;
	Sample _sample;
	Illumination _illumination;
	Detector _detector;
public:
	Spectrum();
	Spectrum(const Sample & s, const Illumination & il, const solid_angle & omega, const Detector & det);
	Spectrum(double ev0, const Compound & c, const solid_angle & omega, const Detector & det);
	// ~Spectrum();
	Spectrum & operator=(const Spectrum & s);	
	
	const std::vector<double> & y_vec;
	const std::vector<int> & Z_vec;
	const std::vector<int> & lines;
	const std::vector<int> & row;
	const Xrf & xrf;
	const Rayleigh & ray;
	const Compton & comp;
	const solid_angle & omega;
	
	const Sample & sample;
	const Illumination & illumination;
	const Detector & detector;
	
	void set_geom(const Sample & s, const Illumination & il, const solid_angle & omega);
	void set_det(const Detector & det);
	
	void calc(bool self_abs = true);
	
	void genspec(bool det_response = true);
	// void genspec_raw(std::vector<double> & ev_vec_, std::vector<double> & y_vec_);	
	void show() const;
	void out(std::ostream & ost) const;
};

	
	
#endif