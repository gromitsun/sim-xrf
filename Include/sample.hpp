//sample.hpp
#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <vector>
#include <iostream>

class Compound
{
private:
	std::vector<int> _Z_vec;
	std::vector<double> _p_vec;
	double _molecular_weight;
	void set_mw();
public:
	Compound();
	Compound(const std::vector<int> & Z_vec_, const std::vector<double> & p_vec_);
	~Compound();
	Compound & operator=(const Compound & c);
	const std::vector<int> & Z_vec;
	const std::vector<double> & p_vec;
	const double & molecular_weight;
	double mac_tot(const double & ev) const;
	double dmac_rayleigh_pol(const double & ev, const double & theta, const double & beta) const;
	double dmac_compton_pol(const double & ev, const double & theta, const double & beta) const;
	void show() const;
	void out(std::ostream & ost) const;
};

class Monolayer : public Compound
{
private:
	double _density;
	double _thickness;
	double _layer;
public:
	Monolayer();
	Monolayer(const Monolayer & ml);
	Monolayer(const std::vector<int> & _Z_vec, const std::vector<double> & _p_vec, const double & density_, const double & thickness_, const double & layer_ = 0);
	~Monolayer();
	Monolayer & operator=(const Monolayer & ml);
	const double & density;
	const double & thickness;
	const double & layer;
	void set_layer(int i);
	void show() const;
	void out(std::ostream & ost) const;
};

class Sample
{
private:
	std::vector<Monolayer> _layer_vec;
	int _nlayers;
	void update();
public:
	Sample();//const std::vector<Monolayer> & layer_vec_);
	~Sample();
	Sample & operator=(const Sample & s);
	const std::vector<Monolayer> & layer_vec;
	const int & nlayers;
	void add_layer(const Monolayer & monolayer_);
	void show() const;
	void out(std::ostream & ost) const;
};

#endif