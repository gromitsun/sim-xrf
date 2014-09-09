/*spectrum.cpp*/

#include "spectrum.hpp"

#include <iostream>
#include <numeric> //std::accumulate
#include <algorithm> //std::find

#include "xraylib.h"
#include "xrf.hpp"
#include "cs.hpp"

Xrf::Xrf() : ev_vec(_ev_vec), y_vec(_y_vec), lines(_lines), Z_vec(_Z_vec), row(_row), y_mat(_y_mat), omega(_omega) {}

Xrf::Xrf(double ev0, const Compound & c)
	: ev_vec(_ev_vec), y_vec(_y_vec), lines(_lines), Z_vec(_Z_vec), row(_row), y_mat(_y_mat), omega(_omega)
{
	std::vector<int> lines_temp;
	std::vector<double> ev_vec_temp, y_vec_temp;
	double weight;
	int Z;
	_row.push_back(0);
	for (int i = 0; i < c.Z_vec.size(); i++)
	{
		Z = c.Z_vec[i];
		weight = c.p_vec[i]*AtomicWeight(Z)/c.molecular_weight;
		lines_temp = mac_xrf(ev0, Z, ev_vec_temp, y_vec_temp, weight);
		if (lines_temp.size() > 0)
		{
			_Z_vec.push_back(Z);
			_lines.insert(_lines.end(),lines_temp.begin(),lines_temp.end());
			_ev_vec.insert(_ev_vec.end(),ev_vec_temp.begin(),ev_vec_temp.end());
			_y_vec.insert(_y_vec.end(),y_vec_temp.begin(),y_vec_temp.end());
			_row.push_back(lines.size());
		}
	}
}

//c++98:
// Xrf::Xrf(double ev0, const Compound & c, const solid_angle omega)
	// : ev_vec(_ev_vec), y_vec(_y_vec), lines(_lines), Z_vec(_Z_vec), row(_row), y_mat(_y_mat), omega(_omega)
// {
	// *this = Xrf(ev0, c);
	// for (int i = 0; i < y_vec.size(); i++)
		// _y_vec[i] *= omega.subtend/(4*Pi);
// }

// c++11:
Xrf::Xrf(double ev0, const Compound & c, const solid_angle omega)
	: Xrf(ev0, c)
{
	for (int i = 0; i < y_vec.size(); i++)
		_y_vec[i] *= omega.subtend/(4*Pi);
}

Xrf::Xrf(const Sample & s, const Illumination & il, const solid_angle & omega_)
	: ev_vec(_ev_vec), y_vec(_y_vec), lines(_lines), Z_vec(_Z_vec), row(_row), y_mat(_y_mat), omega(_omega)
{
	_omega = omega_;
	_row.push_back(0);
	for (std::vector<Monolayer>::const_iterator layer = s.layer_vec.begin(); layer < s.layer_vec.end(); layer++)
	{
		std::cout << "Calculating XRF for layer " << (*layer).layer << std::endl;
		Xrf *temp_p = new Xrf(il.ev0, *layer), & temp = *temp_p;
		std::cout << "Calculating XRF attenuation for layer " << (*layer).layer << std::endl;
		for (int line=0; line < temp.lines.size(); line++)
		{
			for (std::vector<double>::const_iterator theta = omega.theta.begin(); theta < omega.theta.end(); theta++)
			{
				for (std::vector<double>::const_iterator beta = omega.beta.begin(); beta < omega.beta.end(); beta++)
				{
					double psiprime = il.psi_prime(*theta, *beta);
					//attenuation within the layer
					temp._y_mat.push_back(temp.y_vec[line]*atten_mono(il.ev0, temp.ev_vec[line], il.psi, psiprime, *layer)*omega.domega(*theta,*beta)/(4*Pi));
					if (psiprime > 0) //reflection geometry
						//attenuation from layers upstream
						for (std::vector<Monolayer>::const_iterator layer_up = s.layer_vec.begin(); layer_up < layer; layer_up++)
							temp._y_mat.back() *= atten_refl(il.ev0, temp.ev_vec[line], il.psi, psiprime, *layer_up);
					else if (psiprime < 0) //transmission geometry
					{
						//attenuation from layers upstream
						for (std::vector<Monolayer>::const_iterator layer_up = s.layer_vec.begin(); layer_up < layer; layer_up++)
							temp._y_mat.back() *= atten_trans_in(il.ev0, il.psi, *layer_up);
						//attenuation from layers downstream
						for (std::vector<Monolayer>::const_iterator layer_down = layer+1; layer_down < s.layer_vec.end(); layer_down++)
							temp._y_mat.back() *= atten_trans_out(temp.ev_vec[line], psiprime, *layer_down);
					}
				}
			}
		}
		add(temp, true);
		delete temp_p;
	}
	sum();
}

Xrf & Xrf::operator=(const Xrf & x)
{
	if (this != & x)
	{
		_y_mat = x.y_mat;
		_omega = x.omega;
		_ev_vec = x.ev_vec;
		_y_vec = x.y_vec;
		_Z_vec = x.Z_vec;
		_lines = x.lines;
		_row = x.row;
	}
	return *this;
}

void Xrf::add(const Xrf & x, bool mat_only)
{
	int Z, n_pixels;
	std::vector<int>::const_iterator Z_iter;
	if (!x.y_mat.empty())
		n_pixels = x.omega.theta.size()*x.omega.beta.size();
	for (int i = 0; i < x.Z_vec.size(); i++)
	{
		Z = x.Z_vec[i];
		Z_iter = std::find(Z_vec.begin(), Z_vec.end(), Z);
		if (Z_iter != Z_vec.end())
		{
			int p = Z_iter - Z_vec.begin();
			std::cout << "Z = " << Z << ", p = " << p << std::endl;
			if (!mat_only)
			{
				std::vector<double>::const_iterator i2 = x.y_vec.begin()+x.row[i];
				for (std::vector<double>::iterator i1 = _y_vec.begin()+row[p];
					i1 < _y_vec.begin()+row[p+1]; i1++, i2++)
					*i1 += *i2;
			}
			if (n_pixels >= 1)
			{
				std::vector<double>::const_iterator i2 = x.y_mat.begin()+x.row[i]*n_pixels;
				for (std::vector<double>::iterator i1 = _y_mat.begin()+row[p]*n_pixels;
					i1 < _y_mat.begin()+row[p+1]*n_pixels; i1++, i2++)
					*i1 += *i2;
			}
		}
		else
		{
			_Z_vec.push_back(Z);
			_lines.insert(_lines.end(),x.lines.begin()+x.row[i],x.lines.begin()+x.row[i+1]);
			_row.push_back(lines.size());
			_ev_vec.insert(_ev_vec.end(),x.ev_vec.begin()+x.row[i],x.ev_vec.begin()+x.row[i+1]);
			_y_vec.insert(_y_vec.end(),x.y_vec.begin()+x.row[i],x.y_vec.begin()+x.row[i+1]);
			if (n_pixels)
				_y_mat.insert(_y_mat.end(),x.y_mat.begin()+x.row[i]*n_pixels,x.y_mat.begin()+x.row[i+1]*n_pixels);
		}
	}
}

void Xrf::sum()
{
	if (!y_mat.empty())
	{
		_y_vec.clear();
		int n_pixels = omega.theta.size()*omega.beta.size();
		for (int i = 0; i < ev_vec.size(); i++)
			_y_vec.push_back(std::accumulate(y_mat.begin()+i*n_pixels,y_mat.begin()+(i+1)*n_pixels,0.0));
	}
}

void Xrf::show() const
{
	std::vector<int>::const_iterator Z, r, line;
	std::vector<double>::const_iterator ev, y;
	r = row.begin();
	y = y_vec.begin();
	ev = ev_vec.begin();
	line = lines.begin();
	
	for (Z = Z_vec.begin(); Z < Z_vec.end(); Z++)
	{
		std::cout << "******* XRF lines for " << AtomicNumberToSymbol(*Z) << " (" << *Z << ") : *******" << std::endl;
		for (int i=(*r); i<(*(r+1)); i++)
		{
			std::cout << "line " << *(line++) << ", ";
			std::cout << *(ev++) << " eV, ";
			std::cout << "intensity " << *(y++) << std::endl;
		}
		r++;
	}
	std::cout << "******* End *******" << std::endl;
}

Rayleigh::Rayleigh()
	: y_mat(_y_mat), omega(_omega), y(_y), ev0(_ev0)
{}

Rayleigh::Rayleigh(double ev0_, const Compound & c, const solid_angle & omega_, bool calc_sum)
	: y_mat(_y_mat), omega(_omega), y(_y), ev0(_ev0)
{
	_ev0 = ev0_;
	_omega = omega_;
	double theta, beta;
	for (int j = 0; j < omega.theta.size(); j++)
	{
		theta = omega.theta[j];
		for (int k = 0; k < omega.beta.size(); k++)
		{
			beta = omega.beta[k];
			_y_mat.push_back(c.dmac_rayleigh_pol(ev0, theta, beta)*omega.domega(theta, beta));
		}
	}
	if (calc_sum)
		sum();
}

Rayleigh::Rayleigh(const Sample & s, const Illumination & il, const solid_angle & omega_)
	: y_mat(_y_mat), omega(_omega), y(_y), ev0(_ev0)
{
	_ev0 = il.ev0;
	_omega = omega_;
	for (std::vector<Monolayer>::const_iterator layer = s.layer_vec.begin(); layer < s.layer_vec.end(); layer++)
	{
		std::cout << "Calculating Rayleigh scattering for layer " << (*layer).layer << std::endl;
		Rayleigh *temp_p = new Rayleigh(il.ev0, *layer, omega, false), & temp = *temp_p;
		int i=0;
		std::cout << "Calculating Rayleigh scattering attenuation for layer " << (*layer).layer << std::endl;
		for (std::vector<double>::const_iterator theta = omega.theta.begin(); theta < omega.theta.end(); theta++)
		{
			for (std::vector<double>::const_iterator beta = omega.beta.begin(); beta < omega.beta.end(); beta++)
			{
				double psiprime = il.psi_prime(*theta, *beta);
				//attenuation within the layer
				temp._y_mat[i] *= atten_mono(il.ev0, il.ev0, il.psi, psiprime, *layer);
				if (psiprime > 0) //reflection geometry
					//attenuation from layers upstream
					for (std::vector<Monolayer>::const_iterator layer_up = s.layer_vec.begin(); layer_up < layer; layer_up++)
						temp._y_mat[i] *= atten_refl(il.ev0, il.ev0, il.psi, psiprime, *layer_up);
				else if (psiprime < 0) //transmission geometry
				{
					//attenuation from layers upstream
					for (std::vector<Monolayer>::const_iterator layer_up = s.layer_vec.begin(); layer_up < layer; layer_up++)
						temp._y_mat[i]*= atten_trans_in(il.ev0, il.psi, *layer_up);
					//attenuation from layers downstream
					for (std::vector<Monolayer>::const_iterator layer_down = layer+1; layer_down < s.layer_vec.end(); layer_down++)
						temp._y_mat[i] *= atten_trans_out(il.ev0, psiprime, *layer_down);
				}
				i++;
			}
		}
		add(temp, true);
		delete temp_p;
	}
	sum();
}

Rayleigh & Rayleigh::operator=(const Rayleigh & r)
{
	if (this != & r)
	{
		_y_mat = r.y_mat;
		_omega = r.omega;
		_y = r.y;
		_ev0 = r.ev0;
	}
	return *this;
}

void Rayleigh::add(const Rayleigh & x, bool mat_only)
{
	if (!x.y_mat.empty())
	{
		if (y_mat.empty())
			_y_mat = x.y_mat;
		else
			for (int i = 0; i < y_mat.size(); i++)
				_y_mat[i] = x.y_mat[i];
	}
	if (!mat_only)
		_y += x.y;
}

void Rayleigh::sum()
{
	if (!y_mat.empty())
		_y = std::accumulate(y_mat.begin(), y_mat.end(), 0.0);
}

void Rayleigh::show() const
{
	std::cout << "******* Rayleigh scattering *******" << std::endl;
	std::cout << "Peak @ " << ev0 << std::endl;
	std::cout << "Intensity: " << y << std::endl;
	std::cout << "******* End *******" << std::endl;
}


Compton::Compton()
	: y_mat(_y_mat), omega(_omega), y_vec(_y_vec), ev_vec(_ev_vec), ev0(_ev0)
{}

Compton::Compton(double ev0_, const Compound & c, const solid_angle & omega_, bool calc_ev, bool calc_sum)
	: y_mat(_y_mat), omega(_omega), y_vec(_y_vec), ev_vec(_ev_vec), ev0(_ev0)
{
	_ev0 = ev0_;
	_omega = omega_;
	double theta, beta;
	for (int j = 0; j < omega.theta.size(); j++)
	{
		theta = omega.theta[j];
		if (theta < 1e-6)
		{
			_y_mat.resize(y_mat.size()+omega.theta.size(),0);
			if (calc_ev)
				_ev_vec.push_back(ev0);
			if (calc_sum)
				_y_vec.push_back(0);
			continue;
		}
		for (int k = 0; k < omega.beta.size(); k++)
		{
			beta = omega.beta[k];
			_y_mat.push_back(c.dmac_compton_pol(ev0, theta, beta)*omega.domega(theta, beta));
		}
		if (calc_ev)
			_ev_vec.push_back(ev_scattered(theta, ev0));
		if (calc_sum)
			_y_vec.push_back(std::accumulate(y_mat.begin()+j*omega.beta.size(),y_mat.begin()+(j+1)*omega.beta.size(),0.0));
	}
}

Compton::Compton(const Sample & s, const Illumination & il, const solid_angle & omega_)
	: y_mat(_y_mat), omega(_omega), y_vec(_y_vec), ev_vec(_ev_vec), ev0(_ev0)
{
	_ev0 = il.ev0;
	_omega = omega_;
	for (std::vector<Monolayer>::const_iterator layer = s.layer_vec.begin(); layer < s.layer_vec.end(); layer++)
	{
		std::cout << "Calculating Compton scattering for layer " << (*layer).layer << std::endl;
		Compton *temp_p = new Compton(il.ev0, *layer, omega, false, false), & temp = *temp_p;
		std::cout << "Calculating Compton scattering attenuation for layer " << (*layer).layer << std::endl;
		int i = 0;
		for (std::vector<double>::const_iterator theta = omega.theta.begin(); theta < omega.theta.end(); theta++)
		{
			for (std::vector<double>::const_iterator beta = omega.beta.begin(); beta < omega.beta.end(); beta++)
			{
				double psiprime = il.psi_prime(*theta, *beta);
				//attenuation within the layer
				temp._y_mat[i] *= atten_mono(il.ev0, il.ev0, il.psi, psiprime, *layer);
				if (psiprime > 0) //reflection geometry
					//attenuation from layers upstream
					for (std::vector<Monolayer>::const_iterator layer_up = s.layer_vec.begin(); layer_up < layer; layer_up++)
						temp._y_mat[i] *= atten_refl(il.ev0, il.ev0, il.psi, psiprime, *layer_up);
				else if (psiprime < 0) //transmission geometry
				{
					//attenuation from layers upstream
					for (std::vector<Monolayer>::const_iterator layer_up = s.layer_vec.begin(); layer_up < layer; layer_up++)
						temp._y_mat[i] *= atten_trans_in(il.ev0, il.psi, *layer_up);
					//attenuation from layers downstream
					for (std::vector<Monolayer>::const_iterator layer_down = layer+1; layer_down < s.layer_vec.end(); layer_down++)
						temp._y_mat[i] *= atten_trans_out(il.ev0, psiprime, *layer_down);
				}
				i++;
			}
		}
		add(temp, true);
		// temp.sum_ev();
		// for (std::vector<double>::const_iterator xx = y_mat.begin(); xx < y_mat.end(); xx++)
		// {
			// if (*xx < 1e-50)
			// std::cout << *xx << " ";
		// }
		// std::cout << std::endl;
		// int ttt;
		// std::cin >> ttt;
		// temp.show();
		delete temp_p;
	}
	sum_ev();
}


Compton & Compton::operator=(const Compton & c)
{
	if (this != & c)
	{
		_y_mat = c.y_mat;
		_omega = c.omega;
		_y_vec = c.y_vec;
		_ev_vec = c.ev_vec;
		_ev0 = ev0;
	}
	return *this;
}

void Compton::add(const Compton & x, bool mat_only)
{
	if (!x.y_mat.empty())
	{
		if (y_mat.empty())
			_y_mat = x.y_mat;
		else
			for (int i = 0; i < y_mat.size(); i++)
				_y_mat[i] += x.y_mat[i];
	}
	if (!mat_only)
		for (int j = 0; j < x.y_vec.size(); j++)
			_y_vec[j] += x.y_vec[j];
}

void Compton::sum()
{
	if (!y_mat.empty())
	{
		_y_vec.clear();
		for (int j = 0; j < omega.theta.size(); j++)
			_y_vec.push_back(std::accumulate(y_mat.begin()+j*omega.beta.size(),y_mat.begin()+(j+1)*omega.beta.size(),0.0));
	}
}

void Compton::sum_ev()
{
	if (!y_mat.empty())
	{
		_y_vec.clear();
		_ev_vec.clear();
		double theta;
		for (int j = 0; j < omega.theta.size(); j++)
		{
			theta = omega.theta[j];
			_y_vec.push_back(std::accumulate(y_mat.begin()+j*omega.beta.size(),y_mat.begin()+(j+1)*omega.beta.size(),0.0));
			_ev_vec.push_back(ev_scattered(theta, ev0));
		}
	}
}

void Compton::ev()
{
	_ev_vec.clear();
	double theta;
	for (int j = 0; j < omega.theta.size(); j++)
	{
		theta = omega.theta[j];
		_ev_vec.push_back(ev_scattered(theta, ev0));
	}
}

void Compton::show() const
{
	std::cout << "******* Compton scattering *******" << std::endl;
	std::cout << "Energy range " << ev_vec.back() << " - " << ev_vec.front() << std::endl;
	std::cout << "Total intensity: " << std::accumulate(y_vec.begin(), y_vec.end(), 0.0) << std::endl;
	std::cout << "******* End *******" << std::endl;
}

Spectrum::Spectrum() 
	: y_vec(_y_vec), xrf(_xrf), ray(_ray), comp(_comp), omega(_omega), sample(_sample), illumination(_illumination), detector(_detector), Z_vec(xrf.Z_vec), lines(xrf.lines), row(xrf.row) {}

Spectrum::Spectrum(const Sample & s, const Illumination & il, const solid_angle & omega_, const Detector & det) 
	: y_vec(_y_vec), xrf(_xrf), ray(_ray), comp(_comp), omega(_omega), sample(_sample), illumination(_illumination), detector(_detector), Z_vec(xrf.Z_vec), lines(xrf.lines), row(xrf.row)
{
	_omega = omega_;
	_illumination = il;
	_sample = s;
	_detector = det;
	
	_xrf = Xrf(sample, illumination, omega);
	_ray = Rayleigh(sample, illumination, omega);
	_comp = Compton(sample, illumination, omega);
	
	_y_vec.clear();
	_y_vec.resize(detector.channel.n_channels);
	detector.genspec(xrf.ev_vec, xrf.y_vec, _y_vec, false, false);
	detector.genspec(comp.ev_vec, comp.y_vec, _y_vec, false, false);
	detector.genspec(illumination.ev0, ray.y, _y_vec, false, false);
}

Spectrum::Spectrum(double ev0, const Compound & c, const solid_angle & omega_, const Detector & det) 
	: y_vec(_y_vec), xrf(_xrf), ray(_ray), comp(_comp), omega(_omega), sample(_sample), illumination(_illumination), detector(_detector), Z_vec(xrf.Z_vec), lines(xrf.lines), row(xrf.row)
{
	_omega = omega_;
	_detector = det;
	
	_xrf = Xrf(ev0, c, omega);
	_ray = Rayleigh(ev0, c, omega);
	_comp = Compton(ev0, c, omega);
	
	_y_vec.clear();
	_y_vec.resize(detector.channel.n_channels);
	std::cout << "Generating detector binned spectra..." << std::endl;
	detector.genspec(xrf.ev_vec, xrf.y_vec, _y_vec);
	detector.genspec(comp.ev_vec, comp.y_vec, _y_vec);
	detector.genspec(illumination.ev0, ray.y, _y_vec);
}

void Spectrum::show() const
{
	xrf.show();
	ray.show();
	comp.show();
}
