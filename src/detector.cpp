/*detector.cpp*/
#include <iostream>
#include<cmath>

#include "xraylib.h"
#include "detector.hpp"
#include "math.hpp"
#include "constants.hpp"

#define MP 1

Channel::Channel(double ev_offset_, double ev_gain_, int n_channels_) 
	: ev_offset(_ev_offset), ev_gain(_ev_gain), n_channels(_n_channels)
{
	_ev_offset = ev_offset_;
	_ev_gain = ev_gain_;
	_n_channels = n_channels_;
}

Channel & Channel::operator=(const Channel & ch)
{
	if (this != & ch)
	{
		_ev_offset = ch.ev_offset;
		_ev_gain = ch.ev_gain;
		_n_channels = ch.n_channels;
	}
	return *this;
}

int Channel::ev_to_channel(double ev) const
{
	return int((ev - ev_offset)/ev_gain + 0.5);
}

void Channel::bin(const std::vector<double> & ev_raw, const std::vector<double> & y_raw, std::vector<double> & y_binned) const
{
	if (y_binned.size() != n_channels)
	{
		std::cout << "Error: Bin size does not match the number of channels!" << std::endl;
		return;
	}
	#pragma omp parallel for if(MP)
	for (int i=0; i<ev_raw.size(); i++)
		y_binned[ev_to_channel(ev_raw[i])] += y_raw[i];
}

void Channel::bin(const double & ev_raw, const double & y_raw, std::vector<double> & y_binned) const
{
	if (y_binned.size() != n_channels)
	{
		std::cout << "Error: Bin size does not match the number of channels!" << std::endl;
		return;
	}
	y_binned[ev_to_channel(ev_raw)] += y_raw;
}

void Channel::bin(const std::vector<double> & ev_raw, const std::vector<double> & y_raw, std::vector<double> & y_binned, std::vector<double> & y_separate, const std::vector<int> & row) const
{
	if (y_binned.size() != n_channels)
	{
		std::cout << "Error: Bin size does not match the number of channels!" << std::endl;
		return;
	}
	int n_old = y_separate.size();
	y_separate.resize(n_old+n_channels);
	#pragma omp parallel for collapse(2) if(MP)
	for (int i = 0; i < row.size()-1; i++)
		for (int j = row[i]; j < ((row[i+1] > 0) ? row[i+1] : ev_raw.size()); j++)
		{
			int chn = ev_to_channel(ev_raw[j]);
			y_binned[chn] += y_raw[j];
			// for separate output
			y_separate[n_old + i*n_channels + chn] = y_raw[j];
		}
}

void Channel::bin(const double & ev_raw, const double & y_raw, std::vector<double> & y_binned, std::vector<double> & y_separate) const
{
	if (y_binned.size() != n_channels)
	{
		std::cout << "Error: Bin size does not match the number of channels!" << std::endl;
		return;
	}
	int n_old = y_separate.size();
	y_separate.resize(n_old+n_channels);
	int chn = ev_to_channel(ev_raw);
	y_binned[chn] += y_raw;
	y_separate[n_old+chn] = y_raw;
}

void Channel::show() const
{
    out(std::cout);
}

void Channel::out(std::ostream & ost) const
{
	ost << "# Detector channel parameters:" << std::endl;
	ost << "# \tev_offset = " << ev_offset << std::endl;
	ost << "# \tev_gain = " << ev_gain << std::endl;
	ost << "# \tn_channels = " << n_channels << std::endl;
    ost << std::endl;
}


Response::Response(double noise_,
		double fano_,
		double gamma_,
		double fs_,
		double ft_,
		double ev_gain_)
		: noise(_noise), fano(_fano), gamma(_gamma), fs(_fs), ft(_ft), ev_gain(_ev_gain)
{
	_noise = noise_;
	_fano = fano_;
	_gamma = gamma_;
	_fs = fs_;
	_ft = ft_;
	_ev_gain = ev_gain_;
}

Response & Response::operator=(const Response & r)
{
	if (this != & r)
	{
		_fano = r.fano;
		_gamma = r.gamma;
		_fs = r.fs;
		_ft = r.ft;
		_ev_gain = r.ev_gain;
	}
	return *this;
}

void Response::set_gain(const double ev_gain)
{
	_ev_gain = ev_gain;
}

double Response::FWHM(double ev)
{
	return 2.3548*std::sqrt(sq(noise/2.3548)+3.58*fano*ev);
}

void Response::show() const
{
    out(std::cout);
}

void Response::out(std::ostream & ost) const
{
	ost << "# Detector response parameters:" << std::endl;
	ost << "# \tfano = " << fano << std::endl;
	ost << "# \tgamma = " << gamma << std::endl;
	ost << "# \tfs = " << fs << std::endl;
	ost << "# \tft = " << ft << std::endl;
    ost << std::endl;
}



Window::Window(std::string material_ ,
	double thickness_,
	double density_)
	: material(_material), thickness(_thickness), density(_density)
{
	_material = material_;
	_thickness = thickness_;
	_density = density_;
}

Window & Window::operator=(const Window & w)
{
	if (this != & w)
	{
		_material = w.material;
		_thickness = w.thickness;
		_density = w.density;
	}	
	return *this;
}

double Window::transmission(double ev) const
{
	double _mac = CS_Total_CP(material.c_str(), ev/1000.);
	return std::exp(-_mac*density*thickness);
}

void Window::show() const
{
    out(std::cout);
}

void Window::out(std::ostream & ost) const
{
	ost << "# Detector window parameters:" << std::endl;
	ost << "# \tmaterial = " << material << std::endl;
	ost << "# \tthickness = " << thickness << std::endl;
	ost << "# \tdensity = " << density << std::endl;
    ost << std::endl;
}


Detector::Detector(Channel channel_,
			Response response_,
			Window window_)
			: channel(_channel), response(_response), window(_window)
{
	_channel = channel_;
	_response = response_;
	_window = window_;
	_response.set_gain(channel.ev_gain);
}

Detector & Detector::operator=(const Detector & d)
{
	if (this != & d)
	{
		_channel = d.channel;
		_response = d.response;
		_window = d.window;
	}
	return *this;
}

void Detector::genspec(const std::vector<double> & ev_raw, const std::vector<double> & y_raw, std::vector<double> & y_binned, bool det_response, bool det_window) const
{
	if (y_binned.size() != channel.n_channels)
	{
		std::cout << "Error: Bin size does not match the number of channels!" << std::endl;
		return;
	}
	//apply detector response
	if (det_response)
	{
		std::cout << "Generating spectrum with detector response..." << std::endl;
		double sigma, y_temp;
		for (std::vector<double>::const_iterator ev = ev_raw.begin(), y = y_raw.begin();
			ev < ev_raw.end(); ev++, y++)
		{
			sigma = std::sqrt(sq(response.noise/2.3548)+3.58*response.fano*(*ev));
			y_temp = *y;
			//apply detector filter window
			if (det_window)
				y_temp *= window.transmission(*ev);
			#pragma omp parallel for if(MP)
			for (int i = 0; i < y_binned.size(); i++)
			{
				double ev_ch = channel.ev_offset + i*channel.ev_gain;
				//Gaussian
				y_binned[i] += y_temp*channel.ev_gain/(sigma*std::sqrt(2*Pi))*std::exp(-sq((*ev)-ev_ch)/(2*sq(sigma)));
				//Step function
				y_binned[i] += response.fs*y_temp*channel.ev_gain/(2.*(*ev))*std::erfc((ev_ch-(*ev))/(std::sqrt(2)*sigma));
				//Tailing function
				y_binned[i] += response.ft*y_temp*channel.ev_gain/(2.*response.gamma*sigma*std::exp(-1./(2*sq(response.gamma))))*std::exp((ev_ch-(*ev))/(response.gamma*sigma))*std::erfc((ev_ch-(*ev))/(std::sqrt(2)*sigma)+1/(std::sqrt(2)*response.gamma));
			}
		}
	}
	else
	{
		std::cout << "Generating spectrum without detector response..." << std::endl;
		if (det_window)//apply detector filter window
			#pragma omp parallel for if(MP)
			for (int i = 0;	i < ev_raw.size(); i++)
				y_binned[channel.ev_to_channel(ev_raw[i])] += (y_raw[i])*window.transmission(ev_raw[i]);
		else
			channel.bin(ev_raw, y_raw, y_binned);
	}
}

void Detector::genspec(const double & ev_raw, const double & y_raw, std::vector<double> & y_binned, bool det_response, bool det_window) const
{
	if (y_binned.size() != channel.n_channels)
	{
		std::cout << "Error: Bin size does not match the number of channels!" << std::endl;
		return;
	}
	//apply detector response
	if (det_response)
	{
		std::cout << "Generating spectrum with detector response..." << std::endl;
		double sigma, y_temp;
		sigma = std::sqrt(sq(response.noise/2.3548)+3.58*response.fano*ev_raw);
		y_temp = y_raw;
		//apply detector filter window
		if (det_window)
			y_temp *= window.transmission(ev_raw);
		for (int i = 0; i < y_binned.size(); i++)
		{
			double ev_ch = channel.ev_offset + i*channel.ev_gain;
			//Gaussian
			y_binned[i] += y_temp*channel.ev_gain/(sigma*std::sqrt(2*Pi))*std::exp(-sq(ev_raw-ev_ch)/(2*sq(sigma)));
			//Step function
			y_binned[i] += response.fs*y_temp*channel.ev_gain/(2.*ev_raw)*std::erfc((ev_ch-ev_raw)/(std::sqrt(2)*sigma));
			//Tailing function
			y_binned[i] += response.ft*y_temp*channel.ev_gain/(2.*response.gamma*sigma*std::exp(-1./(2*sq(response.gamma))))*std::exp((ev_ch-ev_raw)/(response.gamma*sigma))*std::erfc((ev_ch-ev_raw)/(std::sqrt(2)*sigma)+1/(std::sqrt(2)*response.gamma));
		}
	}
	else
	{
		std::cout << "Generating spectrum without detector response..." << std::endl;
		if (det_window)//apply detector filter window
			y_binned[channel.ev_to_channel(ev_raw)] += y_raw*window.transmission(ev_raw);
		else
			channel.bin(ev_raw, y_raw, y_binned);
	}
}

void Detector::genspec(const std::vector<double> & ev_raw, const std::vector<double> & y_raw, std::vector<double> & y_binned, std::vector<double> & y_separate, const std::vector<int> & row, bool det_response, bool det_window) const
{
	if (y_binned.size() != channel.n_channels)
	{
		std::cout << "Error: Bin size does not match the number of channels!" << std::endl;
		return;
	}
	//apply detector response
	if (det_response)
	{
		std::cout << "Generating spectrum with detector response..." << std::endl;

		std::vector<double>::iterator ys0; // ys0 pointing to the first channel
		int n_old = y_separate.size();
		y_separate.resize(n_old+(row.size()-1)*channel.n_channels);
		ys0 = y_separate.begin() + n_old;

		int i = 0;
		double sigma, y_temp;
		for (std::vector<double>::const_iterator ev = ev_raw.begin(), y = y_raw.begin();
			ev < ev_raw.end(); ev++, y++)
		{
			// for separate output
			if (ev - ev_raw.begin() == row[i+1])
				i++;
			
			sigma = std::sqrt(sq(response.noise/2.3548)+3.58*response.fano*(*ev));
			y_temp = *y;
			//apply detector filter window
			if (det_window)
				y_temp *= window.transmission(*ev);
			#pragma omp parallel for if(MP)
			for (int j = 0; j < y_binned.size(); j++)
			{
				double y_temp1 = 0;
				double ev_ch = channel.ev_offset + j*channel.ev_gain;
				
				//Gaussian
				y_temp1 += y_temp*channel.ev_gain/(sigma*std::sqrt(2*Pi))*std::exp(-sq((*ev)-ev_ch)/(2*sq(sigma)));
				//Step function
				y_temp1 += response.fs*y_temp*channel.ev_gain/(2.*(*ev))*std::erfc((ev_ch-(*ev))/(std::sqrt(2)*sigma));
				//Tailing function
				y_temp1 += response.ft*y_temp*channel.ev_gain/(2.*response.gamma*sigma*std::exp(-1./(2*sq(response.gamma))))*std::exp((ev_ch-(*ev))/(response.gamma*sigma))*std::erfc((ev_ch-(*ev))/(std::sqrt(2)*sigma)+1/(std::sqrt(2)*response.gamma));
				
				// Total output
				y_binned[j] += y_temp1;
				//Separate output
				*(ys0 + i*channel.n_channels + j) += y_temp1;
			}
		}
	}
	else
	{
		std::cout << "Generating spectrum without detector response..." << std::endl;
		std::vector<double>::iterator ys0; // ys0 pointing to the first channel
		int n_old = y_separate.size();
		y_separate.resize(n_old+(row.size()-1)*channel.n_channels);
		ys0 = y_separate.begin() + n_old;
		int i = 0;
		if (det_window)//apply detector filter window
		{
			#pragma omp parallel for collapse(2) if(MP)
			for (int i = 0; i < row.size()-1; i++)
				for (int j = row[i]; j < ((row[i+1] > 0) ? row[i+1] : ev_raw.size()); j++)
				{
					int chn = channel.ev_to_channel(ev_raw[j]);
					double y_temp = (y_raw[j])*window.transmission(ev_raw[j]);
					y_binned[chn] += y_temp;
					// for separate output
					*(ys0 + i*channel.n_channels + chn) = y_temp;
				}
		}
		else
			channel.bin(ev_raw, y_raw, y_binned, y_separate, row);
	}
}

void Detector::genspec(const double & ev_raw, const double & y_raw, std::vector<double> & y_binned, std::vector<double> & y_separate, bool det_response, bool det_window) const
{
	if (y_binned.size() != channel.n_channels)
	{
		std::cout << "Error: Bin size does not match the number of channels!" << std::endl;
		return;
	}
	//apply detector response
	if (det_response)
	{
		std::cout << "Generating spectrum with detector response..." << std::endl;
		std::vector<double>::iterator ys0; // ys0 pointing to the first channel
		int n_old = y_separate.size();
		y_separate.resize(n_old+channel.n_channels);
		ys0 = y_separate.begin() + n_old;
		double sigma, y_temp;
		sigma = std::sqrt(sq(response.noise/2.3548)+3.58*response.fano*ev_raw);
		y_temp = y_raw;
		//apply detector filter window
		if (det_window)
			y_temp *= window.transmission(ev_raw);
		#pragma omp parallel for if(MP)
		for (int i = 0; i < y_binned.size(); i++)
		{
			double y_temp1 = 0;
			double ev_ch = channel.ev_offset + i*channel.ev_gain;
			//Gaussian
			y_temp1 += y_temp*channel.ev_gain/(sigma*std::sqrt(2*Pi))*std::exp(-sq(ev_raw-ev_ch)/(2*sq(sigma)));
			//Step function
			y_temp1 += response.fs*y_temp*channel.ev_gain/(2.*ev_raw)*std::erfc((ev_ch-ev_raw)/(std::sqrt(2)*sigma));
			//Tailing function
			y_temp1 += response.ft*y_temp*channel.ev_gain/(2.*response.gamma*sigma*std::exp(-1./(2*sq(response.gamma))))*std::exp((ev_ch-ev_raw)/(response.gamma*sigma))*std::erfc((ev_ch-ev_raw)/(std::sqrt(2)*sigma)+1/(std::sqrt(2)*response.gamma));
			
			// Total output
			y_binned[i] += y_temp1;
			//Separate output
			*(ys0+i) = y_temp1;
		}
	}
	else
	{
		std::cout << "Generating spectrum without detector response..." << std::endl;
		std::vector<double>::iterator ys0; // ys0 pointing to the first channel
		int n_old = y_separate.size();
		y_separate.resize(n_old+channel.n_channels);
		ys0 = y_separate.begin() + n_old;
		
		if (det_window)//apply detector filter window
		{
			int chn = channel.ev_to_channel(ev_raw);
			double y_temp = y_raw*window.transmission(ev_raw);
			y_binned[chn] += y_temp;
			*(ys0+chn) = y_temp;
		}
		else
			channel.bin(ev_raw, y_raw, y_binned, y_separate);
	}
}

void Detector::show() const
{
    out(std::cout);
}

void Detector::out(std::ostream & ost) const
{
    ost << "# ======================================== #" << std::endl;
	ost << "# # # # # Detector parameters # # # # #" << std::endl;
	channel.out(ost);
	response.out(ost);
	window.out(ost);
	ost << "# # # # # End of detector parameters # # # # #" << std::endl;
    ost << "# ======================================== #" << std::endl;
}