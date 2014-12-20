/*input.cpp*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "input.hpp"
#include "constants.hpp"
#include "math.hpp"

inline std::string trim_comment(const std::string & s, const std::string & delimiter="#")
{
    if (s.empty())
        return s;
    else
        return s.substr(0, s.find(delimiter));
}

inline std::string trim_right(
  const std::string & s,
  const std::string & delimiters = " \f\n\r\t\v" )
{
	if (s.empty())
		return s;
	return s.substr(0, s.find_last_not_of(delimiters) + 1);
}

inline std::string trim_left(
  const std::string & s,
  const std::string & delimiters = " \f\n\r\t\v" )
{
	if (s.empty())
		return s;
	return s.substr(s.find_first_not_of(delimiters));
}

inline std::string trim(
  const std::string & s,
  const std::string & delimiters = " \f\n\r\t\v" )
{
	return trim_left(trim_right(trim_comment(s), delimiters), delimiters);
}

inline std::vector<int> parse_vec_int(const std::string & s)
{
	std::vector<int> v;
	std::istringstream is(s);
	std::string t;
	while(std::getline(is, t, ','))
	// {
		// std::cout << trim(t) << std::endl;
		v.push_back(std::stoi(trim(t)));
	// }
	return v;
}

inline std::vector<double> parse_vec_double(const std::string & s)
{
	std::vector<double> v;
	std::istringstream is(s);
	std::string t;
	while(std::getline(is, t, ','))
		v.push_back(std::stod(trim(t)));
	return v;
}

inline bool parse_parameter(const std::string & line, std::string & name, std::string & value, const std::string & sep = "=")
{
	std::size_t pos = line.find(sep);
	if (pos == std::string::npos)
		return false;
	else
	{
		name = trim(line.substr(0, pos-1));
		value = trim(line.substr(pos+1));
		return true;
	}
}

inline double deg_to_rad(const double & deg)
{
	return deg*Pi/180;
}


void readfile(std::string filename, Sample & sp, Illumination & il, solid_angle & omega, Detector & det)
{
	std::ifstream fin;
	std::string s;
	fin.open(filename.c_str());

	std::size_t pos;
	while (getline(fin, s))
	{
		if (trim(s).empty())
			continue;
		else if (s.find("Sample:") != std::string::npos)
		{
			std::cout << "Reading in sample data..." << std::endl;
			// read sample
			while (getline(fin, s) && !trim(s).empty())
			{
				int layer = std::stoi(trim(s.substr(s.find("layer")+5, s.find(":")-6)));
				std::cout << "Reading layer #" << layer << "..." << std::endl;
				//read layer
				std::vector<int> Z_vec;
				std::vector<double> p_vec;
				double thickness, density;
				while (getline(fin, s) && !trim(s).empty())
				{
					std::string name, value;
					if (parse_parameter(s, name, value))
					{
						if (name == "Z")
						{
							Z_vec = parse_vec_int(value);					
						}
						else if (name == "p")
						{
							p_vec = parse_vec_double(value);
						}
						else if (name == "thickness")
						{
							thickness = std::stod(value);
						}
						else if (name == "density")
						{
							density = std::stod(value);
						}
					}
				}
				sp.add_layer(Monolayer(Z_vec, p_vec, density, thickness, layer));
			}
			
		}
		else if (s.find("Detector:") != std::string::npos)
		{
			std::cout << "Reading in detector data..." << std::endl;
			// read detector
			double ev_offset, ev_gain, n_channels, 
				noise, fano, gamma, fs, ft, 
				thickness, density;
			std::string material;
			while (getline(fin, s) && !trim(s).empty())
			{
				std::string name, value;
				if(parse_parameter(s, name, value))
				{
					if (name == "ev_offset")
					{
						ev_offset = std::stod(value);
					}
					else if (name == "ev_gain")
					{
						ev_gain = std::stod(value);
					}
					else if (name == "n_channels")
					{
						n_channels = std::stod(value);
					}
					else if (name == "noise")
					{
						noise = std::stod(value);
					}
					else if (name == "fano")
					{
						fano = std::stod(value);
					}
					else if (name == "gamma")
					{
						gamma = std::stod(value);
					}
					else if (name == "fs")
					{
						fs = std::stod(value);
					}
					else if (name == "ft")
					{
						ft = std::stod(value);
					}
					else if (name == "material")
					{
						material = value;
					}
					else if (name == "thickness")
					{
						thickness = std::stod(value);
					}
					else if (name == "density")
					{
						density = std::stod(value);
					}
				}
				
				det = Detector(Channel(ev_offset, ev_gain, n_channels), 
					Response(noise, fano, gamma, fs, ft, ev_gain), 
					Window(material, thickness, density));
			}
		}
		else if (s.find("Illumination:") != std::string::npos)
		{
			std::cout << "Reading in illumination data..." << std::endl;
			// read illumination
			double ev0, psi, alpha, n_photons=-1, beam_cross_section=-1;
			while (getline(fin, s) && !trim(s).empty())
			{
				std::string name, value;
				if(parse_parameter(s, name, value))
				{
					if (name == "ev0")
					{
						ev0 = std::stod(value);
					}
					else if (name == "psi")
					{
						psi = deg_to_rad(std::stod(value));
					}
					else if (name == "alpha")
					{
						alpha = deg_to_rad(std::stod(value));
					}
					else if (name == "n_photons")
					{
						n_photons = std::stod(value);
					}
					else if (name == "beam_cross_section")
					{
						beam_cross_section = std::stod(value);
					}
					else if (name == "beam_diameter")
					{
						beam_cross_section = sq(std::stod(value)/2.)*Pi;
					}
				}
				
			}
			il = Illumination(ev0, psi, alpha, n_photons, beam_cross_section);
		}
		
		else if (s.find("Solid angle:") != std::string::npos)
		{
			std::cout << "Reading in solid angle data..." << std::endl;
			// read solid angle
			std::vector<double> ar;
			double theta_inc, beta_inc;
			while (getline(fin, s) && !trim(s).empty())
			{
				std::string name, value;
				if(parse_parameter(s, name, value))
				{
					if (name == "angle_range")
					{
						ar = parse_vec_double(value);
						for (std::vector<double>::iterator i = ar.begin(); i < ar.end(); i++)
							*i = deg_to_rad(*i);
					}
					else if (name == "theta_inc")
					{
						theta_inc = deg_to_rad(std::stod(value));
					}
					else if (name == "beta_inc")
					{
						beta_inc = deg_to_rad(std::stod(value));
					}
				}
			}
			omega = solid_angle(ar.data(), theta_inc, beta_inc);
		}
	}
	std::cout << "Done!" << std::endl;
}


