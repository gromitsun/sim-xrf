/*detector.hpp*/
#ifndef DETECTOR_HPP
#define DETECTOR_HPP

#include<vector>
#include<string>


class Channel
{
private:
	double _ev_offset;
	double _ev_gain;
	int _n_channels;
public:
	Channel(double ev_offset_ = 0, double ev_gain_ = 10, int n_channels_ = 2048);
	// ~Channel();
	Channel & operator=(const Channel & ch);
	const double & ev_offset;
	const double & ev_gain;
	const int & n_channels;
	int ev_to_channel(double ev) const;
	void bin(const std::vector<double> & ev_raw, const std::vector<double> & y_raw, std::vector<double> & y_binned) const;
	void bin(const double & ev_raw, const double & y_raw, std::vector<double> & y_binned) const;
};

class Response
{
private:
	double _noise;
	double _fano;
	double _gamma;
	double _fs;
	double _ft;
	double _ev_gain;
public:
	Response(double noise_ = 100,
		double fano_ = 0.114,
		double gamma_ = 2.5,
		double fs_ = 0.03,
		double ft_ = 0.02,
		double ev_gain_ = 0);
	// ~Response();
	Response & operator=(const Response & r);
	const double & noise;
	const double & fano;
	const double & gamma;
	const double & fs;
	const double & ft;
	const double & ev_gain;
	
	void set_gain(const double ev_gain);
	double FWHM(double ev);
};

		
class Window
{
private:
	std::string _material;
	double _thickness;
	double _density;
public:
	Window(std::string material_ = "Be",
	double thickness_ = 24e-4,
	double density_ = 1);
	// ~Window();
	Window & operator=(const Window & w);
	const std::string & material;
	const double & thickness;
	const double & density;
	double transmission(double ev) const;
};
		
class Detector
{
private:
	Channel _channel;
	Response _response;
	Window _window;
public:
	Detector(Channel channel_ = Channel(), Response response_ = Response(), Window window_ = Window());
	// ~Detector();
	Detector & operator=(const Detector & d);
	const Channel & channel;
	const Response & response;
	const Window & window;
	void genspec(const std::vector<double> & ev_raw, const std::vector<double> & y_raw, std::vector<double> & y_binned, bool det_response = true, bool det_window = true) const;
	void genspec(const double & ev_raw, const double & y_raw, std::vector<double> & y_binned, bool det_response = true, bool det_window = true) const;
};

#endif