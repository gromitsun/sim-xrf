/*input.hpp*/
#ifndef INPUT_HPP
#define INPUT_HPP

#include "detector.hpp"
#include "geometry.hpp"
#include "sample.hpp"

void readfile(std::string filename, Sample & sp, Illumination & il, solid_angle & omega, Detector & det);

#endif