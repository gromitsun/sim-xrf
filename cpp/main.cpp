#include <iostream>
#include <fstream>

#include "spectrum.hpp"
#include "input.hpp"

int main()
{
	Sample sp;
	Illumination il;
	solid_angle omega;
	Detector det;
	readfile("input.txt", sp, il, omega, det);
	sp.show();
	omega.show();
	Spectrum spec(sp, il, omega, det);
	spec.show();
	return 1;
}