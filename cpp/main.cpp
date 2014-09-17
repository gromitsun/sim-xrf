#include <iostream>
#include <fstream>
#include <string>

#include "spectrum.hpp"
#include "input.hpp"

int main(int argc, char ** argv)
{
	std::string input_file = "input.txt", output_file = "out.txt";
	if (argc > 2)
		output_file = argv[2];
	else if (argc > 1)
		input_file = argv[1];
		
	Sample sp;
	Illumination il;
	solid_angle omega;
	Detector det;
	
	// Read input file
	readfile(input_file, sp, il, omega, det);
	
	// Show configurations
	sp.show();
	omega.show();
	
	// Calculate spectrum
	Spectrum spec(sp, il, omega, det, true, true, true);
	
	// Show results on screen
	spec.show();
	
	// Save results to file
	std::ofstream fout;
	fout.open(output_file);
	
	spec.out(fout);
	
	return 1;
}