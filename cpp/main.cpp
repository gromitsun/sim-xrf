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
	std::cout << "Reading file \"" << input_file << "\"..." << std::endl;
	readfile(input_file, sp, il, omega, det);
	
	// Show configurations
	sp.show();
	omega.show();
	det.show();

	// Calculate spectrum
	std::cout << "Calculating spectrum..." << std::endl;
	Spectrum spec(sp, il, omega, det, true, true, true);

	// Show results on screen
	spec.show();

	// Save results to file
	std::cout << "Writing output into \"" << output_file << "\"..." << std::endl;
	std::ofstream fout;
	fout.open(output_file);
	
	spec.out(fout);
	fout << std::endl;
	det.out(fout);
	std::cout << "Done!" << std::endl;
	
	return 1;
}