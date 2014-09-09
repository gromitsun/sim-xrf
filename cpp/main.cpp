#include <iostream>
#include <fstream>

#include "spectrum.hpp"

int main()
{
	using namespace std;
	double ar[] = {0.5, 1, 0, 1.57};
	solid_angle omega(ar, 0.02, 0.02 );
	omega.show();
	
	
	
	Sample s;
	vector<int> *Z = new vector<int>;
	vector<double> *p = new vector<double>;
	(*Z).push_back(26);
	(*Z).push_back(6);
	(*Z).push_back(16);
	(*Z).push_back(12);
	(*Z).push_back(8);
	(*p).push_back(1);
	(*p).push_back(10);
	(*p).push_back(2);
	(*p).push_back(20);
	(*p).push_back(30);
	
	Compound c((*Z), (*p));
	
	
	s.add_layer(* (new Monolayer((*Z),(*p), 1.5, 2e-5, 1)));
	
	Z = new vector<int>;
	p = new vector<double>;
	(*Z).push_back(29);
	(*Z).push_back(6);
	(*Z).push_back(16);
	(*Z).push_back(12);
	(*Z).push_back(8);
	(*p).push_back(1);
	(*p).push_back(10);
	(*p).push_back(2);
	(*p).push_back(20);
	(*p).push_back(30);
	
	s.add_layer(* (new Monolayer((*Z),(*p), 1.2, 3e-5, 2)));
	
	Z = new vector<int>;
	p = new vector<double>;
	(*Z).push_back(30);
	(*Z).push_back(6);
	(*Z).push_back(16);
	(*Z).push_back(12);
	(*Z).push_back(8);
	(*p).push_back(1);
	(*p).push_back(10);
	(*p).push_back(2);
	(*p).push_back(20);
	(*p).push_back(30);
	s.add_layer(* (new Monolayer((*Z),(*p), 2.2, 1.5e-5, 3)));
	
	
	Illumination il;
	Detector det;
	// for (int i=0; i<s.layer_vec[0].(*Z)[i]; i++)
		// cout << s.layer_vec[0].(*Z)[i] << " ";
	// cout << endl;
	Spectrum spec(s, il, omega, det);
	// Spectrum spec(1e4, c, omega, det);
	// Xrf xrf(1e4, c, omega);
	// xrf.show();
	
	ofstream fout;
	fout.open("out.txt");
	for (int i=0; i<spec.y_vec.size(); i++)
		fout << spec.y_vec[i] << " ";
	fout << endl;
	spec.out(fout);
	fout.close();
	
	spec.show();
	cout << spec.xrf.lines.size() << endl;

	return 1;
}