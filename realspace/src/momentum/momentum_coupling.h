#ifndef MOMENTUM_COUPLING_HEADER
#define MOMENTUM_COUPLING_HEADER

#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void load_fftw_complex(	std::vector< std::vector< std::vector< std::vector< fftw_complex* > > > >&out,
											 	std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > >&disps,
											 	std::string file);
double interp_4point(double x, double y, double v1, double v2, double v3, double v4);

class Momentum_coupling{

	private:

		// fftw global info, to prevent having to reload *fftw.dat file for every call...
		// for coupling between o1 and o2 from sheet s1 to sheet s2, at x,y of type (real/cpx) call:
		//  >> fftw_data[s1][s2][o1][o2][y + x*y_s][type]
		// where y_s = length_y*L_y, is the y-dimension stride in the 1D vectorization of 2D FFTW output
		std::vector<std::vector< std::vector< std::vector< fftw_complex* > > > > fftw_data;
		std::vector<std::vector< std::vector< std::vector< std::vector<double> > > > > fftw_disps;
		int L_x;
		int L_y;
		int length_x;
		int length_y;

	public:

		Momentum_coupling();
		~Momentum_coupling();

		double interp_fft(double x_input, double y_input, int s1, int s2, int o1, int o2, int entry);
		double interp_fft_v(double x_input, double y_input, int s1, int s2, int o1, int o2, int entry);
		double get_fft_orb_disps(int s1, int s2, int o1, int o2, int dir);
		void fft_setup(int L_x_in, int L_y_in, int length_x_in, int length_y_in, std::string fftw_file_in);

};

#endif
