
#include "momentum_coupling.h"
#include "materials/materials.h"
#include <cmath>
#include <stdio.h>
// assumes sheet 1 to sheet 2 orientation, and sheet 2 is rotated by theta

const double pi2o3 	= M_PI*2/3;
const double pi1o3 	= M_PI/3;
const double pi2	= M_PI/2;
const double pi6 	= M_PI/6;
const double r_cut_graphene = 8;
const double r_cut2_graphene = 7;

// assuming (x,y) aims from sheet 1 to sheet 2 site
// theta rotates sheet counter-clockwise

Momentum_coupling::Momentum_coupling() {
}

Momentum_coupling::~Momentum_coupling() {
}

void load_fftw_complex(std::vector< std::vector<fftw_complex*> > &out, std::string file) {

	std::ifstream fin(file.c_str());
	int num_orb_1;
	int num_orb_2;

	fin >> num_orb_1;
	fin >> num_orb_2;

	for (int o1 = 0; o1 < num_orb_1; ++o1){

		std::vector<fftw_complex*> temp_vec;

		for (int o2 = 0; o2 < num_orb_2; ++o2){

			int x_size;
			int y_size;

			int file_o1;
			int file_o2;

			int file_type;

			fin >> file_o1;
			fin >> file_o2;
			fin >> x_size;
			fin >> y_size;
			fin >> file_type;

			if (file_o1 != o1){
				printf("Warning! orbit mismatch in *fft.dat file \n");
				break;
			}
			if (file_o2 != o2){
				printf("Warning! orbit mismatch in *fft.dat file \n");
				break;
			}
			if (file_type != 0){
				printf("Warning! orbit mismatch in *fft.dat file \n");
				break;
			}

			// debugging print statement
			//printf("load: [x_size, y_size] = [%d, %d] \n",x_size, y_size);


			fftw_complex* temp_out;
			temp_out = (fftw_complex*) fftw_malloc(x_size*y_size*2*sizeof(fftw_complex));

			for (int i = 0; i < 2*x_size; i++)
				for (int j = 0; j < y_size; j++)
					fin >> temp_out[j + i*y_size][0];

			// following prints loaded matrix to terminal, for debugging
			/*
			for (int i = 0; i < 2*x_size; i++){
				for (int j = 0; j < y_size; j++){
					printf("%lf, ",temp_out[j + i*y_size][0]);
				}
				printf(" \n");
			}
			*/

			// Now get complex data for this orbit pairing

			fin >> file_o1;
			fin >> file_o2;
			fin >> x_size;
			fin >> y_size;
			fin >> file_type;

			if (file_o1 != o1){
				printf("Warning! orbit mismatch in *fft.dat file \n");
				break;
			}
			if (file_o2 != o2){
				printf("Warning! orbit mismatch in *fft.dat file \n");
				break;
			}
			if (file_type != 1){
				printf("Warning! orbit mismatch in *fft.dat file \n");
				break;
			}

			for (int i = 0; i < 2*x_size; i++)
				for (int j = 0; j < y_size; j++)
					fin >> temp_out[j + i*y_size][1];

			temp_vec.push_back(temp_out);


		}

		out.push_back(temp_vec);

	}

	fin.close();
}

double interp_4point(double x, double y, double v1, double v2, double v3, double v4) {
   double value = v1*(1-x)*(1-y) + v2*x*(1-y)+ v3*(1-x)*y + v4*x*y;
   return value;
}

double Momentum_coupling::interp_fft(double x_input, double y_input, int o1, int o2, int entry) {

	// (x_input, y_input) is location at which you want to know the (interpolated) value of the fftw_complex data
	// o1, o2 are two two orbitals whose interlayer coupling you are computing
	// entry is 0 for real and 1 for imaginary part

	double momentum_inner_cutoff = 5.0;
	double momentum_outer_cutoff = 4.0;

	double r_sq = (x_input)*(x_input) + (y_input)*(y_input);

	if (r_sq > momentum_outer_cutoff*momentum_outer_cutoff){
		return 0.0;
	} else {

		int x_s = length_x*L_x;
		int y_s = length_y*L_y;
		double x = x_input*L_x/M_PI+x_s;
		double y = y_input*L_y/M_PI;

		// Keep track of odd/even nature of discrete FT symmetry
		double sign = 1;

		if (y < 0){

			y = -y; // we use the y-symmetry in a FFT of purely real input data

			if (x >= 1){
				x =  2*x_s - x; // also need to flip the x-axis about its center (the first x row stays the same though)!
			}

			if (entry == 1){
				sign = -1; // complex part is odd, real part is even
			}


		}

		if (x < 0 || x > 2*x_s - 1 || y > y_s - 1){ // here we do "- 1" to prevent wrap-around errors (i.e. interpolating at [2*x_s,y] would sample [0,y+1] for right-hand points!!
			return 0;
		}

		int x_int = int(x);
		int y_int = int(y);

		double value = interp_4point(x-x_int,y-y_int, fftw_data[o1][o2][y_int + x_int*y_s][entry], fftw_data[o1][o2][y_int + (x_int+1)*y_s][entry], fftw_data[o1][o2][y_int+1 + x_int*y_s][entry],fftw_data[o1][o2][y_int+1+(x_int+1)*y_s][entry]);

		if (r_sq < momentum_inner_cutoff*momentum_inner_cutoff){
			return sign*value;
		} else {
			double r_sq_diff = r_sq - momentum_inner_cutoff*momentum_inner_cutoff;
			return sign*value*exp(-2.0*r_sq_diff);
		}
	}

}

// verbose version of above call, for debugging purposes
double Momentum_coupling::interp_fft_v(double x_input, double y_input, int o1, int o2, int entry) {

	// (x_input, y_input) is location at which you want to know the (interpolated) value of the fftw_complex data
	// o1, o2 are two two orbitals whose interlayer coupling you are computing
	// entry is 0 for real and 1 for imaginary part

	fftw_complex* data;
	data = fftw_data[o1][o2];

	int x_s = length_x*L_x;
	int y_s = length_y*L_y;
	double x = x_input*L_x/M_PI+x_s;
	double y = y_input*L_y/M_PI;

	printf("[x_s, y_s] = [%d, %d] \n",x_s,y_s);
	printf("[x,y] = [%lf, %lf] \n",x,y);

	if (y < 0)
		y = -y; // we use the y-symmetry in a FFT of purely real input data
	if (x < 0 || x > 2*x_s - 1 || y > y_s - 1) // here we do "- 1" to prevent wrap-around errors (i.e. interpolating at [2*x_s,y] would sample [0,y+1] for right-hand points!!
		return 0;



	int x_int = int(x);
	int y_int = int(y);

	printf("[x_int, y_int] = [%d, %d] \n",x_int, y_int);
	printf("[data] = [%lf, %lf, %lf, %lf] \n",data[y_int + x_int*y_s][entry], data[y_int + (x_int+1)*y_s][entry], data[y_int+1 + x_int*y_s][entry],data[y_int+1+(x_int+1)*y_s][entry]);

	double value = interp_4point(x-x_int,y-y_int, data[y_int + x_int*y_s][entry], data[y_int + (x_int+1)*y_s][entry], data[y_int+1 + x_int*y_s][entry],data[y_int+1+(x_int+1)*y_s][entry]);
	return value;

}

void Momentum_coupling::fft_setup(int L_x_in, int L_y_in, int length_x_in, int length_y_in, std::string fftw_file_in){

	L_x = L_x_in;
	L_y = L_y_in;
	length_x = length_x_in;
	length_y = length_y_in;

	load_fftw_complex(fftw_data, fftw_file_in);
}
