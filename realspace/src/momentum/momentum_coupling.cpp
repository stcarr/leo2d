
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

void load_fftw_complex(std::vector< std::vector< std::vector< std::vector< fftw_complex* > > > >&out,
											 std::vector< std::vector< std::vector< std::vector< std::vector<double> > > > >&disps,
											 std::string file){

	std::ifstream fin(file.c_str());


	int num_sheets;
	std::vector<int> num_orbs;
	fin >> num_sheets;

	out.resize(num_sheets);
	disps.resize(num_sheets);
	num_orbs.resize(num_sheets);


	for (int s = 0; s < num_sheets; ++s){
		out[s].resize(num_sheets);
		disps[s].resize(num_sheets);
		fin >> num_orbs[s];
	}

	for (int s = 0; s < num_sheets-1; ++s){
		for (int direction = 0; direction < 2; ++direction){
			// direction = 0 -> lower sheet coupling to upper sheet
			// direction = 1 -> upper sheet coupling to lower sheet
			int s1, s2;
			if (direction == 0){
				s1 = s;
				s2 = s+1;
			} else {
				s1 = s+1;
				s2 = s;
			}

			int num_orb_1 = num_orbs[s1];
			int num_orb_2 = num_orbs[s2];

			disps[s1][s2].resize(num_orb_1);

			std::vector< std::vector<fftw_complex*> > temp_o1_data;

			for (int o1 = 0; o1 < num_orb_1; ++o1){

				disps[s1][s2][o1].resize(num_orb_2);
				std::vector<fftw_complex*> temp_o2_data;

				for (int o2 = 0; o2 < num_orb_2; ++o2){

					int file_s1;
					int file_s2;

					int file_o1;
					int file_o2;

					double orb_disp_x;
					double orb_disp_y;

					int x_size;
					int y_size;

					int file_type;

					fin >> file_s1;
					fin >> file_s2;
					fin >> file_o1;
					fin >> file_o2;
					fin >> orb_disp_x;
					fin >> orb_disp_y;
					fin >> x_size;
					fin >> y_size;
					fin >> file_type;

					/*
					printf("expect: [%d, %d, %d, %d, ?, ?, ?, ?, 0] \n",s1,s2,o1,o2);
					printf("file: [%d, %d, %d, %d, %lf, %lf, %d, %d, %d] \n",
									file_s1,file_s2,file_o1,file_o2, orb_disp_x, orb_disp_y, x_size, y_size, file_type);
					*/

					if (file_s1 != s1){
						printf("Warning! sheet1 mismatch in *fft.dat file \n");
						break;
					}
					if (file_s2 != s2){
						printf("Warning! sheet2 mismatch in *fft.dat file \n");
						break;
					}
					if (file_o1 != o1){
						printf("Warning! orbit1 mismatch in *fft.dat file \n");
						break;
					}
					if (file_o2 != o2){
						printf("Warning! orbit2 mismatch in *fft.dat file \n");
						break;
					}
					if (file_type != 0){
						printf("Warning! value type mismatch in *fft.dat file \n");
						break;
					}

					disps[s1][s2][o1][o2].resize(2);
					disps[s1][s2][o1][o2][0] = orb_disp_x;
					disps[s1][s2][o1][o2][1] = orb_disp_y;

					// debugging print statement
					// printf("load: [x_size, y_size] = [%d, %d] \n",x_size, y_size);


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

					fin >> file_s1;
					fin >> file_s2;
					fin >> file_o1;
					fin >> file_o2;
					fin >> orb_disp_x;
					fin >> orb_disp_y;
					fin >> x_size;
					fin >> y_size;
					fin >> file_type;

					/*
					printf("expect: [%d, %d, %d, %d, ?, ?, ?, ?, 1] \n",s1,s2,o1,o2);
					printf("file: [%d, %d, %d, %d, %lf, %lf, %d, %d, %d] \n",
									file_s1,file_s2,file_o1,file_o2, orb_disp_x, orb_disp_y, x_size, y_size, file_type);
					*/

					if (file_s1 != s1){
						printf("Warning! sheet1 mismatch in *fft.dat file \n");
						break;
					}
					if (file_s2 != s2){
						printf("Warning! sheet2 mismatch in *fft.dat file \n");
						break;
					}
					if (file_o1 != o1){
						printf("Warning! orbit1 mismatch in *fft.dat file \n");
						break;
					}
					if (file_o2 != o2){
						printf("Warning! orbit2 mismatch in *fft.dat file \n");
						break;
					}
					if (file_type != 1){
						printf("Warning! value type mismatch in *fft.dat file \n");
						break;
					}

					for (int i = 0; i < 2*x_size; i++)
						for (int j = 0; j < y_size; j++)
							fin >> temp_out[j + i*y_size][1];

					temp_o2_data.push_back(temp_out);


				} // end of o2 loop

				temp_o1_data.push_back(temp_o2_data);

			} // end of o1 loop

			out[s1][s2] = temp_o1_data;

		} // end of direction loop
	} // end of sheet loop

	fin.close();
}

double interp_4point(double x, double y, double v1, double v2, double v3, double v4) {
   double value = v1*(1-x)*(1-y) + v2*x*(1-y)+ v3*(1-x)*y + v4*x*y;
   return value;
}

double Momentum_coupling::interp_fft(double x_input, double y_input, int s1, int s2, int o1, int o2, int entry) {

	//printf("x_input = %lf, y_input = %lf, s1 = %d, s2 = %d, o1 = %d, o2 = %d, entry = %d \n",x_input, y_input,s1,s2,o1,o2,entry);
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

		double value = interp_4point(	x-x_int,
																	y-y_int,
																	fftw_data[s1][s2][o1][o2][y_int + x_int*y_s][entry],
																	fftw_data[s1][s2][o1][o2][y_int + (x_int+1)*y_s][entry],
																	fftw_data[s1][s2][o1][o2][y_int+1 + x_int*y_s][entry],
																	fftw_data[s1][s2][o1][o2][y_int+1+(x_int+1)*y_s][entry]
																	);

		if (r_sq < momentum_inner_cutoff*momentum_inner_cutoff){
			return sign*value;
		} else {
			double r_sq_diff = r_sq - momentum_inner_cutoff*momentum_inner_cutoff;
			return sign*value*exp(-2.0*r_sq_diff);
		}
	}

}

// verbose version of above call, for debugging purposes
double Momentum_coupling::interp_fft_v(double x_input, double y_input, int s1, int s2, int o1, int o2, int entry) {

	// (x_input, y_input) is location at which you want to know the (interpolated) value of the fftw_complex data
	// o1, o2 are two two orbitals whose interlayer coupling you are computing
	// entry is 0 for real and 1 for imaginary part

	fftw_complex* data;
	data = fftw_data[s1][s2][o1][o2];

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

double Momentum_coupling::get_fft_orb_disps(int s1, int s2, int o1, int o2, int dir){
	return fftw_disps[s1][s2][o1][o2][dir];
}


void Momentum_coupling::fft_setup(int L_x_in, int L_y_in, int length_x_in, int length_y_in, std::string fftw_file_in){

	L_x = L_x_in;
	L_y = L_y_in;
	length_x = length_x_in;
	length_y = length_y_in;
	load_fftw_complex(fftw_data, fftw_disps, fftw_file_in);
}
