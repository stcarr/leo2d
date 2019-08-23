/*
 * File:   hstruct.cpp
 * Author: Stephen Carr
 *
 * Created on January 13, 2016, 3:16 PM
 */

#include "locality.h"
#include "tools/numbers.h"
#include "params/param_tools.h"
#include "params/job_params.h"
#include "materials/read_mat.h"

//#include "transport/ballistic.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <time.h>
#include <random>
#include <cmath>

//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>

using namespace numbers;

Locality::Locality(std::vector<Sdata> sdata_in,std::vector<double> heights_in,std::vector<double> angles_in) {
    
    // Set all of the geometry options for matrix construction
    
    job_name = "HSTRUCT_TEST";
    
    sdata = sdata_in;
    heights = heights_in;
    angles = angles_in;
    
}

Locality::Locality(const Locality& orig) {
    
}

Locality::~Locality() {
    
}

void Locality::setup(Job_params opts_in){
    
    // Controls run-specific parameters of the solver methods
    opts = opts_in;
    job_name = opts.getString("job_name");
    int mat_from_file = opts.getInt("mat_from_file");
    
    if (rank == root){
        opts.printParams();
        if (mat_from_file == 1){
            loadedMatData = ReadMat::loadMat("test_data");
            loadedMatData.MPI_Bcast_root(root);
        }
    } else {
        if (mat_from_file == 1){
            // wait for command to read file
            loadedMatData.MPI_Bcast(root);
        }
    }
    
    int boundary_condition = opts.getInt("boundary_condition");
    
    // update supercell if BC
    if (boundary_condition == 1){
        Locality::setupSupercell();
    }
    
}

void Locality::setupSupercell(){
    
    int type = opts.getInt("supercell_type");
    int mat_from_file = opts.getInt("mat_from_file");
    
    // Standard fixed supercell definition
    if (type == 0) {
        std::vector< std::vector<double> > sc = opts.getDoubleMat("supercell");
        for (int i = 0; i < sdata.size(); ++i){
            double theta = angles[i];
            std::vector< std::vector<double> > sc_here;
            sc_here.resize(2);
            sc_here[0].resize(2);
            sc_here[1].resize(2);
            
            sc_here[0][0] =  sc[0][0]*cos(theta) + sc[0][1]*sin(theta);
            sc_here[0][1] = -sc[0][0]*sin(theta) + sc[0][1]*cos(theta);
            sc_here[1][0] =  sc[1][0]*cos(theta) + sc[1][1]*sin(theta);
            sc_here[1][1] = -sc[1][0]*sin(theta) + sc[1][1]*cos(theta);
            sdata[i].supercell = sc_here;
            
        }
    } else if (type == 1) { // (M,N) Supercell type, in the /src/matrix_gen.cpp file
        std::vector<int> sc_groups;
        
        if ((int)sdata.size() < 3){
            sc_groups.resize(2);
            sc_groups[0] = 0;
            sc_groups[1] = 1;
            opts.setParam("sc_groups",sc_groups);
        } else {
            sc_groups = opts.getIntVec("sc_groups");
            
            printf("sc_groups = ");
            for (int i = 0; i < sc_groups.size(); ++i){
                printf("%d , ",sc_groups[i]);
            }
            printf("\n");
        }
        
        int M = opts.getInt("m_supercell");
        int N = opts.getInt("n_supercell");
        
        double theta = acos((N*N + 4*N*M + M*M)/(2.0*(N*N + N*M + M*M)));
        if (M < N){
            theta = -theta;
        }
        printf("supercell theta = %lf degrees (acos(%lf) )\n",360.0*theta/(2.0*M_PI), (N*N + 4*N*M + M*M)/(2.0*(N*N + N*M + M*M)));
        // we assume unitCell is same for both sheets...
        for (int i = 0; i < sdata.size(); ++i){
            
            std::array<std::array<double, 2>, 2> base_unitCell;
            if (mat_from_file == 0){
                base_unitCell = Materials::lattice(sdata[0].mat);
            } else {
                base_unitCell = ReadMat::getLattice(loadedMatData, 0);
            }
            
            int A1_num_a1;
            int A1_num_a2;
            int A2_num_a1;
            int A2_num_a2;
            
            int sc_group_here = sc_groups[i];
            
            if (sc_group_here == 0){
                angles[i] = 0;
                A1_num_a1 = N;
                A1_num_a2 = M;
                A2_num_a1 = -M;
                A2_num_a2 = (M+N);
            } else if (sc_group_here == 1){
                angles[i] = theta;
                A1_num_a1 = M;
                A1_num_a2 = N;
                A2_num_a1 = -N;
                A2_num_a2 = (M+N);
            }
            
            std::vector< std::vector<double> > unitCell;
            unitCell.resize(2);
            unitCell[0].resize(2);
            unitCell[1].resize(2);
            
            // Not updating the twist-angle is correct, as this will make sure the supercell is always
            // saved in coordinates relative to the individual sheet, not the overall heterostructure.
            double theta_here = 0.0;
            
            unitCell[0][0] = cos(theta_here)*base_unitCell[0][0] - sin(theta_here)*base_unitCell[0][1];
            unitCell[0][1] = sin(theta_here)*base_unitCell[0][0] + cos(theta_here)*base_unitCell[0][1];
            unitCell[1][0] = cos(theta_here)*base_unitCell[1][0] - sin(theta_here)*base_unitCell[1][1];
            unitCell[1][1] = sin(theta_here)*base_unitCell[1][0] + cos(theta_here)*base_unitCell[1][1];
            
            std::vector< std::vector<double> > sc_here;
            sc_here.resize(2);
            sc_here[0].resize(2);
            sc_here[1].resize(2);
            
            sc_here[0][0] = A1_num_a1*unitCell[0][0] + A1_num_a2*unitCell[1][0];
            sc_here[0][1] = A1_num_a1*unitCell[0][1] + A1_num_a2*unitCell[1][1];
            sc_here[1][0] = A2_num_a1*unitCell[0][0] + A2_num_a2*unitCell[1][0];
            sc_here[1][1] = A2_num_a1*unitCell[0][1] + A2_num_a2*unitCell[1][1];
            
            std::vector< std::vector<int> > sc_stride_here;
            sc_stride_here.resize(2);
            sc_stride_here[0].resize(2);
            sc_stride_here[1].resize(2);
            
            sc_stride_here[0][0] = A1_num_a1;
            sc_stride_here[0][1] = A1_num_a2;
            sc_stride_here[1][0] = A2_num_a1;
            sc_stride_here[1][1] = A2_num_a2;
            
            /*
             std::vector< std::vector<double> > local_sc_here;
             local_sc_here.resize(2);
             local_sc_here[0].resize(2);
             local_sc_here[1].resize(2);
             
             double theta_here = -angles[i];
             
             local_sc_here[0][0] = cos(theta_here)*sc_here[0][0] - sin(theta_here)*sc_here[0][1];
             local_sc_here[0][1] = sin(theta_here)*sc_here[0][0] + cos(theta_here)*sc_here[0][1];
             local_sc_here[1][0] = cos(theta_here)*sc_here[1][0] - sin(theta_here)*sc_here[1][1];
             local_sc_here[1][1] = sin(theta_here)*sc_here[1][0] + cos(theta_here)*sc_here[1][1];
             */
            printf("unitCell  = [%lf %lf; %lf %lf]\n",unitCell[0][0],unitCell[0][1],unitCell[1][0],unitCell[1][1]);
            printf("supercell = [%lf %lf; %lf %lf]\n", sc_here[0][0], sc_here[0][1], sc_here[1][0], sc_here[1][1]);
            //printf("local_supercell = [%lf %lf; %lf %lf]\n", local_sc_here[0][0], local_sc_here[0][1], local_sc_here[1][0], local_sc_here[1][1]);
            printf("sc_stride = [%d %d; %d %d]\n", sc_stride_here[0][0], sc_stride_here[0][1], sc_stride_here[1][0], sc_stride_here[1][1]);
            
            
            sdata[i].supercell = sc_here;
            sdata[i].supercell_stride = sc_stride_here;
            
            // Since sheet[0] has 0 twist, we can use it's supercell as the supercell for hstruct!!
            if (i == 0){
                opts.setParam("supercell",sc_here);
                
                int matrix_pos_save = opts.getInt("matrix_pos_save");
                if (matrix_pos_save == 1){
                    
                    string job_name = opts.getString("job_name");
                    const char* extension = "_supercell.dat";
                    
                    ofstream outFile;
                    outFile.open ((job_name + extension).c_str());
                    outFile <<     "0, 0, 0, " << sc_here[0][0] << ", " << sc_here[0][1] << ", 0\n" <<
                    "0, 0, 0, " << sc_here[1][0] << ", " << sc_here[1][1] << ", 0\n" <<
                    sc_here[0][0] << ", " << sc_here[0][1] << ", 0, " <<
                    sc_here[0][0]+sc_here[1][0] << ", " << sc_here[0][1]+sc_here[1][1] << ", 0\n" <<
                    sc_here[1][0] << ", " << sc_here[1][1] << ", 0, " <<
                    sc_here[0][0]+sc_here[1][0] << ", " << sc_here[0][1]+sc_here[1][1] << ", 0\n";
                    outFile.close();
                }
            }
            
        }
        
    } else if (type == 2) { // (X,Y) Grid type
        
        int X = opts.getInt("x_supercell");
        int Y = opts.getInt("y_supercell");
        
        std::array<std::array<double, 2>, 2> unitCell;
        if (mat_from_file == 0){
            unitCell = Materials::lattice(sdata[0].mat);
        } else {
            unitCell = ReadMat::getLattice(loadedMatData, 0);
        }
        
        std::vector< std::vector<double> > sc_here;
        sc_here.resize(2);
        sc_here[0].resize(2);
        sc_here[1].resize(2);
        
        sc_here[0][0] = X*unitCell[0][0];
        sc_here[0][1] = X*unitCell[0][1];
        sc_here[1][0] = Y*unitCell[1][0];
        sc_here[1][1] = Y*unitCell[1][1];
        
        std::vector< std::vector<int> > sc_stride_here;
        sc_stride_here.resize(2);
        sc_stride_here[0].resize(2);
        sc_stride_here[1].resize(2);
        
        sc_stride_here[0][0] = X;
        sc_stride_here[0][1] = 0;
        sc_stride_here[1][0] = 0;
        sc_stride_here[1][1] = Y;
        
        for (int i = 0; i < (int)sdata.size(); ++i){
            
            double theta = angles[i];
            std::vector< std::vector<double> > sc_rot;
            sc_rot.resize(2);
            sc_rot[0].resize(2);
            sc_rot[1].resize(2);
            
            sc_rot[0][0] =  sc_here[0][0]*cos(theta) + sc_here[0][1]*sin(theta);
            sc_rot[0][1] = -sc_here[0][0]*sin(theta) + sc_here[0][1]*cos(theta);
            sc_rot[1][0] =  sc_here[1][0]*cos(theta) + sc_here[1][1]*sin(theta);
            sc_rot[1][1] = -sc_here[1][0]*sin(theta) + sc_here[1][1]*cos(theta);
            
            std::vector< std::vector<int> > sc_stride_rot;
            sc_stride_rot.resize(2);
            sc_stride_rot[0].resize(2);
            sc_stride_rot[1].resize(2);
            
            sc_stride_rot[0][0] =  sc_stride_here[0][0]*cos(theta) + sc_stride_here[0][1]*sin(theta);
            sc_stride_rot[0][1] = -sc_stride_here[0][0]*sin(theta) + sc_stride_here[0][1]*cos(theta);
            sc_stride_rot[1][0] =  sc_stride_here[1][0]*cos(theta) + sc_stride_here[1][1]*sin(theta);
            sc_stride_rot[1][1] = -sc_stride_here[1][0]*sin(theta) + sc_stride_here[1][1]*cos(theta);
            
            //printf("unitCell  %d = [%lf %lf; %lf %lf]\n",i, unitCell[0][0],unitCell[0][1],unitCell[1][0],unitCell[1][1]);
            printf("supercell %d = [%lf %lf; %lf %lf]\n", i, sc_rot[0][0], sc_rot[0][1], sc_rot[1][0], sc_rot[1][1]);
            printf("sc_stride %d = [%d %d; %d %d]\n", i, sc_stride_rot[0][0], sc_stride_rot[0][1], sc_stride_rot[1][0], sc_stride_rot[1][1]);
            
            sdata[i].supercell = sc_rot;
            sdata[i].supercell = sc_rot;
            sdata[i].supercell_stride = sc_stride_here;
        }
        
        // Since sheet[0] has 0 twist, we can use it's supercell as the supercell for hstruct!!
        opts.setParam("supercell",sc_here);
        
        int matrix_pos_save = opts.getInt("matrix_pos_save");
        if (matrix_pos_save == 1){
            
            string job_name = opts.getString("job_name");
            const char* extension = "_supercell.dat";
            
            ofstream outFile;
            outFile.open ((job_name + extension).c_str());
            outFile <<     "0, 0, 0, " << sc_here[0][0] << ", " << sc_here[0][1] << ", 0\n" <<
            "0, 0, 0, " << sc_here[1][0] << ", " << sc_here[1][1] << ", 0\n" <<
            sc_here[0][0] << ", " << sc_here[0][1] << ", 0, " <<
            sc_here[0][0]+sc_here[1][0] << ", " << sc_here[0][1]+sc_here[1][1] << ", 0\n" <<
            sc_here[1][0] << ", " << sc_here[1][1] << ", 0, " <<
            sc_here[0][0]+sc_here[1][0] << ", " << sc_here[0][1]+sc_here[1][1] << ", 0\n";
            outFile.close();
        }
        
    } else if (type == 3) { // (M,N) type, trilayer, helical/chiral case
        std::vector<int> sc_groups;
        
        if ((int)sdata.size() < 4){
            sc_groups.resize(3); // matches the number of layers
            sc_groups[0] = 0;
            sc_groups[1] = 1;
            sc_groups[2] = 2;
            opts.setParam("sc_groups",sc_groups);
        } else {
            sc_groups = opts.getIntVec("sc_groups");
            
            printf("sc_groups = ");
            for (int i = 0; i < sc_groups.size(); ++i){
                printf("%d , ",sc_groups[i]);
            }
            printf("\n");
        }
        
        int M = opts.getInt("m_supercell");
        int N = opts.getInt("n_supercell");
        
        double theta = acos((N*N + 4*N*M + M*M)/(2.0*(N*N + N*M + M*M)));
        if (M < N){
            theta = -theta;
        }
        printf("supercell theta = %lf degrees (acos(%lf) )\n",360.0*theta/(2.0*M_PI), (N*N + 4*N*M + M*M)/(2.0*(N*N + N*M + M*M)));
        // we assume unitCell is same for both sheets...
        for (int i = 0; i < sdata.size(); ++i){
            
            std::array<std::array<double, 2>, 2> base_unitCell;
            if (mat_from_file == 0){
                base_unitCell = Materials::lattice(sdata[0].mat);
            } else {
                base_unitCell = ReadMat::getLattice(loadedMatData, 0);
            }
            
            int A1_num_a1;
            int A1_num_a2;
            int A2_num_a1;
            int A2_num_a2; // so far same as bilayer
            
            int sc_group_here = sc_groups[i];
            
            if (sc_group_here == 0){
                angles[i] = -theta;
                A1_num_a1 = pow(N,2) - pow(M,2);
                A1_num_a2 = pow(M,2) + 2*M*N;
                A2_num_a1 = -(pow(M,2) + 2*M*N);
                A2_num_a2 = pow(N,2) + 2*M*N;
            } else if (sc_group_here == 1){
                angles[i] = 0;
                A1_num_a1 = 0;
                A1_num_a2 = pow(M,2) + (M*N) + pow(N,2);
                A2_num_a1 = -(pow(N,2) + (M*N) + pow(M,2));
                A2_num_a2 = pow(N,2)+ M*N + pow(M,2);
            } else if (sc_group_here == 2){
                angles[i] = theta;
                A1_num_a1 = pow(M,2) - pow(N,2);
                A1_num_a2 = pow(N,2) + 2*M*N;
                A2_num_a1 = -(pow(N,2) + (2*M*N));
                A2_num_a2 = pow(M,2) + (2*M*N);
            }
            
            std::vector< std::vector<double> > unitCell;
            unitCell.resize(2);
            unitCell[0].resize(2);
            unitCell[1].resize(2);
            
            // Not updating the twist-angle is correct, as this will make sure the supercell is always
            // saved in coordinates relative to the individual sheet, not the overall heterostructure.
            double theta_here = 0.0;
            
            unitCell[0][0] = cos(theta_here)*base_unitCell[0][0] - sin(theta_here)*base_unitCell[0][1];
            unitCell[0][1] = sin(theta_here)*base_unitCell[0][0] + cos(theta_here)*base_unitCell[0][1];
            unitCell[1][0] = cos(theta_here)*base_unitCell[1][0] - sin(theta_here)*base_unitCell[1][1];
            unitCell[1][1] = sin(theta_here)*base_unitCell[1][0] + cos(theta_here)*base_unitCell[1][1];
            
            std::vector< std::vector<double> > sc_here;
            sc_here.resize(2);
            sc_here[0].resize(2);
            sc_here[1].resize(2);
            
            sc_here[0][0] = A1_num_a1*unitCell[0][0] + A1_num_a2*unitCell[1][0];
            sc_here[0][1] = A1_num_a1*unitCell[0][1] + A1_num_a2*unitCell[1][1];
            sc_here[1][0] = A2_num_a1*unitCell[0][0] + A2_num_a2*unitCell[1][0];
            sc_here[1][1] = A2_num_a1*unitCell[0][1] + A2_num_a2*unitCell[1][1];
            
            std::vector< std::vector<int> > sc_stride_here;
            sc_stride_here.resize(2);
            sc_stride_here[0].resize(2);
            sc_stride_here[1].resize(2);
            
            sc_stride_here[0][0] = A1_num_a1;
            sc_stride_here[0][1] = A1_num_a2;
            sc_stride_here[1][0] = A2_num_a1;
            sc_stride_here[1][1] = A2_num_a2;
            
            printf("unitCell  = [%lf %lf; %lf %lf]\n",unitCell[0][0],unitCell[0][1],unitCell[1][0],unitCell[1][1]);
            printf("supercell = [%lf %lf; %lf %lf]\n", sc_here[0][0], sc_here[0][1], sc_here[1][0], sc_here[1][1]);
            
            printf("sc_stride = [%d %d; %d %d]\n", sc_stride_here[0][0], sc_stride_here[0][1], sc_stride_here[1][0], sc_stride_here[1][1]);
            
            
            sdata[i].supercell = sc_here;
            sdata[i].supercell_stride = sc_stride_here;
            
            // Since sheet[0] has 0 twist, we can use it's supercell as the supercell for hstruct!!
            if (i == 0){
                opts.setParam("supercell",sc_here);
                
                int matrix_pos_save = opts.getInt("matrix_pos_save");
                if (matrix_pos_save == 1){
                    
                    string job_name = opts.getString("job_name");
                    const char* extension = "_supercell.dat";
                    
                    ofstream outFile;
                    outFile.open ((job_name + extension).c_str());
                    outFile <<     "0, 0, 0, " << sc_here[0][0] << ", " << sc_here[0][1] << ", 0\n" <<
                    "0, 0, 0, " << sc_here[1][0] << ", " << sc_here[1][1] << ", 0\n" <<
                    sc_here[0][0] << ", " << sc_here[0][1] << ", 0, " <<
                    sc_here[0][0]+sc_here[1][0] << ", " << sc_here[0][1]+sc_here[1][1] << ", 0\n" <<
                    sc_here[1][0] << ", " << sc_here[1][1] << ", 0, " <<
                    sc_here[0][0]+sc_here[1][0] << ", " << sc_here[0][1]+sc_here[1][1] << ", 0\n";
                    outFile.close();
                }
            }
            
        }
        
    }
    
}


void Locality::getVacanciesFromFile(std::vector<std::vector<int> > &v, std::vector<std::vector<int> > &t, Job_params opts_in){
    
    int solver_type = opts_in.getInt("solver_type");
    std::string vacancy_file = opts.getString("vacancy_file");
    
    // MLMC method
    if (solver_type == 3){
        
        std::string line;
        std::ifstream in_file;
        in_file.open(vacancy_file.c_str());
        
        int jobID = -1;
        int clusterID = -1;
        std::vector<int> temp_v;
        std::vector<int> temp_ids;
        
        if (in_file.is_open())
        {
            while ( getline(in_file,line) ) {
                
                std::istringstream in_line(line);
                std::string in_string;
                
                while ( getline(in_line, in_string, ' ') )    {
                    
                    if (in_string == "JOBID"){
                        
                        // push back the previous job if this is not the first one
                        if (jobID != -1){
                            //printf("jobID = %d, clusterID = %d \n",jobID,clusterID);
                            temp_ids.push_back(jobID);
                            temp_ids.push_back(clusterID);
                            t.push_back(temp_ids);
                            v.push_back(temp_v);
                            
                            temp_ids.clear();
                            temp_v.clear();
                        }
                        
                        // get jobID
                        getline(in_line,in_string,' ');
                        getline(in_line,in_string,' ');
                        jobID = atoi(in_string.c_str());
                        
                    }
                    
                    if (in_string == "CLUSTERID"){
                        
                        // get clusterID
                        getline(in_line,in_string,' ');
                        getline(in_line,in_string,' ');
                        clusterID = atoi(in_string.c_str());
                    }
                    
                    if (in_string == "VAC"){
                        
                        // gets the "=" sign
                        getline(in_line,in_string,' ');
                        
                        // gets all vacancies (NO SPACES!)
                        std::string v_line;
                        getline(in_line,v_line,' ');
                        std::istringstream in_v_line(v_line);
                        
                        while ( getline(in_v_line, in_string, ',') )    {
                            temp_v.push_back(atoi(in_string.c_str()) - 1);
                        }
                        
                    }
                    
                }
            }
        }
        
        // and make sure to push back the last job!
        temp_ids.push_back(jobID);
        temp_ids.push_back(clusterID);
        t.push_back(temp_ids);
        v.push_back(temp_v);
    }
    
    if (solver_type == 4) {
        
        std::string line;
        std::ifstream in_file;
        in_file.open("vacancies.dat");
        if (in_file.is_open())
        {
            while ( getline(in_file,line) ) {
                
                std::vector<int> temp_v;
                std::vector<int> temp_t;
                
                std::istringstream in_line(line);
                
                std::string in_string;
                
                while ( getline(in_line, in_string, '=') )    {
                    
                    if (in_string == "JOB"){
                        
                        std::string v_line;
                        std::string t_line;
                        
                        getline(in_file,v_line);
                        getline(in_file,t_line);
                        
                        std::istringstream in_v_line(v_line);
                        std::istringstream in_t_line(t_line);
                        
                        while ( getline(in_v_line, in_string, ',') )    {
                            temp_v.push_back(atoi(in_string.c_str()) - 1);
                        }
                        while ( getline(in_t_line, in_string, ',') )    {
                            temp_t.push_back(atoi(in_string.c_str()) - 1);
                        }
                        
                        v.push_back(temp_v);
                        t.push_back(temp_t);
                        
                    }
                    
                }
            }
        }
    }
}

double Locality::getLocalTheta(int k1, int* index_to_grid, double* index_to_pos, int max){
    
    
    // we don't have access to the grid_to_index array anymore
    // so we have to try to find neighbours by guessing their index:
    
    // first, the starting info
    int i1 = index_to_grid[k1*4 + 0];
    int j1 = index_to_grid[k1*4 + 1];
    int o1 = index_to_grid[k1*4 + 2];
    int s1 = index_to_grid[k1*4 + 3];
    
    double x1 = index_to_pos[k1*3 + 0];
    double y1 = index_to_pos[k1*3 + 1];
    double z1 = index_to_pos[k1*3 + 2];
    
    // now loop over 4 possible nearest neighbors
    // these are in the same "j" column of the tight-binding model)
    int k2, i2, j2, o2, s2;
    int found_neighbor = 0;
    
    for (int dk = -3; dk < 4; dk = dk + 2) {
        int temp_k = k1 + dk;
        if (temp_k >= 0 && temp_k < max){
            int temp_i = index_to_grid[temp_k*4 + 0];
            int temp_j = index_to_grid[temp_k*4 + 1];
            int temp_o = index_to_grid[temp_k*4 + 2];
            int temp_s = index_to_grid[temp_k*4 + 3];
            
            if (isNearestNeighbor(i1, j1, o1, s1, temp_i, temp_j, temp_o, temp_s)){
                k2 = temp_k;
                i2 = temp_i;
                j2 = temp_j;
                o2 = temp_o;
                s2 = temp_s;
                found_neighbor = 1;
            }
        }
        
        if (found_neighbor == 1)
            break;
    }
    
    double theta;
    if (found_neighbor == 1){
        double x2 = index_to_pos[k2*3 + 0];
        double y2 = index_to_pos[k2*3 + 1];
        double z2 = index_to_pos[k2*3 + 2];
        
        double dx = x2 - x1;
        double dy = y2 - y1;
        
        if (o1 == 0){
            theta = atan2(dy, dx) - PI_6;
        } else if (o1 == 1){
            theta = atan2(dy, dx) - PI_2;
        }
        //printf("found neighbour with r = %lf (%lf, %lf) \n",sqrt(dx*dx + dy*dy),dx,dy);
        
    } else {
        // we couldn't find a nearest neighbor, so just throw the global twist...
        // (should only occurs on a few edge orbitals, so not very important)
        theta = angles[s1];
    }
    
    return theta;
    
}

int Locality::isNearestNeighbor(int i1, int j1, int o1, int s1, int i2, int j2, int o2, int s2){
    
    // have to be on same sheet
    if (s1 == s2){
        // have to have conjugate orbs
        if ((o1 == 0 && o2 == 1) || (o1 == 1 && o2 ==0)){
            
            // hardcoded for graphene... should use Material::(type) function call instead...
            // in same unit cell, or one of the nearby unit-cells with right orbital coupling
            if (i1 == i2 && j1 == j2){
                return 1;
            } else if (o1 == 0){
                if ( (i2 - i1) == -1 && (j2 - j1) == 0 )
                    return 1;
                
                if ( (i2 - i1) == 0  && (j2 - j1) == -1 )
                    return 1;
                
            } else if (o1 == 1){
                if ( (i2 - i1) == 1  && (j2 - j1) == 0 )
                    return 1;
                
                if ( (i2 - i1) == 0  && (j2 - j1) == 1 )
                    return 1;
            }
            
        }
    }
    
    return 0;
    
}

int Locality::initMPI(int argc, char** argv){
    
    // Start  MPI on each rank
    
    MPI_Init(&argc,&argv);
    
    root = 0;
    print_rank = 1;
    
    size = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    
    if (size == 1)
        return -1;
    else
        return 0;
}

void Locality::constructGeom(){
    
    // The root node (rank == 0) uses the hstruct and sheet objects to figure out the geometery and tight-binding model sparse matrix structure
    // This information is then pased to all worker nodes
    time(&constructStart);
    
    int solver_type = opts.getInt("solver_type");
    int solver_space = opts.getInt("solver_space");
    int boundary_condition = opts.getInt("boundary_condition");
    int strain_type = opts.getInt("strain_type");
    int fft_from_file = opts.getInt("fft_from_file");
    int nShifts = opts.getInt("nShifts");
    int mat_from_file = opts.getInt("mat_from_file");
    
    int max_pairs;
    
    int* inter_pairs_i;
    int* inter_pairs_j;
    int* inter_supercell_vecs_i;
    int* inter_supercell_vecs_j;
    
    int* intra_pairs_i;
    int* intra_pairs_j;
    double* intra_pairs_t;
    int* intra_supercell_vecs_i;
    int* intra_supercell_vecs_j;
    
    int* index_to_grid_i;
    int* index_to_grid_j;
    int* index_to_grid_l;
    int* index_to_grid_s;
    
    double* index_to_pos_x;
    double* index_to_pos_y;
    double* index_to_pos_z;
    
    std::vector< std::vector<double> > shift_configs;
    double* shift_configs_x;
    double* shift_configs_y;
    
    // vacancy and target indices
    
    std::vector<std::vector<int> > v_work;
    std::vector<std::vector<int> > target_indices;
    
    if (rank == root){
        
        //LoadedMat lmat = ReadMat::loadMat("test_data");
        
        // Build Hstruct object
        std::vector<Sheet> sheets;
        
        printf("building sheets\n");
        for (int i = 0; i < sdata.size(); ++i){
            if (mat_from_file == 0){
                sheets.push_back(Sheet(sdata[i]));
            } else {
                sheets.push_back( Sheet(sdata[i],loadedMatData,i) );
            }
        }
        
        printf("rank %d building Hstruct. \n", rank);
        Hstruct h(sheets,angles,heights,solver_space,opts);
        
        if (mat_from_file > 0){
            h.setLoadedMatData(loadedMatData);
        }
        
        // Broadcast "index to grid" mapping information
        max_index = h.getMaxIndex();
        
        // These are Twist or Strain solver types, so we only need to get targets
        if (solver_type == 1 || solver_type == 2 || solver_type == 5 || solver_type == 6){
            target_indices = h.getTargetList(opts);
        }
        
        // MLMC of vacancy defects, loads vacancies in job creation loop later
        if (solver_type == 3){
            
            // Old code, no longer needed as all vacancy loading in rootChebSolve() now
            /*
             target_indices = h.getTargetList(opts);
             std::vector<int> temp_v_work;
             temp_v_work.push_back(-1);
             v_work.push_back(temp_v_work);
             //v_work = h.getVacancyList(center_index[0],nShifts);
             */
        }
        
        // load vacancies from file
        if (solver_type == 4){
            getVacanciesFromFile(v_work, target_indices, opts);
        }
        
        
        MPI::COMM_WORLD.Bcast(&max_index, 1, MPI_INT, root);
        
        std::vector<std::vector<int> > index_vec = h.getIndexArray();
        
        index_to_grid_i = new int[max_index];
        index_to_grid_j = new int[max_index];
        index_to_grid_l = new int[max_index];
        index_to_grid_s = new int[max_index];
        
        for (int k = 0; k < max_index; ++k){
            
            index_to_grid_i[k] = index_vec[k][0];
            index_to_grid_j[k] = index_vec[k][1];
            index_to_grid_l[k] = index_vec[k][2];
            index_to_grid_s[k] = index_vec[k][3];
            
        }
        
        // Construct shift configs if using configuration strain (strain_type == 2 or 5)
        if (strain_type == 2 || strain_type == 5){
            
            printf("Building shift_configs. \n");
            h.getShiftConfigs(shift_configs, opts);
            shift_configs_x = new double[max_index];
            shift_configs_y = new double[max_index];
            
            for (int k = 0; k < max_index; ++k){
                shift_configs_x[k] = shift_configs[k][0];
                shift_configs_y[k] = shift_configs[k][1];
            }
            
        }
        
        
        // Construct and prepare the pairs arrays for broadcasting
        
        printf("Building inter and intra pairs. \n");
        std::vector<std::vector<int> > inter_pairs_vec;
        std::vector<std::vector<int> > inter_supercell_vecs;
        
        h.getInterPairs(inter_pairs_vec,inter_supercell_vecs,opts);
        //printf("inter pairs done \n");
        std::vector<int> intra_pairs_vec_i;
        std::vector<int> intra_pairs_vec_j;
        std::vector<double> intra_pairs_vec_t;
        std::vector<std::vector<int> > intra_supercell_vecs;
        
        h.getIntraPairs(intra_pairs_vec_i, intra_pairs_vec_j, intra_pairs_vec_t, intra_supercell_vecs, opts);
        
        //printf("Inter and intra pair construction complete. \n");
        max_inter_pairs = static_cast<int>(inter_pairs_vec.size());
        max_intra_pairs = static_cast<int>(intra_pairs_vec_i.size());
        
        printf("Heterostructure has %d orbitals. \n", max_index);
        printf("%d entries expected from intra. \n", max_intra_pairs);
        printf("%d entries expected from inter. \n", max_inter_pairs);
        
        if (solver_space == 1){
            if (fft_from_file == 0){
                // Generate FFT of interlayer coupling for Momentum space
                int n_x = opts.getInt("fft_n_x");
                int n_y = opts.getInt("fft_n_y");
                int L_x = opts.getInt("fft_L_x");
                int L_y = opts.getInt("fft_L_y");
                int length_x = opts.getInt("fft_length_x");
                int length_y = opts.getInt("fft_length_y");
                std::string fft_file = opts.getString("fft_file");
                //
                
                printf("Making *fft.dat file. \n");
                
                h.makeInterFFTFile(n_x, n_y, L_x, L_y, length_x, length_y, fft_file);
                
                // fft_data is accessed via fft_data[orbital_1][orbital_2][position][real/cpx]
                
                // !! also when putting data[o1][o2] into H, will need to multiply every interlayer term by sqrt(Area(RL_1)*Area(RL_2))
                // !! Using real-space area: This term is ~ 4*pi^2/sqrt(Area1*Area2) ~ 7.4308 for tBLG
                
            } else {
                printf("Warning: Not creating FFT file, assuming a *fft.dat file exists in current directory. \n");
            }
        }
        
        MPI::COMM_WORLD.Bcast(&max_inter_pairs, 1, MPI_INT, root);
        MPI::COMM_WORLD.Bcast(&max_intra_pairs, 1, MPI_INT, root);
        
        inter_pairs_i = new int[max_inter_pairs];
        inter_pairs_j = new int[max_inter_pairs];
        if (boundary_condition == 1){
            inter_supercell_vecs_i = new int[max_inter_pairs];
            inter_supercell_vecs_j = new int[max_inter_pairs];
            
        }
        
        for(int x = 0; x < max_inter_pairs; ++x){
            inter_pairs_i[x] = inter_pairs_vec[x][0];
            inter_pairs_j[x] = inter_pairs_vec[x][1];
            if (boundary_condition == 1){
                inter_supercell_vecs_i[x] = inter_supercell_vecs[x][0];
                inter_supercell_vecs_j[x] = inter_supercell_vecs[x][1];
            }
        }
        
        intra_pairs_i = new int[max_intra_pairs];
        intra_pairs_j = new int[max_intra_pairs];
        intra_pairs_t = new double[max_intra_pairs];
        if (boundary_condition == 1){
            intra_supercell_vecs_i = new int[max_intra_pairs];
            intra_supercell_vecs_j = new int[max_intra_pairs];
            
        }
        for(int x = 0; x < max_intra_pairs; ++x){
            
            intra_pairs_i[x] = intra_pairs_vec_i[x];
            intra_pairs_j[x] = intra_pairs_vec_j[x];
            intra_pairs_t[x] = intra_pairs_vec_t[x];
            if (boundary_condition == 1){
                intra_supercell_vecs_i[x] = intra_supercell_vecs[x][0];
                intra_supercell_vecs_j[x] = intra_supercell_vecs[x][1];
            }
        }
        
        // Get and prepare the index_to_pos array for broadcasting
        printf("Building indexToPos arrays. \n");
        index_to_pos_x = new double[max_index];
        index_to_pos_y = new double[max_index];
        index_to_pos_z = new double[max_index];
        
        h.getIndexToPos(index_to_pos_x,0);
        h.getIndexToPos(index_to_pos_y,1);
        h.getIndexToPos(index_to_pos_z,2);
        printf("indexToPos complete. \n");
        
        // POSITION DEBUG PRINT
        /*
         
         for (int k = 0; k < max_index; ++k)
         printf("%lf, %lf, %lf \n",index_to_pos_x[k],index_to_pos_y[k],index_to_pos_z[k]);
         */
        //
        
    }
    
    if (rank != root){
        
        // Allocate memory to receive pair and "index to grid" information
        MPI::COMM_WORLD.Bcast(&max_index, 1, MPI_INT, root);
        MPI::COMM_WORLD.Bcast(&max_inter_pairs, 1, MPI_INT, root);
        MPI::COMM_WORLD.Bcast(&max_intra_pairs, 1, MPI_INT, root);
        
        index_to_grid_i = new int[max_index];
        index_to_grid_j = new int[max_index];
        index_to_grid_l = new int[max_index];
        index_to_grid_s = new int[max_index];
        
        inter_pairs_i = new int[max_inter_pairs];
        inter_pairs_j = new int[max_inter_pairs];
        if (boundary_condition == 1){
            inter_supercell_vecs_i = new int[max_inter_pairs];
            inter_supercell_vecs_j = new int[max_inter_pairs];
            
        }
        
        intra_pairs_i = new int[max_intra_pairs];
        intra_pairs_j = new int[max_intra_pairs];
        intra_pairs_t = new double[max_intra_pairs];
        if (boundary_condition == 1){
            intra_supercell_vecs_i = new int[max_intra_pairs];
            intra_supercell_vecs_j = new int[max_intra_pairs];
            
        }
        
        index_to_pos_x = new double[max_index];
        index_to_pos_y = new double[max_index];
        index_to_pos_z = new double[max_index];
        
        if (strain_type == 2){
            shift_configs_x = new double[max_index];
            shift_configs_y = new double[max_index];
        }
        
    }
    
    MPI::COMM_WORLD.Bcast(index_to_grid_i, max_index, MPI_INT, root);
    MPI::COMM_WORLD.Bcast(index_to_grid_j, max_index, MPI_INT, root);
    MPI::COMM_WORLD.Bcast(index_to_grid_l, max_index, MPI_INT, root);
    MPI::COMM_WORLD.Bcast(index_to_grid_s, max_index, MPI_INT, root);
    
    MPI::COMM_WORLD.Bcast(inter_pairs_i, max_inter_pairs, MPI_INT, root);
    MPI::COMM_WORLD.Bcast(inter_pairs_j, max_inter_pairs, MPI_INT, root);
    
    if (boundary_condition == 1){
        MPI::COMM_WORLD.Bcast(inter_supercell_vecs_i, max_inter_pairs, MPI_INT, root);
        MPI::COMM_WORLD.Bcast(inter_supercell_vecs_j, max_inter_pairs, MPI_INT, root);
        
    }
    
    MPI::COMM_WORLD.Bcast(intra_pairs_i, max_intra_pairs, MPI_INT, root);
    MPI::COMM_WORLD.Bcast(intra_pairs_j, max_intra_pairs, MPI_INT, root);
    MPI::COMM_WORLD.Bcast(intra_pairs_t, max_intra_pairs, MPI_DOUBLE, root);
    if (boundary_condition == 1){
        MPI::COMM_WORLD.Bcast(intra_supercell_vecs_i, max_intra_pairs, MPI_INT, root);
        MPI::COMM_WORLD.Bcast(intra_supercell_vecs_j, max_intra_pairs, MPI_INT, root);
        
    }
    
    MPI::COMM_WORLD.Bcast(index_to_pos_x, max_index, MPI_DOUBLE, root);
    MPI::COMM_WORLD.Bcast(index_to_pos_y, max_index, MPI_DOUBLE, root);
    MPI::COMM_WORLD.Bcast(index_to_pos_z, max_index, MPI_DOUBLE, root);
    int* index_to_grid = new int[max_index*4];
    
    if (strain_type == 2){
        MPI::COMM_WORLD.Bcast(shift_configs_x, max_index, MPI_DOUBLE, root);
        MPI::COMM_WORLD.Bcast(shift_configs_y, max_index, MPI_DOUBLE, root);
    }
    
    for (int k = 0; k < max_index; ++k){
        index_to_grid[k*4 + 0] = index_to_grid_i[k];
        index_to_grid[k*4 + 1] = index_to_grid_j[k];
        index_to_grid[k*4 + 2] = index_to_grid_l[k];
        index_to_grid[k*4 + 3] = index_to_grid_s[k];
    }
    
    delete[] index_to_grid_i;
    delete[] index_to_grid_j;
    delete[] index_to_grid_l;
    delete[] index_to_grid_s;
    
    int* inter_pairs = new int[2*max_inter_pairs];
    
    std::vector< std::vector<int> > inter_sc_vecs;
    if (boundary_condition == 1){
        inter_sc_vecs.resize(max_inter_pairs);
        
    }
    for (int x = 0; x < max_inter_pairs; ++x){
        inter_pairs[x*2 + 0] = inter_pairs_i[x];
        inter_pairs[x*2 + 1] = inter_pairs_j[x];
        if (boundary_condition == 1){
            inter_sc_vecs[x].resize(2);
            inter_sc_vecs[x][0] = inter_supercell_vecs_i[x];
            inter_sc_vecs[x][1] = inter_supercell_vecs_j[x];
        }
    }
    
    delete inter_pairs_i;
    delete inter_pairs_j;
    if (boundary_condition == 1){
        delete inter_supercell_vecs_i;
        delete inter_supercell_vecs_j;
        
    }
    
    int* intra_pairs = new int[2*max_intra_pairs];
    
    
    std::vector< std::vector<int> > intra_sc_vecs;
    if (boundary_condition == 1){
        intra_sc_vecs.resize(max_intra_pairs);
        
    }
    
    for (int x = 0; x < max_intra_pairs; ++x){
        intra_pairs[x*2 + 0] = intra_pairs_i[x];
        intra_pairs[x*2 + 1] = intra_pairs_j[x];
        if (boundary_condition == 1){
            intra_sc_vecs[x].resize(2);
            intra_sc_vecs[x][0] = intra_supercell_vecs_i[x];
            intra_sc_vecs[x][1] = intra_supercell_vecs_j[x];
        }
    }
    
    delete intra_pairs_i;
    delete intra_pairs_j;
    if (boundary_condition == 1){
        delete intra_supercell_vecs_i;
        delete intra_supercell_vecs_j;
        
    }
    
    double* index_to_pos = new double[3*max_index];
    
    for (int k = 0; k < max_index; ++k){
        index_to_pos[k*3 + 0] = index_to_pos_x[k];
        index_to_pos[k*3 + 1] = index_to_pos_y[k];
        index_to_pos[k*3 + 2] = index_to_pos_z[k];
        
    }
    
    delete index_to_pos_x;
    delete index_to_pos_y;
    delete index_to_pos_z;
    
    
    if (strain_type == 2){
        shift_configs.resize(max_index);
        for (int k = 0; k < max_index; ++k){
            shift_configs[k].resize(2);
            shift_configs[k][0] = shift_configs_x[k];
            shift_configs[k][1] = shift_configs_y[k];
        }
        delete shift_configs_x;
        delete shift_configs_y;
    }
    
    time(&constructEnd);
    
    // Call the next method, which will send out jobs depending on the solver type to each worker and have them construct a matrix specific to that job.
    constructMatrix(index_to_grid,index_to_pos,inter_pairs,inter_sc_vecs,intra_pairs,intra_pairs_t,intra_sc_vecs,shift_configs,v_work,target_indices);
    
    delete index_to_grid;
    delete index_to_pos;
    delete inter_pairs;
    delete intra_pairs;
    delete intra_pairs_t;
    
}

void Locality::constructMatrix(int* index_to_grid,double* index_to_pos, int* inter_pairs, std::vector< std::vector<int> > inter_sc_vecs,
                               int* intra_pairs, double* intra_pairs_t, std::vector< std::vector<int> > intra_sc_vecs,
                               std::vector< std::vector<double> > shift_configs, std::vector< std::vector<int> > v_work,
                               std::vector< std::vector<int> > target_indices){
    
    time(&solveStart);
    
    int solver_type = opts.getInt("solver_type");
    //printf("Entering construct matrix \n");
    //if (rank == print_rank)
    //printf("rank %d entering constructMatrix(). \n", rank);
    
    // 1: Chebyshev polynomial sampling of eigenvalue spectrum (DOS)
    if(solver_type < 7){
        if (rank == root) {
            rootChebSolve(index_to_grid, index_to_pos,
                          inter_pairs, inter_sc_vecs,
                          intra_pairs, intra_pairs_t, intra_sc_vecs,
                          v_work, target_indices);
        } else {
            workerChebSolve(index_to_grid,index_to_pos,inter_pairs,inter_sc_vecs,intra_pairs,intra_pairs_t,intra_sc_vecs,shift_configs);
        }
    }
    
    time(&solveEnd);
}

void Locality::sendRootWork(Job_params jobIn, int target_r){
    
    jobIn.sendParams(target_r);
    
}

void Locality::recursiveShiftCalc(std::vector<Job_params>& jobArray, std::vector< std::vector<double> > shifts, int solver_type, int nShifts, int maxJobs, int num_sheets, int num_shift_sheets, std::vector<int> shift_sheets, std::vector< std::vector<int> > target_indices) {
    
    // Uniform sample over a grid
    if (solver_type == 1){
        
        if (num_shift_sheets == 0){
            
            int jobID = (int)jobArray.size() + 1;
            
            int n_targets = (int)target_indices[0].size();
            std::vector<int> targets;
            targets.resize(n_targets);
            for (int t = 0; t < n_targets; ++t){
                targets[t] = target_indices[0][t];
            }
            
            Job_params tempJob(opts);
            tempJob.setParam("shifts",shifts);
            tempJob.setParam("jobID",jobID);
            tempJob.setParam("max_jobs",maxJobs);
            tempJob.setParam("num_targets",n_targets);
            tempJob.setParam("target_list",targets);
            jobArray.push_back(tempJob);
            
        } else {
            
            int tar_sheet = shift_sheets[num_shift_sheets-1];
            
            for (int i = 0; i < nShifts; ++i){
                for (int j = 0; j < nShifts; ++j){
                    
                    double x = (1.0/(double) (nShifts))*i;
                    double y = (1.0/(double) (nShifts))*j;
                    
                    shifts[tar_sheet][0] = x;
                    shifts[tar_sheet][1] = y;
                    shifts[tar_sheet][2] = 0;
                    
                    recursiveShiftCalc(jobArray, shifts, solver_type, nShifts, maxJobs, num_sheets, num_shift_sheets-1, shift_sheets, target_indices);
                    
                    
                }
            }
        }
        
    }
    
}

void Locality::rootChebSolve(int* index_to_grid, double* index_to_pos,
                             int* inter_pairs, std::vector< std::vector<int> > inter_sc_vecs,
                             int* intra_pairs, double* intra_pairs_t, std::vector< std::vector<int> > intra_sc_vecs,
                             std::vector< std::vector<int> > v_work, std::vector< std::vector<int> > target_indices) {
    
    int solver_type = opts.getInt("solver_type");
    int observable_type = opts.getInt("observable_type");
    int solver_space = opts.getInt("solver_space");
    int diagonalize = opts.getInt("diagonalize");
    int d_weights = opts.getInt("d_weights");
    int d_vecs = opts.getInt("d_vecs");
    int d_cond = opts.getInt("d_cond");
    int k_sampling = opts.getInt("k_sampling");
    int kpm_trace = opts.getInt("kpm_trace");
    
    int mlmc = opts.getInt("mlmc");
    
    int maxJobs;
    int nShifts;
    
    int currentJob = 0;
    int length;
    
    std::vector<Job_params> jobArray;
    Mlmc_handler mlmc_h;
    
    
    if (mlmc == 1){
        mlmc_h.setup(opts);
    }
    
    int num_sheets = (int)sdata.size();
    
    // Make job arrays
    
    if (solver_space == 0) {
        
        // Uniform sample over a grid
        if (solver_type == 1) {
            nShifts = opts.getInt("nShifts");
            int num_shift_sheets = opts.getInt("num_shift_sheets");
            std::vector<int> shift_sheets = opts.getIntVec("shift_sheets");
            
            maxJobs = pow(nShifts*nShifts,num_shift_sheets);
            
            std::vector< std::vector<double> > shifts;
            shifts.resize(num_sheets);
            for (int s = 0; s < num_sheets; ++s){
                shifts[s].resize(3);
                shifts[s][0] = 0;
                shifts[s][1] = 0;
                shifts[s][2] = 0;
            }
            
            recursiveShiftCalc(jobArray,shifts, solver_type, nShifts, maxJobs, num_sheets, num_shift_sheets, shift_sheets, target_indices);
            
        }
        
        // Linecuts through the unit cell
        if (solver_type == 2){
            
            nShifts = opts.getInt("nShifts");
            maxJobs = nShifts;
            
            int num_lc_points = opts.getInt("num_lc_points");
            std::vector< std::vector<double> > lc_points = opts.getDoubleMat("lc_points");
            
            for (int i = 0; i < maxJobs; ++i){
                
                double x = (1.0/((double) maxJobs))*i;
                
                double total_x = x*(num_lc_points-1);
                
                int lc_index = (int)floor(total_x);
                double eff_x = fmod(total_x,1);
                
                std::vector< std::vector<double> > shifts;
                shifts.resize(num_sheets);
                
                for(int s = 0; s < num_sheets; ++s){
                    shifts[s].resize(3);
                    shifts[s][0] = 0;
                    shifts[s][1] = 0;
                    shifts[s][2] = 0;
                }
                
                shifts[num_sheets-1][0] = lc_points[lc_index][0] + eff_x*(lc_points[lc_index+1][0] - lc_points[lc_index][0]);
                shifts[num_sheets-1][1] = lc_points[lc_index][1] + eff_x*(lc_points[lc_index+1][1] - lc_points[lc_index][1]);
                shifts[num_sheets-1][2] = 0;
                
                // Custom shifts for making "perfect" 2H aligned bilayers in TMDCs
                // Assumes the top layer has 180 rotation, and bottom layer has 0 rotation
                //shifts[1][0] = 2.0/3.0;
                //shifts[1][1] = 1.0/3.0;
                
                printf("shifts[%d] = [%lf %lf] \n",num_sheets-1,shifts[num_sheets-1][0],shifts[num_sheets-1][1]);
                
                // for debugging: check the linecutting algorithm
                //printf("job %d has shift = [%lf, %lf] (lc_index = %d, eff_x = %lf) \n",i,shifts[num_sheets-1][0],shifts[num_sheets-1][1],lc_index,eff_x);
                
                int n_targets = (int)target_indices[0].size();
                std::vector<int> targets;
                targets.resize(n_targets);
                
                for (int t = 0; t < n_targets; ++t){
                    targets[t] = target_indices[0][t];
                }
                
                Job_params tempJob(opts);
                tempJob.setParam("shifts",shifts);
                tempJob.setParam("jobID",i+1);
                tempJob.setParam("max_jobs",maxJobs);
                tempJob.setParam("num_targets",n_targets);
                tempJob.setParam("target_list",targets);
                jobArray.push_back(tempJob);
                
            }
        }
        
        // Vacancy solvers
        
        // With MLMC
        
        if (solver_type == 3){
            std::vector< std::vector<int> > mlmc_vacs;
            std::vector< std::vector<int> > mlmc_ids;
            
            getVacanciesFromFile(mlmc_vacs, mlmc_ids, opts);
            
            int maxJobs = (int)mlmc_vacs.size();
            
            //printf("maxJobs = %d \n",maxJobs);
            
            for (int i = 0; i < maxJobs; ++i){
                
                // no shifts
                std::vector< std::vector<double> > shifts;
                shifts.resize(num_sheets);
                
                for(int s = 0; s < num_sheets ; ++s){
                    shifts[s].resize(3);
                    shifts[s][0] = 0;
                    shifts[s][1] = 0;
                    shifts[s][2] = 0;
                }
                
                int n_vac = (int)mlmc_vacs[i].size();
                
                std::vector<int> vacancies;
                vacancies.resize(n_vac);
                for (int v = 0; v < n_vac; ++v){
                    vacancies[v] = mlmc_vacs[i][v];
                }
                
                // catch the zero vacancy case!
                // vacancy file has -1 as empty list element
                // but we subtract 1 during fileread (thus -2)
                if (n_vac == 1 && vacancies[0] == -2){
                    n_vac = 0;
                }
                
                int n_targets = (int)mlmc_ids[i].size();
                std::vector<int> targets;
                targets.resize(n_targets);
                
                for (int t = 0; t < n_targets; ++t){
                    targets[t] = mlmc_ids[i][t];
                }
                
                Job_params tempJob(opts);
                tempJob.setParam("shifts",shifts);
                //tempJob.setParam("target_list",targets,n_targets);
                tempJob.setParam("num_vacancies",n_vac);
                tempJob.setParam("vacancy_list",vacancies);
                tempJob.setParam("target_list",targets);
                tempJob.setParam("jobID",mlmc_ids[i][0]);
                tempJob.setParam("mlmcID",mlmc_ids[i][0]);
                tempJob.setParam("mlmc_clusterID",mlmc_ids[i][1]);
                tempJob.setParam("max_jobs",maxJobs);
                jobArray.push_back(tempJob);
                
            }
            
        }
        
        // Without MLMC
        
        if (solver_type == 4){
            
            maxJobs = (int)v_work.size();
            
            for (int i = 0; i < maxJobs; ++i){
                
                std::vector< std::vector<double> > shifts;
                shifts.resize(num_sheets);
                for(int s = 0; s < num_sheets ; ++s){
                    shifts.resize(3);
                    shifts[s][0] = 0;
                    shifts[s][1] = 0;
                    shifts[s][2] = 0;
                }
                
                int n_vac = (int)v_work[i].size();
                std::vector<int> vacancies;
                vacancies.resize(n_vac);
                for (int v = 0; v < n_vac; ++v){
                    vacancies[v] = v_work[i][v];
                }
                
                int n_targets = (int)target_indices[i].size();
                std::vector<int> targets;
                targets.resize(n_targets);
                for (int t = 0; t < n_targets; ++t){
                    targets[t] = target_indices[i][t];
                }
                
                
                Job_params tempJob(opts);
                tempJob.setParam("shifts",shifts);
                tempJob.setParam("target_list",targets);
                tempJob.setParam("vacancy_list",vacancies);
                tempJob.setParam("jobID",i+1);
                tempJob.setParam("max_jobs",maxJobs);
                jobArray.push_back(tempJob);
            }
        }
        
        // Strain solver(s)
        
        if (solver_type == 5){
            
            int strain_type = opts.getInt("strain_type");
            
            if (strain_type == 3) { // realspace strain type
                
                nShifts = opts.getInt("nShifts");
                maxJobs = nShifts;
                
                for (int i = 0; i < maxJobs; ++i){
                    
                    double max_shift = 20;
                    double x = 2*max_shift*((i/(double)maxJobs) - 0.5);
                    
                    std::vector< std::vector<double> > shifts;
                    shifts.resize(num_sheets);
                    
                    for(int s = 0; s < num_sheets; ++s){
                        shifts[s].resize(3);
                        shifts[s][0] = 0;
                        shifts[s][1] = 0;
                        shifts[s][2] = 0;
                    }
                    
                    int n_targets = (int)target_indices[0].size();
                    std::vector<int> targets;
                    targets.resize(n_targets);
                    
                    for (int t = 0; t < n_targets; ++t){
                        targets[t] = target_indices[0][t];
                    }
                    
                    
                    
                    Job_params tempJob(opts);
                    tempJob.setParam("strain_shift",x);
                    tempJob.setParam("shifts",shifts);
                    tempJob.setParam("jobID",i+1);
                    tempJob.setParam("max_jobs",maxJobs);
                    tempJob.setParam("num_targets",n_targets);
                    tempJob.setParam("target_list",targets);
                    jobArray.push_back(tempJob);
                }
                
            } else {
                
                maxJobs = (int)target_indices.size();
                
                for (int i = 0; i < maxJobs; ++i){
                    
                    std::vector< std::vector<double> > shifts;
                    shifts.resize(num_sheets);
                    for(int s = 0; s < num_sheets ; ++s){
                        shifts[s].resize(3);
                        shifts[s][0] = 0;
                        shifts[s][1] = 0;
                        shifts[s][2] = 0;
                    }
                    
                    int n_targets = (int)target_indices[i].size();
                    std::vector<int> targets;
                    targets.resize(n_targets);
                    for (int t = 0; t < n_targets; ++ t){
                        targets[t] = target_indices[i][t];
                    }
                    
                    Job_params tempJob(opts);
                    tempJob.setParam("shifts",shifts);
                    tempJob.setParam("target_list",targets);
                    tempJob.setParam("jobID",i+1);
                    tempJob.setParam("max_jobs",maxJobs);
                    tempJob.setParam("num_targets",n_targets);
                    tempJob.setParam("target_list",targets);
                    jobArray.push_back(tempJob);
                }
            }
        }
        
        // runs one job with a fixed shift
        
        if (solver_type == 6){
            maxJobs = 1;
            std::vector< std::vector<double> > shifts;
            shifts.resize(num_sheets);
            for(int s = 0; s < num_sheets ; ++s){
                shifts[s].resize(3);
                if (s == 0){
                    shifts[s][0] = 0;
                    shifts[s][1] = 0;
                    shifts[s][2] = 0;
                } else {
                    shifts[s][0] = 1.0/3.0;
                    shifts[s][1] = 1.0/3.0;
                    shifts[s][2] = 0;
                }
            }
            
            int n_targets = (int)target_indices[0].size();
            std::vector<int> targets;
            targets.resize(n_targets);
            for (int t = 0; t < n_targets; ++ t){
                targets[t] = target_indices[0][t];
            }
            
            Job_params tempJob(opts);
            tempJob.setParam("shifts",shifts);
            tempJob.setParam("target_list",targets);
            tempJob.setParam("jobID",1);
            tempJob.setParam("max_jobs",maxJobs);
            tempJob.setParam("num_targets",n_targets);
            tempJob.setParam("target_list",targets);
            jobArray.push_back(tempJob);
            
        }
        
        if (k_sampling == 1){
            
            int k_type = opts.getInt("k_type");
            
            if (k_type == 0){ // grid sampling
                int num_k1 = opts.getInt("num_k1");
                int num_k2 = opts.getInt("num_k2");
                std::vector< std::vector<double> > a = opts.getDoubleMat("supercell");
                std::vector< std::vector<double> > b = getReciprocal(a);
                
                int k_jobID = 1;
                std::vector<Job_params> k_jobArray;
                std::ofstream outFile_k;
                const char* ext = "_kgrid.dat";
                outFile_k.open ((job_name + ext).c_str());
                
                for (int i = 0; i < (int)jobArray.size(); ++i){
                    for (int k1 = 0; k1 < num_k1; ++k1){
                        for (int k2 = 0; k2 < num_k2; ++k2){
                            
                            Job_params tempJob(jobArray[i]);
                            //tempJob.setParam("origJobID",tempJob.getInt("jobID"));
                            tempJob.setParam("jobID",k_jobID);
                            
                            std::vector<double> k_vec;
                            k_vec.resize(2);
                            k_vec[0] = (k1/(double)num_k1)*b[0][0] + (k2/(double)num_k2)*b[1][0];
                            k_vec[1] = (k1/(double)num_k1)*b[0][1] + (k2/(double)num_k2)*b[1][1];
                            std::vector<double> k_vec_grid;
                            k_vec_grid.resize(2);
                            k_vec_grid[0] = (k1/(double)num_k1);
                            k_vec_grid[1] = (k2/(double)num_k2);
                            tempJob.setParam("k_vec",k_vec);
                            tempJob.setParam("k_vec_grid",k_vec_grid);
                            tempJob.setParam("unit_cell",a);
                            
                            outFile_k << k_vec[0] << "," << k_vec[1] << "\n";
                            
                            k_jobArray.push_back(tempJob);
                            ++k_jobID;
                        }
                    }
                }
                jobArray = k_jobArray;
            
                
                // Sample around twisted supercell BZ
            } else if (k_type == 1){
                
                double num_k1 = opts.getInt("num_k1");
                maxJobs = num_k1;
                
                int type = opts.getInt("supercell_type");
                
                int base_s1 = 0;
                int base_s2 = 1;
                if (type == 1){
                    std::vector<int> sc_groups = opts.getIntVec("sc_groups");
                    int first_type = sc_groups[0];
                    for (int idx = 1; idx < sc_groups.size(); ++idx){
                        if (sc_groups[idx] != first_type){
                            base_s2 = idx;
                            break;
                        }
                    }
                }
                
                std::vector< std::vector<double> > b1 = getReciprocal(base_s1);
                std::vector< std::vector<double> > b2 = getReciprocal(base_s2); // reciprocal lattice of the monolayers
                
                
                printf("b1 = [%lf %lf; %lf %lf] \n",b1[0][0],b1[0][1],b1[1][0],b1[1][1]);
                printf("b2 = [%lf %lf; %lf %lf] \n",b2[0][0],b2[0][1],b2[1][0],b2[1][1]);
                
                //printf("shift = [%lf, %lf] \n",shifts[0],shifts[1]);
                
                double k_1[2];
                double k_2[2];
                
                // k_1 = K of layer 1, k_2 = K of layer 2
                
                k_1[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/2)*b1[1][0] + sin(M_PI/2)*b1[1][1]);
                k_1[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/2)*b1[1][0] + cos(M_PI/2)*b1[1][1]);
                
                k_2[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/2)*b2[1][0] + sin(M_PI/2)*b2[1][1]);
                k_2[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/2)*b2[1][0] + cos(M_PI/2)*b2[1][1]);
                
                double k[2];
                double gamma[2];
                double m[2];
                
                k[0] = k_1[0];
                k[1] = k_1[1];
                
                double d = (1.0/2.0)*sqrt((k_2[0] - k_1[0])*(k_2[0] - k_1[0]) + (k_2[1] - k_1[1])*(k_2[1] - k_1[1]));
                double x_dir[2];
                double y_dir[2];
                
                y_dir[0] = (k_2[0] - k_1[0])/(2.0*d);
                y_dir[1] = (k_2[1] - k_1[1])/(2.0*d);
                x_dir[0] = cos(-M_PI/2)*y_dir[0] - sin(-M_PI/2)*y_dir[1];
                x_dir[1] = sin(-M_PI/2)*y_dir[0] + cos(-M_PI/2)*y_dir[0];
                
                gamma[0] = k[0] + d*y_dir[0] + sqrt(3)*d*x_dir[0];
                gamma[1] = k[1] + d*y_dir[1] + sqrt(3)*d*x_dir[1];
                
                m[0] = k[0] + d*y_dir[0];
                m[1] = k[1] + d*y_dir[1];
                
                printf("k = [%lf, %lf], gamma = [%lf, %lf], m = [%lf, %lf] \n",k[0],k[1],gamma[0],gamma[1],m[0],m[1]);
                
                // save k line data
                int k_jobID = 1;
                std::vector<Job_params> k_jobArray;
                std::ofstream outFile_k;
                const char* ext = "_kline.dat";
                outFile_k.open ((job_name + ext).c_str());
                
                for (int j = 0; j < (int) jobArray.size(); ++j){
                    for (int i = 0; i < maxJobs; ++i){
                        //double x = (1.0/((double) maxJobs))*i;
                        //double x = 0.5 + 2.0*(1.0/((double) maxJobs))*(i - maxJobs/2.0);
                        //double x = .3333 + (1.0/((double) maxJobs))*(i-maxJobs/2)/(20);
                        double c = (3.0/((double) maxJobs))*i;
                        
                        double k_x = 0;
                        double k_y = 0;
                        
                        if (c <= 1) {
                            k_x = (1.0-c)*k[0] + (c-0.0)*gamma[0];
                            k_y = (1.0-c)*k[1] + (c-0.0)*gamma[1];
                        } else if (c <= 2) {
                            k_x = (2.0-c)*gamma[0] + (c-1.0)*m[0];
                            k_y = (2.0-c)*gamma[1] + (c-1.0)*m[1];
                        } else {
                            k_x = (3.0-c)*m[0] + (c-2.0)*k[0];
                            k_y = (3.0-c)*m[1] + (c-2.0)*k[1];
                        }
                        
                        Job_params tempJob(jobArray[j]);
                        //tempJob.setParam("origJobID",tempJob.getInt("jobID"));
                        tempJob.setParam("jobID",k_jobID);
                        
                        std::vector<double> k_vec;
                        k_vec.resize(2);
                        k_vec[0] = k_x;
                        k_vec[1] = k_y;
                        printf("k = [%lf, %lf]\n",k_vec[0],k_vec[1]);
                        tempJob.setParam("k_vec",k_vec);
                        
                        outFile_k << k_vec[0] << "," << k_vec[1] << "\n";
                        
                        k_jobArray.push_back(tempJob);
                        ++k_jobID;
                        
                    }
                }
                
                jobArray = k_jobArray;
                
                // Sample around Monolayer BZ
            }else if (k_type == 2){
                
                double num_k1 = opts.getInt("num_k1");
                maxJobs = num_k1;
                
                std::vector< std::vector<double> > b1 = getReciprocal(0);
                
                printf("b1 = [%lf %lf; %lf %lf] \n",b1[0][0],b1[0][1],b1[1][0],b1[1][1]);
                
                //printf("shift = [%lf, %lf] \n",shifts[0],shifts[1]);
                
                double gamma[2];
                double m[2];
                double k[2];
                
                double dot = b1[0][0]*b1[1][0] + b1[0][1]*b1[1][1];    //dot product
                double det = b1[0][0]*b1[1][1] - b1[1][0]*b1[0][1];    //determinant
                double phi = atan2(det, dot);                                           // angle between the two
                
                gamma[0] = 0.0;
                gamma[1] = 0.0;
                m[0] = b1[0][0]/2.0;
                m[1] = b1[0][1]/2.0;
                k[0] = (m[0]*cos(phi/2.0) - m[1]*sin(phi/2.0))/cos(phi/2.0);
                k[1] = (m[0]*sin(phi/2.0) + m[1]*cos(phi/2.0))/cos(phi/2.0);
                /*
                 double k[2];
                 double k_prime[2];
                 
                 k[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/2)*b1[1][0] + sin(M_PI/2)*b1[1][1]);
                 k[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/2)*b1[1][0] + cos(M_PI/2)*b1[1][1]);
                 
                 k_prime[0] = cos(M_PI/3)*k[0] - sin(M_PI/3)*k[1];
                 k_prime[1] = sin(M_PI/3)*k[0] + cos(M_PI/3)*k[1];
                 
                 double gamma[2];
                 double m[2];
                 
                 double d = (1.0/2.0)*sqrt((k_prime[0] - k[0])*(k_prime[0] - k[0]) + (k_prime[1] - k[1])*(k_prime[1] - k[1]));
                 double x_dir[2];
                 double y_dir[2];
                 
                 y_dir[0] = (k_prime[0] - k[0])/(2.0*d);
                 y_dir[1] = (k_prime[1] - k[1])/(2.0*d);
                 x_dir[0] = cos(-M_PI/2)*y_dir[0] - sin(-M_PI/2)*y_dir[1];
                 x_dir[1] = sin(-M_PI/2)*y_dir[0] + cos(-M_PI/2)*y_dir[0];
                 
                 gamma[0] = k[0] + d*y_dir[0] + sqrt(3)*d*x_dir[0];
                 gamma[1] = k[1] + d*y_dir[1] + sqrt(3)*d*x_dir[1];
                 
                 m[0] = k[0] + d*y_dir[0];
                 m[1] = k[1] + d*y_dir[1];
                 */
                
                printf("k = [%lf, %lf], gamma = [%lf, %lf], m = [%lf, %lf] \n",k[0],k[1],gamma[0],gamma[1],m[0],m[1]);
                
                int k_jobID = 1;
                std::vector<Job_params> k_jobArray;
                for (int j = 0; j < (int) jobArray.size(); ++j){
                    for (int i = 0; i < maxJobs; ++i){
                        //double x = (1.0/((double) maxJobs))*i;
                        //double x = 0.5 + 2.0*(1.0/((double) maxJobs))*(i - maxJobs/2.0);
                        //double x = .3333 + (1.0/((double) maxJobs))*(i-maxJobs/2)/(20);
                        double c = (3.0/((double) maxJobs))*i;
                        
                        double k_x = 0;
                        double k_y = 0;
                        
                        /*
                         if (c <= 1) {
                         k_x = (1.0-c)*k[0] + (c-0.0)*gamma[0];
                         k_y = (1.0-c)*k[1] + (c-0.0)*gamma[1];
                         } else if (c <= 2) {
                         k_x = (2.0-c)*gamma[0] + (c-1.0)*m[0];
                         k_y = (2.0-c)*gamma[1] + (c-1.0)*m[1];
                         } else {
                         k_x = (3.0-c)*m[0] + (c-2.0)*k[0];
                         k_y = (3.0-c)*m[1] + (c-2.0)*k[1];
                         }
                         */
                        
                        if (c <= 1) {
                            k_x = (1.0-c)*gamma[0] + (c-0.0)*m[0];
                            k_y = (1.0-c)*gamma[1] + (c-0.0)*m[1];
                        } else if (c <= 2) {
                            k_x = (2.0-c)*m[0] + (c-1.0)*k[0];
                            k_y = (2.0-c)*m[1] + (c-1.0)*k[1];
                        } else {
                            k_x = (3.0-c)*k[0] + (c-2.0)*gamma[0];
                            k_y = (3.0-c)*k[1] + (c-2.0)*gamma[1];
                        }
                        
                        Job_params tempJob(jobArray[j]);
                        //tempJob.setParam("origJobID",tempJob.getInt("jobID"));
                        tempJob.setParam("jobID",k_jobID);
                        
                        std::vector<double> k_vec;
                        k_vec.resize(2);
                        k_vec[0] = k_x;
                        k_vec[1] = k_y;
                        printf("k = [%lf, %lf]\n",k_vec[0],k_vec[1]);
                        tempJob.setParam("k_vec",k_vec);
                        
                        k_jobArray.push_back(tempJob);
                        ++k_jobID;
                        
                    }
                }
                
                jobArray = k_jobArray;
                // zoe: sampling around the trilayer supercell BZ (commensurated only)
            }else if (k_type == 3){
                
                int num_k1 = opts.getInt("num_k1");
                maxJobs = 3*num_k1;
                
                int type = opts.getInt("supercell_type");   // type 1: bilayer; type 3: trilayer
                
                if ( (type != 1) && (type != 3)){
                    printf("You need to have a supercell for this type of the k sampling!");
                }
                
                
                int M = opts.getInt("m_supercell");
                int N = opts.getInt("n_supercell");
                
                double theta = acos((N*N + 4*N*M + M*M)/(2.0*(N*N + N*M + M*M)));
                
                std::vector<int> sc_groups = opts.getIntVec("sc_groups");
                
                int A11;
                int A12;
                int A21;
                int A22;
                
                if (sc_groups.size() == 2){
                    
                    if(M>N){
                        
                        A11 = N;
                        A12 = M;
                        A21 = -M;
                        A22 = (M+N);
                        
                    } else {
                        
                        A11 = M;
                        A12 = N;
                        A21 = -N;
                        A22 = (M+N);
                    }
                    
                    
                } else if (sc_groups.size() == 3) {
                    
                    A11 = 0;
                    A12 = pow(M,2) + (M*N) + pow(N,2);
                    A21 = -(pow(N,2) + (M*N) + pow(M,2));
                    A22 = pow(N,2)+ M*N + pow(M,2);
                    
                }
                
                // Unit cell vector
                std::vector< std::vector<double> > unitCell = sdata[0].a;
                
                std::vector< std::vector<double> > sc_here;
                sc_here.resize(2);
                sc_here[0].resize(2);
                sc_here[1].resize(2);
                
                sc_here[0][0] = A11*unitCell[0][0] + A12*unitCell[1][0];
                sc_here[0][1] = A11*unitCell[0][1] + A12*unitCell[1][1];
                sc_here[1][0] = A21*unitCell[0][0] + A22*unitCell[1][0];
                sc_here[1][1] = A21*unitCell[0][1] + A22*unitCell[1][1];
                
                std::vector< std::vector<double> > b_sc = getReciprocal(sc_here);
                
                // define high symmetry points
                double k[2];
                double m[2];
                double gamma[2];
                
                double r3 = sqrt(3);
                double alpha = M_PI/6;
                std::vector< std::vector<double> > rot30;
                rot30.resize(2);
                rot30[0].resize(2);
                rot30[0].resize(2);
                
                rot30[0] = {cos(alpha), sin(alpha)};
                rot30[1] = {-sin(alpha), cos(alpha)};
                
                
                for(int i = 0; i < 2; i++){
                    m[i] = 0.5*(b_sc[0][i] + b_sc[1][i]);
                }
                
                
                for(int i = 0; i < 2; i++){
                    k[i] = 0;
                    for(int j = 0; j < 2; j++){
                        k[i] = k[i] + rot30[i][j] * b_sc[1][j] / r3;
                    }
                }
                
                gamma[0] = 0;
                gamma[1] = 0;
                
                printf("\n");
                printf("b_sc = [%lf %lf; %lf %lf] \n",b_sc[0][0],b_sc[0][1],b_sc[1][0],b_sc[1][1]);
                std::cout << "K = [" << k[0] << "," << k[1] << "]\n";
                std::cout << "M = [" << m[0] << "," << m[1] << "]\n";
                std::cout << "Gamma = [" << gamma[0] << "," << gamma[1] << "]\n";
                
                
                double s1 = (k[1] - gamma[1]) / (k[0] - gamma[0]);
                double s2 = (m[1] - k[1]) / (m[0] - k[0]);
                double s3 = (gamma[1] - m[1]) / (gamma[0] - m[0]);
                
                double b1 = (k[0]*gamma[1] - gamma[0]*k[1]) / (k[0] - gamma[0]);
                double b2 = (m[0]*k[1] - k[0]*m[1]) / (m[0] - k[0]);
                double b3 = (gamma[0]*m[1] - m[0]*gamma[1]) / (gamma[0] - m[0]);
                
                double l1 = k[0]-gamma[0];
                double l2 = m[0]-k[0] ;
                double l3 = gamma[0] - m[0];
                
                // create and save k line data
                int k_jobID = 1;
                std::vector<Job_params> k_jobArray;
                
                std::ofstream outFile_k;
                const char* ext = "_kline.dat";
                outFile_k.open ((job_name + ext).c_str());
                
                double kx[maxJobs];
                double ky[maxJobs];
                
                for (int j = 0; j < (int) jobArray.size(); ++j){
                    for (int i = 0; i < maxJobs; ++i){
                        
                        if (i < num_k1){
                            // Line 1
                            kx[i] = gamma[0] + i * l1 / num_k1;
                            ky[i] = gamma[1] + s1 * i * l1 / num_k1 + b1;
                            
                        } else if (i < num_k1*2){
                            // Line 2
                            int n = i - num_k1;
                            kx[i] = k[0] + n * l2 / num_k1;
                            ky[i] = k[1] + s2 * n * l2 / num_k1;
                        } else {
                            // Line 3
                            int n = i - 2*num_k1;
                            if (m[0] ==m [1]){
                                double len = gamma[1] - m[1];
                                kx[i] = m[0];
                                ky[i] = m[1] +  n * len / num_k1;
                            } else{
                                kx[i] = m[0] + n * l3 / num_k1;
                                ky[i] = m[1] + s3 * n * l3 / num_k1 + b3;
                            }
                        }
                        
                        
                        Job_params tempJob(jobArray[j]);
                        
                        tempJob.setParam("jobID",k_jobID);
                        
                        std::vector<double> k_vec;
                        k_vec.resize(2);
                        k_vec[0] = kx[i];
                        k_vec[1] = ky[i];
                        printf("k = [%lf, %lf]\n",k_vec[0],k_vec[1]);
                        tempJob.setParam("k_vec",k_vec);
                        
                        outFile_k << k_vec[0] << "," << k_vec[1] << "\n";
                        
                        k_jobArray.push_back(tempJob);
                        ++k_jobID;
                        
                    }
                }
                
                jobArray = k_jobArray;
                
            } else if (k_type == 4){ // Connecting the three Dirac points of the trilayer
                int num_k1 = opts.getInt("num_k1");
                maxJobs = 2*num_k1;
                
                int type = opts.getInt("supercell_type");   // type 3: trilayer
                
                if ( (type != 3)){
                    printf("You don't have a trilayer supercell!");
                }
                
                int M = opts.getInt("m_supercell");
                int N = opts.getInt("n_supercell");
                
                double theta = acos((N*N + 4*N*M + M*M)/(2.0*(N*N + N*M + M*M)));
                
                std::vector< std::vector<double> > b1 = getReciprocal(0);
                std::vector< std::vector<double> > b2 = getReciprocal(1);
                std::vector< std::vector<double> > b3 = getReciprocal(2); // reciprocal lattice of the monolayers
                
                printf("b1 = [%lf %lf; %lf %lf] \n",b1[0][0],b1[0][1],b1[1][0],b1[1][1]);
                printf("b2 = [%lf %lf; %lf %lf] \n",b2[0][0],b2[0][1],b2[1][0],b2[1][1]);
                printf("b3 = [%lf %lf; %lf %lf] \n",b3[0][0],b3[0][1],b3[1][0],b3[1][1]);
                
                // K points of each layer
                double k_1[2];
                double k_2[2];
                double k_3[2];
                
                for(int i = 0; i < 2; i++){
                    k_1[i] = 1.0/3.0*(2.0*b1[0][i]+b1[1][i]);
                    k_2[i] = 1.0/3.0*(2.0*b2[0][i]+b2[1][i]);
                    k_3[i] = 1.0/3.0*(2.0*b3[0][i]+b3[1][i]);
                }

                printf("K_L1 = [%lf %lf] \n", k_1[0], k_1[1]);
                printf("K_L2 = [%lf %lf] \n", k_2[0], k_2[1]);
                printf("K_L3 = [%lf %lf] \n", k_3[0], k_3[1]);
                
                // create and save k line data
                int k_jobID = 1;
                std::vector<Job_params> k_jobArray;
                
                std::ofstream outFile_k;
                const char* ext = "_kline.dat";
                outFile_k.open ((job_name + ext).c_str());
                
                double s1 = (k_2[1]-k_1[1])/(k_2[0]-k_1[0]);
                double s2 = (k_3[1]-k_2[1])/(k_3[0]-k_2[0]);
                double m1 = k_2[1] - s1*k_2[0];
                double m2 = k_3[1] - s2*k_3[0];
                double l1 = k_2[0]-k_1[0];
                double l2 = k_3[0]-k_2[0];
                printf("[%lf %lf] \n", s2, m2);
                
                double kx[maxJobs+1];
                double ky[maxJobs+1]; // add one because we need to include K_L3
                
                for (int j = 0; j < (int) jobArray.size(); ++j){
                    for (int i = 0; i < maxJobs+1; ++i){
                        
                        if (i < num_k1){
                            // Line 1, K_L1 to K_L2
                            kx[i] = k_1[0] + i * l1 / num_k1;
                            ky[i] = s1 * kx[i] + m1;
                            
                        } else if (i < maxJobs){
                            // Line 2, K_L2 to K_L3
                            int n = i - num_k1;
                            kx[i] = k_2[0] + n * l2 / num_k1;
                            ky[i] = s2 * kxLq:Quanta:Quanta:Quanta[i] + m2;
                            
                        } else{
                            kx[i] = k_3[0];
                            ky[i] = k_3[1];
                        }
                    
                        Job_params tempJob(jobArray[j]);
                        
                        tempJob.setParam("jobID",k_jobID);
                        
                        std::vector<double> k_vec;
                        k_vec.resize(2);
                        k_vec[0] = kx[i];
                        k_vec[1] = ky[i];
                        printf("k = [%lf, %lf]\n",k_vec[0],k_vec[1]);
                        tempJob.setParam("k_vec",k_vec);
                        
                        outFile_k << k_vec[0] << "," << k_vec[1] << "\n";
                        
                        k_jobArray.push_back(tempJob);
                        ++k_jobID;
                    }
                }
                jobArray = k_jobArray;
            }
            
        }
        
        // random sampling of the trace
        if (kpm_trace == 1){
            
            int num_trace_samps = opts.getInt("kpm_trace_samps");
            if (num_trace_samps == -1) {
                
                // should get SMALLER as matrix size gets bigger (need less samples)
                num_trace_samps = (int) 1.0e6 / ( (double) max_index ); // just a guess
                
                // some edge cases!
                if (num_trace_samps < 1)
                    num_trace_samps = 1;
                
                // if matrix is too small, don't need too many samples
                if (num_trace_samps > 100)
                    num_trace_samps = 100;
                
            }
            
            
            // setup random number generator
            std::default_random_engine generator;
            std::normal_distribution<double> distribution(0.0,1.0); // mean = 0, std = 1
            
            // loop over existing jobs
            int trace_jobID = 1;
            std::vector<Job_params> trace_jobArray;
            for (int i = 0; i < (int)jobArray.size(); ++i){
                
                int local_max_index;
                // take care of vacancy case
                if (solver_type == 3){
                    local_max_index = max_index - jobArray[i].getInt("num_vacancies");
                } else {
                    local_max_index = max_index;
                }
                
                Job_params tempJob(jobArray[i]);
                //tempJob.setParam("origJobID",tempJob.getInt("jobID"));
                tempJob.setParam("jobID",trace_jobID);
                
                std::vector< std::vector<double> > trace_vectors;
                trace_vectors.resize(num_trace_samps);
                for (int s = 0; s < num_trace_samps; ++s){
                    
                    trace_vectors[s].resize(local_max_index);
                    
                    for (int rand_idx = 0; rand_idx < local_max_index; ++rand_idx){
                        trace_vectors[s][rand_idx] = distribution(generator); // random gaussian distribution
                    }
                    
                }
                
                tempJob.setParam("kpm_trace_vectors",trace_vectors);
                tempJob.setParam("num_targets",num_trace_samps);
                trace_jobArray.push_back(tempJob);
                ++trace_jobID;
                
            }
            jobArray = trace_jobArray;
            
        }
        
    } else if (solver_space == 1) {
        
        int mom_vf_only = opts.getInt("mom_vf_only");
        int num_mom_groups = opts.getInt("num_mom_groups");
        if (num_mom_groups != 2){
            printf("!!LEO2D ERROR!! Momentum-space only supported for NUM_MOM_GROUPS = 2 \n");
            throw std::runtime_error("Locality::rootChevSolve Momentum-space only supported for NUM_MOM_GROUPS = 2!!");
        }
        std::vector< std::vector<int> > mom_groups = opts.getIntMat("mom_groups");
        
        printf("comparing sheet %d to sheet %d \n",mom_groups[0][0]+1,mom_groups[1][0]+1);
        
        std::vector< std::vector<double> > b1 = getReciprocal(mom_groups[0][0]);
        std::vector< std::vector<double> > b2 = getReciprocal(mom_groups[1][0]);
        
        // jobs that try to compute fermi-velocity only
        if (mom_vf_only == 1){
            
            // hard-coded control of slope method (type = 0) and gamme-point gap method (type = 1)
            //int mom_vf_type = 0;
            
            int mom_vf_type = opts.getInt("mom_vf_type");
            printf("mom_vf_type = %d \n",mom_vf_type);
            
            if (mom_vf_type == 0){
                maxJobs = 2;
            }else{
                maxJobs = 1;
            }
            
            double k_1[2];
            double k_2[2];
            
            k_1[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/2)*b1[1][0] + sin(M_PI/2)*b1[1][1]);
            k_1[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/2)*b1[1][0] + cos(M_PI/2)*b1[1][1]);
            
            k_2[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/2)*b2[1][0] + sin(M_PI/2)*b2[1][1]);
            k_2[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/2)*b2[1][0] + cos(M_PI/2)*b2[1][1]);
            
            printf("k_1 = [%lf, %lf], k_2 = [%lf, %lf] \n",k_1[0], k_1[1], k_2[0], k_2[1]);
            
            double k[2];
            double gamma[2];
            double m[2];
            
            double gamma_2[2];
            
            k[0] = k_1[0];
            k[1] = k_1[1];
            
            // d is the distance between k_1 and k_2
            double d = (1.0/2.0)*sqrt((k_2[0] - k_1[0])*(k_2[0] - k_1[0]) + (k_2[1] - k_1[1])*(k_2[1] - k_1[1]));
            double x_dir[2];
            double y_dir[2];
            
            y_dir[0] = (k_2[0] - k_1[0])/(2.0*d);
            y_dir[1] = (k_2[1] - k_1[1])/(2.0*d);
            x_dir[0] = cos(-M_PI/2)*y_dir[0] - sin(-M_PI/2)*y_dir[1];
            x_dir[1] = sin(-M_PI/2)*y_dir[0] + cos(-M_PI/2)*y_dir[0];
            
            gamma[0] = k[0] + d*y_dir[0] + sqrt(3)*d*x_dir[0];
            gamma[1] = k[1] + d*y_dir[1] + sqrt(3)*d*x_dir[1];
            
            gamma_2[0] = k[0] + d*y_dir[0] - sqrt(3)*d*x_dir[0];
            gamma_2[1] = k[1] + d*y_dir[1] - sqrt(3)*d*x_dir[1];
            
            m[0] = k[0] + d*y_dir[0];
            m[1] = k[1] + d*y_dir[1];
            
            //printf("k = [%lf, %lf], gamma = [%lf, %lf], m = [%lf, %lf] \n",k[0],k[1],gamma[0],gamma[1],m[0],m[1]);
            
            for (int i = 0; i < maxJobs; ++i){
                
                
                //double c = (3.0/((double) maxJobs))*i;
                
                //double c = (6.0/((double) maxJobs))*i;
                double c = 0.0;
                if(mom_vf_type == 0){
                    c = 0.0;
                } else{
                    c = 1.0;
                }
                
                if (i == 1){
                    // only go a little along K-Gamma line to compute the fermi-velocity
                    c = 0.06;
                }
                
                double shift_x = 0;
                double shift_y = 0;
                
                if (c <= 1) {
                    shift_x = (1.0-c)*k[0] + (c-0.0)*gamma[0];
                    shift_y = (1.0-c)*k[1] + (c-0.0)*gamma[1];
                } else if (c <= 2) {
                    shift_x = (2.0-c)*gamma[0] + (c-1.0)*m[0];
                    shift_y = (2.0-c)*gamma[1] + (c-1.0)*m[1];
                } else if (c <= 3) {
                    shift_x = (3.0-c)*m[0] + (c-2.0)*k[0];
                    shift_y = (3.0-c)*m[1] + (c-2.0)*k[1];
                } else if (c <= 4) {
                    shift_x = (4.0-c)*k[0] + (c-3.0)*gamma_2[0];
                    shift_y = (4.0-c)*k[1] + (c-3.0)*gamma_2[1];
                } else if (c <= 5) {
                    shift_x = (5.0-c)*gamma_2[0] + (c-4.0)*m[0];
                    shift_y = (5.0-c)*gamma_2[1] + (c-4.0)*m[1];
                } else {
                    shift_x = (6.0-c)*m[0] + (c-5.0)*k[0];
                    shift_y = (6.0-c)*m[1] + (c-5.0)*k[1];
                }
                
                std::vector< std::vector<double> > shifts;
                shifts.resize(num_sheets);
                
                
                for(int s = 0; s < num_sheets; ++s){
                    shifts[s].resize(3);
                    shifts[s][0] = shift_x;
                    shifts[s][1] = shift_y;
                    shifts[s][2] = 0;
                }
                
                int n_targets = (int)target_indices[0].size();
                std::vector<int> targets;
                targets.resize(n_targets);
                for (int t = 0; t < n_targets; ++t){
                    targets[t] = target_indices[0][t];
                }
                
                Job_params tempJob(opts);
                tempJob.setParam("shifts",shifts);
                tempJob.setParam("jobID",i+1);
                tempJob.setParam("max_jobs",maxJobs);
                tempJob.setParam("num_targets",n_targets);
                tempJob.setParam("target_list",targets);
                
                jobArray.push_back(tempJob);
                
                
            }
        }
        
        // Square sampling
        else if (solver_type == 1){
            nShifts = opts.getInt("nShifts");
            maxJobs = nShifts*nShifts;
            
            // Changes to make us cut through K and K' points (i.e. Gamma -> K -> K' -> Gamma)
            //b1[1][0] = -b1[1][0];
            //b2[1][0] = -b2[1][0];
            
            //printf("b1 = [%lf %lf; %lf %lf] \n",b1[0][0],b1[0][1],b1[1][0],b1[1][1]);
            //printf("b2 = [%lf %lf; %lf %lf] \n",b2[0][0],b2[0][1],b2[1][0],b2[1][1]);
            
            //printf("shift = [%lf, %lf] \n",shifts[0],shifts[1]);
            
            double d_vec1[2];
            double d_vec2[2];
            
            // Here we assume b1 ~ b2, the small angle limit of similar (or even identical) input lattices
            
            d_vec1[0] = b2[0][0] - b1[0][0];
            d_vec1[1] = b2[0][1] - b1[0][1];
            
            d_vec2[0] = b2[1][0] - b1[1][0];
            d_vec2[1] = b2[1][1] - b1[1][1];
            
            //printf("d_vec1 = [%lf, %lf], d_vec2 = [%lf, %lf] \n",d_vec1[0],d_vec1[1],d_vec2[0],d_vec2[1]);
            
            // p1 and p2 are the K points of layer 1 and 2, respectively
            
            double k_1[2];
            double k_2[2];
            
            
            k_1[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/2)*b1[1][0] + sin(M_PI/2)*b1[1][1]);
            k_1[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/2)*b1[1][0] + cos(M_PI/2)*b1[1][1]);
            
            k_2[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/2)*b2[1][0] + sin(M_PI/2)*b2[1][1]);
            k_2[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/2)*b2[1][0] + cos(M_PI/2)*b2[1][1]);
            
            /*
             p1[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/6)*b1[1][0] + sin(M_PI/6)*b1[1][1]);
             p1[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/6)*b1[1][0] + cos(M_PI/6)*b1[1][1]);
             
             p2[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/6)*b2[1][0] + sin(M_PI/6)*b2[1][1]);
             p2[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/6)*b2[1][0] + cos(M_PI/6)*b2[1][1]);
             */
            
            printf("k1 = [%lf, %lf], k2 = [%lf, %lf] \n",k_1[0],k_1[1],k_2[0],k_2[1]);
            printf("d_vec1 = [%lf, %lf], d_vec2 =[%lf, %lf] \n",d_vec1[0],d_vec1[1],d_vec2[0],d_vec2[1]);
            for (int i = 0; i < nShifts; ++i){
                double x = (1.0/((double) nShifts))*i;
                for (int j = 0; j < nShifts; ++j){
                    double y = (1.0/((double) nShifts))*j;
                    
                    std::vector< std::vector<double> > shifts;
                    shifts.resize(num_sheets);
                    for(int s = 0; s < num_sheets; ++s){
                        shifts[s].resize(3);
                        shifts[s][0] = x*d_vec1[0] + y*d_vec2[0] + k_1[0];
                        shifts[s][1] = x*d_vec1[1] + y*d_vec2[1] + k_1[1];
                        shifts[s][2] = 0;
                    }
                    //printf("k = [%lf, %lf] \n",shifts[0][0],shifts[0][1]);
                    
                    int n_targets = (int)target_indices[0].size();
                    std::vector<int> targets;
                    targets.resize(n_targets);
                    for (int t = 0; t < n_targets; ++t){
                        targets[t] = target_indices[0][t];
                    }
                    
                    Job_params tempJob(opts);
                    tempJob.setParam("shifts",shifts);
                    tempJob.setParam("jobID",i*nShifts + j + 1);
                    tempJob.setParam("max_jobs",maxJobs);
                    tempJob.setParam("num_targets",n_targets);
                    tempJob.setParam("target_list",targets);
                    jobArray.push_back(tempJob);
                }
                
            }
            
        }
        
        // Linecut sampling
        else if (solver_type == 2){
            nShifts = opts.getInt("nShifts");
            maxJobs = nShifts;
            
            // cut through K and K' points (i.e. Gamma -> K -> K' -> Gamma)
            
            //printf("b1 = [%lf %lf; %lf %lf] \n",b1[0][0],b1[0][1],b1[1][0],b1[1][1]);
            //printf("b2 = [%lf %lf; %lf %lf] \n",b2[0][0],b2[0][1],b2[1][0],b2[1][1]);
            
            //printf("shift = [%lf, %lf] \n",shifts[0],shifts[1]);
            
            double k_1[2];
            double k_2[2];
            
            
            // k_1 = K of layer 1, k_2 = K of layer 2
            /*
             k_1[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/6)*b1[1][0] + sin(M_PI/6)*b1[1][1]);
             k_1[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/6)*b1[1][0] + cos(M_PI/6)*b1[1][1]);
             
             k_2[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/6)*b2[1][0] + sin(M_PI/6)*b2[1][1]);
             k_2[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/6)*b2[1][0] + cos(M_PI/6)*b2[1][1]);
             */
            
            k_1[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/2)*b1[1][0] + sin(M_PI/2)*b1[1][1]);
            k_1[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/2)*b1[1][0] + cos(M_PI/2)*b1[1][1]);
            
            k_2[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/2)*b2[1][0] + sin(M_PI/2)*b2[1][1]);
            k_2[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/2)*b2[1][0] + cos(M_PI/2)*b2[1][1]);
            
            double k[2];
            double gamma[2];
            double m[2];
            
            double gamma_2[2];
            
            k[0] = k_1[0];
            k[1] = k_1[1];
            
            // d is the distance between k_1 and k_2
            double d = (1.0/2.0)*sqrt((k_2[0] - k_1[0])*(k_2[0] - k_1[0]) + (k_2[1] - k_1[1])*(k_2[1] - k_1[1]));
            double x_dir[2];
            double y_dir[2];
            
            y_dir[0] = (k_2[0] - k_1[0])/(2.0*d);
            y_dir[1] = (k_2[1] - k_1[1])/(2.0*d);
            x_dir[0] = cos(-M_PI/2)*y_dir[0] - sin(-M_PI/2)*y_dir[1];
            x_dir[1] = sin(-M_PI/2)*y_dir[0] + cos(-M_PI/2)*y_dir[0];
            
            gamma[0] = k[0] + d*y_dir[0] + sqrt(3)*d*x_dir[0];
            gamma[1] = k[1] + d*y_dir[1] + sqrt(3)*d*x_dir[1];
            
            gamma_2[0] = k[0] + d*y_dir[0] - sqrt(3)*d*x_dir[0];
            gamma_2[1] = k[1] + d*y_dir[1] - sqrt(3)*d*x_dir[1];
            
            m[0] = k[0] + d*y_dir[0];
            m[1] = k[1] + d*y_dir[1];
            
            printf("k = [%lf, %lf], gamma = [%lf, %lf], m = [%lf, %lf], gamma_2 = [%lf, %lf] \n",k[0],k[1],gamma[0],gamma[1],m[0],m[1],gamma_2[0],gamma_2[1]);
            
            for (int i = 0; i < maxJobs; ++i){
                //double x = (1.0/((double) maxJobs))*i;
                //double x = 0.5 + 2.0*(1.0/((double) maxJobs))*(i - maxJobs/2.0);
                //double x = .3333 + (1.0/((double) maxJobs))*(i-maxJobs/2)/(20);
                
                
                //double c = (3.0/((double) maxJobs))*i;
                
                double c = (6.0/((double) maxJobs))*i;
                
                double shift_x = 0;
                double shift_y = 0;
                
                /*
                 if (c <= 1) {
                 shift_x = (1.0-c)*k[0] + (c-0.0)*gamma[0];
                 shift_y = (1.0-c)*k[1] + (c-0.0)*gamma[1];
                 } else if (c <= 2) {
                 shift_x = (2.0-c)*gamma[0] + (c-1.0)*m[0];
                 shift_y = (2.0-c)*gamma[1] + (c-1.0)*m[1];
                 } else {
                 shift_x = (3.0-c)*m[0] + (c-2.0)*k[0];
                 shift_y = (3.0-c)*m[1] + (c-2.0)*k[1];
                 }
                 */
                
                if (c <= 1) {
                    shift_x = (1.0-c)*k[0] + (c-0.0)*gamma[0];
                    shift_y = (1.0-c)*k[1] + (c-0.0)*gamma[1];
                } else if (c <= 2) {
                    shift_x = (2.0-c)*gamma[0] + (c-1.0)*m[0];
                    shift_y = (2.0-c)*gamma[1] + (c-1.0)*m[1];
                } else if (c <= 3) {
                    shift_x = (3.0-c)*m[0] + (c-2.0)*k[0];
                    shift_y = (3.0-c)*m[1] + (c-2.0)*k[1];
                } else if (c <= 4) {
                    shift_x = (4.0-c)*k[0] + (c-3.0)*gamma_2[0];
                    shift_y = (4.0-c)*k[1] + (c-3.0)*gamma_2[1];
                } else if (c <= 5) {
                    shift_x = (5.0-c)*gamma_2[0] + (c-4.0)*m[0];
                    shift_y = (5.0-c)*gamma_2[1] + (c-4.0)*m[1];
                } else {
                    shift_x = (6.0-c)*m[0] + (c-5.0)*k[0];
                    shift_y = (6.0-c)*m[1] + (c-5.0)*k[1];
                }
                
                std::vector< std::vector<double> > shifts;
                shifts.resize(num_sheets);
                
                //printf("k[%d] = [%lf, %lf] \n",i,shift_x, shift_y);
                
                
                
                for(int s = 0; s < num_sheets; ++s){
                    shifts[s].resize(3);
                    shifts[s][0] = shift_x;
                    shifts[s][1] = shift_y;
                    //shifts[s][0] = x*b1[0][0] + (1-x)*b1[1][0];
                    //shifts[s][1] = x*b1[0][1] + (1-x)*b1[1][1];
                    shifts[s][2] = 0;
                }
                
                int n_targets = (int)target_indices[0].size();
                std::vector<int> targets;
                targets.resize(n_targets);
                for (int t = 0; t < n_targets; ++t){
                    targets[t] = target_indices[0][t];
                }
                
                Job_params tempJob(opts);
                tempJob.setParam("shifts",shifts);
                tempJob.setParam("jobID",i+1);
                tempJob.setParam("max_jobs",maxJobs);
                tempJob.setParam("num_targets",n_targets);
                tempJob.setParam("target_list",targets);
                
                jobArray.push_back(tempJob);
                
                
                
            }
        }
        
        else if (solver_type == 3 || solver_type == 4 || solver_type == 5){
            printf("!!WARNING!!: Momentum-space mode (solver_space = M) is NOT compatible with vacancy/strain solvers. \n");
        }
        
    }
    
    // Now do the job send/recv loop
    
    MPI::Status status;
    
    //std::vector<std::vector<double> > result_array;
    //result_array.resize(maxJobs);
    
    maxJobs = (int) jobArray.size();
    std::vector<Job_params> result_array;
    result_array.resize(maxJobs);
    
    // try to give each worker its first job
    for (int r = 1; r < size; ++r) {
        if (currentJob < maxJobs) {
            sendRootWork(jobArray[currentJob], r);
            ++currentJob;                    // one more job sent out!
            printf("rank %d has sent job to worker rank %d (%d/%d) \n", rank, r, currentJob, maxJobs);
            
        }
        
    }
    
    
    // printf("rank %d has sent first batch of work... \n", rank);
    
    // Receive results and dispense new work
    
    while (currentJob < maxJobs) {
        
        int source;
        Job_params temp_result;
        
        source = temp_result.recvSpool();
        
        int jobTag = temp_result.getInt("jobID");
        
        if (mlmc == 0){
            result_array[jobTag-1] = (temp_result);
        } else if (mlmc == 1){
            mlmc_h.process(temp_result);
        }
        
        sendRootWork(jobArray[currentJob], source);
        ++currentJob;                        // one more job sent out!
        printf("rank %d has sent job to worker rank %d (%d/%d) \n", rank, source, currentJob, maxJobs);
    }
    
    // Receive final work
    // and tell workers to exit workerMatrixSolve()
    
    int r_done = 0;
    
    while (r_done < size-1){
        
        int source;
        Job_params temp_result;
        
        source = temp_result.recvSpool();
        
        int jobTag = temp_result.getInt("jobID");
        
        if (mlmc == 0){
            result_array[jobTag-1] = temp_result;
        } else if (mlmc == 1){
            mlmc_h.process(temp_result);
        }
        
        printf("rank %d has received final work from rank %d. \n", rank, source);
        
        Job_params tempJob;
        
        // solver_type == -1 is the STOPTAG for workers
        tempJob.setParam("solver_type", -1);
        
        sendRootWork(tempJob, source);
        ++r_done;
        
    }
    
    if (mlmc == 0){
        
        std::cout << "Saving " << job_name << ".cheb to disk. \n";
        std::ofstream outFile;
        const char* extension =".cheb";
        outFile.open( (job_name + extension).c_str() );
        
        for (int job = 0; job < maxJobs; ++job){
            Param_tools::save(result_array[job], outFile);
        }
        outFile.close();
        
        int wan_save = opts.getInt("wan_save");
        if (wan_save == 1){
            
            std::ofstream valsOutFile;
            const char* extension_vals = ".vals";
            valsOutFile.open( (job_name + extension_vals).c_str() );
            
            std::ofstream vecsOutFile;
            const char* extension_vecs = ".vecs";
            vecsOutFile.open( (job_name + extension_vecs).c_str() );
            
            Param_tools::wannier_raw_eig_save(result_array, valsOutFile, vecsOutFile);
            valsOutFile.close();
            vecsOutFile.close();
            
            std::ofstream winOutFile;
            const char* extension_win =".win";
            winOutFile.open( (job_name + extension_win).c_str() );
            Param_tools::wannier_win_save(result_array, winOutFile);
            winOutFile.close();
            
            std::ofstream eigOutFile;
            const char* extension_eig =".eig";
            eigOutFile.open( (job_name + extension_eig).c_str() );
            Param_tools::wannier_eig_save(result_array, eigOutFile);
            eigOutFile.close();
            
            std::ofstream mmnOutFile;
            const char* extension_mmn =".mmn";
            mmnOutFile.open( (job_name + extension_mmn).c_str() );
            Param_tools::wannier_mmn_save(result_array, mmnOutFile);
            mmnOutFile.close();
            
        }
        
        printTiming(result_array);
        
    } else if (mlmc == 1){
        mlmc_h.save();
    }
    
    
}

void Locality::workerChebSolve(int* index_to_grid, double* index_to_pos,
                               int* inter_pairs, std::vector< std::vector<int> > inter_sc_vecs,
                               int* intra_pairs, double* intra_pairs_t, std::vector< std::vector<int> > intra_sc_vecs,
                               std::vector< std::vector<double> > shift_configs) {
    
    MPI::Status status;
    // ---------------------
    // Enter MPI worker loop
    // ---------------------
    
    // Load FFT data once if its momentum-space code...
    if (opts.getInt("solver_space") == 1){
        // Load FFTW data into Momentum_coupling object
        
        // The following 7 variables should eventually be taken as input parameters in Loc_params.cpp from hstruct.in
        // They define the settings used to generate the interlayer_fft.dat file
        
        int n_x = opts.getInt("fft_n_x");
        int n_y = opts.getInt("fft_n_y");
        int L_x = opts.getInt("fft_L_x");
        int L_y = opts.getInt("fft_L_y");
        int length_x = opts.getInt("fft_length_x");
        int length_y = opts.getInt("fft_length_y");
        std::string fft_file = opts.getString("fft_file");
        //
        
        fftw_inter.fft_setup(L_x,L_y,length_x,length_y,fft_file);
    }
    
    // For testing how LAPACK returns eigenvectors...
    /*
     double vals[9];
     vals[0] = 1.0;
     vals[1] = 6.0;
     vals[2] = 3.0;
     vals[3] = 6.0;
     vals[4] = 2.0;
     vals[5] = 1.0;
     vals[6] = 3.0;
     vals[7] = 1.0;
     vals[8] = 3.0;
     DMatrix test_mat;
     test_mat.setup(3,3,vals);
     test_mat.debugPrint();
     
     std::vector<double> vls;
     vls.resize(3);
     DMatrix vcs;
     
     test_mat.eigenSolve(vls, vcs, 'V', 'A', 0, 3);
     vcs.debugPrint();
     */
    
    while (1) {
        
        Job_params jobIn;
        jobIn.recvParams(root);
        int solver_type = jobIn.getInt("solver_type");
        
        // If worker gets solver_type = -1 it ends this loop
        if (solver_type == -1) {
            //printf("rank %d received STOPTAG. \n", rank);
            return;
        }
        
        Job_params results_out(jobIn);
        
        int jobID = jobIn.getInt("jobID");
        int max_jobs = jobIn.getInt("max_jobs");
        int magOn = jobIn.getInt("magOn");
        int observable_type = jobIn.getInt("observable_type");
        int solver_space = jobIn.getInt("solver_space");
        int diagonalize = jobIn.getInt("diagonalize");
        int boundary_condition = jobIn.getInt("boundary_condition");
        int chiral_on = jobIn.getInt("chiral_on");
        int d_kpm_dos = jobIn.getInt("d_kpm_dos");
        int d_vecs =  jobIn.getInt("d_vecs");
        int d_cond =  jobIn.getInt("d_cond");
        int d_weights = jobIn.getInt("d_weights");
        int k_sampling = jobIn.getInt("k_sampling");
        int ballistic_transport = jobIn.getInt("ballistic_transport");
        int kpm_trace = jobIn.getInt("kpm_trace");
        
        int complex_matrix = 0;
        
        if (magOn == 1 || solver_space == 1 || k_sampling == 1 || chiral_on == 1 || ballistic_transport == 1){
            complex_matrix = 1;
        }
        
        results_out.setParam("complex_matrix",complex_matrix);
        
        // If not STOPTAG, start timing the solver
        time_t tempStart;
        time(&tempStart);
        solverTimes.push_back(tempStart);
        
        clock_t timeStart;
        timeStart = clock();
        //results_out.saveTiming(timeStart, "START");
        
        double vacancy_chance = jobIn.getDouble("vacancy_chance");
        int num_targets = jobIn.getInt("num_targets");
        std::vector<int> target_list;
        
        if (num_targets > 0){
            target_list = jobIn.getIntVec("target_list");
        }
        
        int poly_order = jobIn.getInt("poly_order");
        
        std::vector< std::vector<double> > shifts = jobIn.getDoubleMat("shifts");
        //printf("rank %d received a job (%d/%d) \n", rank, jobID,max_jobs);
        
        // ------------------------------------------------------
        // Determine the work-specific positions of every orbital
        // ------------------------------------------------------
        
        // i2pos = "index to position"
        double i2pos[max_index*3];
        std::vector< std::vector<double> > s_configs;
        std::vector< std::vector< std::vector<double> > > strain;
        
        
        setConfigPositions(i2pos, index_to_pos, index_to_grid, s_configs, shift_configs, strain, jobIn);
        
        clock_t timePos;
        timePos = clock();
        Param_tools::saveTiming(results_out, ( ((double)timePos) - ((double)timeStart) )/CLOCKS_PER_SEC, "POS_UPDATE");
        
        // Keep track of variables for vacancy simulation
        
        int local_max_index = max_index;
        
        std::vector<int> current_index_reduction(max_index+1,0);    // tells us how to relabel our indices in the reduced size matrix
        int vacancy_counter = 0;                                    // total number of vacancies
        
        // Vacancy creation loop
        
        int num_vacancies = jobIn.getInt("num_vacancies");
        std::vector<int> vac_list;
        if (num_vacancies > 0){
            vac_list = jobIn.getIntVec("vacancy_list");
        }
        
        if (solver_type == 3 || solver_type == 4){
            for (int i = 0; i < max_index; ++i){
                
                current_index_reduction[i] = vacancy_counter;
                
                for (int j = 0; j < num_vacancies; ++j){
                    if(i == vac_list[j]){
                        vacancy_counter++;
                        break;
                        
                    }
                }
                
            }
            
            current_index_reduction[max_index] = vacancy_counter;
            local_max_index = max_index - vacancy_counter;
            
            /*
             printf("First few vacancies = [");
             for (int j = 0; j < 6; ++j){
             if (num_vacancies > j){
             printf("%d, ",vac_list[j]);
             }
             }
             if (num_vacancies > 5){
             printf("...] \n");
             } else{
             printf("] \n");
             }
             */
        }
        
        
        // -----------------------
        // Build the Sparse matrix
        // -----------------------
        SpMatrix H;
        SpMatrix dxH;
        SpMatrix dyH;
        
        SpMatrix J_B_x;
        SpMatrix J_B_y;
        SpMatrix J_T_x;
        SpMatrix J_T_y;
        SpMatrix J_TB_x;
        SpMatrix J_TB_y;
        
        double* alpha_0_x_arr = new double[num_targets*local_max_index];
        double* alpha_0_y_arr = new double[num_targets*local_max_index];
        
        if (solver_space == 0){
            if (complex_matrix == 0){
                generateRealH(H, dxH, dyH, alpha_0_x_arr, alpha_0_y_arr, jobIn, index_to_grid, i2pos,
                              inter_pairs, inter_sc_vecs, intra_pairs, intra_pairs_t,
                              intra_sc_vecs, strain, current_index_reduction, local_max_index);
            } else if (complex_matrix == 1) {
                generateCpxH(H, dxH, dyH,
                             J_B_x, J_B_y,J_T_x, J_T_y, J_TB_x, J_TB_y,
                             alpha_0_x_arr, alpha_0_y_arr, jobIn, index_to_grid, i2pos,
                             inter_pairs, inter_sc_vecs, intra_pairs, intra_pairs_t,
                             intra_sc_vecs, strain, current_index_reduction, local_max_index);
            }
        } else if (solver_space == 1){
            generateMomH(H, jobIn, index_to_grid, i2pos, inter_pairs, intra_pairs, intra_pairs_t, current_index_reduction, local_max_index);
        }
        
        clock_t timeMat;
        timeMat = clock();
        Param_tools::saveTiming(results_out, ( ((double)timeMat) - ((double)timePos) )/CLOCKS_PER_SEC, "MATRIX_GEN");
        
        // ---------------------
        // Chebyshev Computation
        // ---------------------
        
        if (ballistic_transport == 1){
            //Ballistic::runBallisticTransport(results_out, index_to_grid,i2pos, H);
            
            // !!!!
            // Have to hack in a file save here as sending large data over MPI crashes on local machine...
            std::cout << "Saving " << results_out.getString("job_name") << ".cheb to disk. \n";
            std::ofstream outFile;
            const char* extension =".cheb";
            outFile.open( (results_out.getString("job_name") + extension).c_str() );
            Param_tools::save(results_out, outFile);
            outFile.close();
            // end hackiness
            // !!!!
            
        } else if (diagonalize == 0){
            if (observable_type == 0){
                
                std::vector< std::vector<double> > cheb_coeffs;
                cheb_coeffs.resize(num_targets);
                
                
                for (int t = 0; t < num_targets; ++t){
                    cheb_coeffs[t].resize(poly_order);
                }
                
                if (kpm_trace == 0) {
                    computeDosKPM(cheb_coeffs, H, jobIn, current_index_reduction, local_max_index);
                } else {
                    computeDosTraceKPM(cheb_coeffs, H, jobIn, current_index_reduction, local_max_index);
                }
                
                results_out.setParam("cheb_coeffs",cheb_coeffs);
                
                if (kpm_trace == 0){
                    if (opts.getInt("dos_transform") == 1){
                        Param_tools::densityTransform(results_out);
                    }
                }
                
                // need a different transformation technique for kpm trace (average over moments before idct!)
                if (kpm_trace == 1) {
                    int num_vac = jobIn.getInt("num_vacancies");
                    double num_dofs = local_max_index - num_vac;
                    
                    Param_tools::densityTransformTrace(results_out, 1.0/num_dofs);
                }
                
            } else if (observable_type == 1){
                
                std::vector< std::vector<double> > cheb_coeffs;
                cheb_coeffs.resize(num_targets);
                
                for (int t = 0; t < num_targets; ++t){
                    cheb_coeffs[t].resize(poly_order*poly_order);
                }
                
                computeCondKPM(cheb_coeffs, H, dxH, jobIn, current_index_reduction, local_max_index, alpha_0_x_arr);
                
                results_out.setParam("kpm_M_xx",cheb_coeffs);
                if (opts.getInt("cond_transform") == 1){
                    Param_tools::conductivityTransform(results_out);
                }
                
            }
            
        } else if (diagonalize == 1) {
            int d_type = opts.getInt("d_type");
            
            if (complex_matrix == 0){
                
                std::vector<double> eigenvalue_array;
                if (d_type == 2){
                    eigenvalue_array.resize(8);
                } else {
                    eigenvalue_array.resize(local_max_index);
                }
                
                DMatrix eigvecs;
                
                DMatrix kpm_dos;
                DMatrix M_xx;
                DMatrix M_yy;
                DMatrix M_xy;
                
                computeEigen(eigenvalue_array, eigvecs, kpm_dos, M_xx, M_yy, M_xy, H, dxH, dyH, jobIn, current_index_reduction, local_max_index);
                
                std::vector< std::vector<double> > eigenweights;
                std::vector< std::vector<double> > eigenvectors;
                std::vector< std::vector<double> > cond_xx;
                std::vector< std::vector<double> > cond_yy;
                std::vector< std::vector<double> > cond_xy;
                
                results_out.setParam("eigenvalues",eigenvalue_array);
                
                double* eigenvector_array;
                eigenvector_array = eigvecs.getValPtr();
                
                if (d_kpm_dos == 1){
                    
                    std::vector<double> kpm_dos_out;
                    std::vector< std::vector<double> > kpm_dos_in;
                    
                    double* val_kpm_dos;
                    val_kpm_dos = kpm_dos.getValPtr();
                    
                    
                    for (int i = 0; i < poly_order; ++i){
                        std::vector<double> temp_vec;
                        
                        for (int j = 0; j < poly_order; ++j) {
                            temp_vec.push_back(val_kpm_dos[i + j*poly_order]);
                        }
                        
                        kpm_dos_in.push_back(temp_vec);
                    }
                    
                    Job_params dummy_result(results_out);
                    dummy_result.setParam("kpm_dos_mat",kpm_dos_in);
                    Param_tools::matrixResponseTransform(dummy_result, "kpm_dos_mat");
                    kpm_dos_in = dummy_result.getDoubleMat("kpm_dos_mat");
                    
                    for (int i = 0; i < poly_order; ++i){
                        kpm_dos_out.push_back(kpm_dos_in[i][i]);
                    }
                    
                    results_out.setParam("kpm_dos",kpm_dos_out);
                    
                }
                
                if (d_weights == 1){
                    
                    for (int t = 0; t < num_targets; ++t){
                        
                        int tar_here = target_list[t];
                        std::vector<double> temp_vec;
                        
                        for (int i = 0; i < local_max_index; ++i){
                            double w = eigenvector_array[i*local_max_index + tar_here];
                            temp_vec.push_back(w*w);
                        }
                        
                        eigenweights.push_back(temp_vec);
                    }
                    
                    results_out.setParam("eigenweights",eigenweights);
                    
                }
                
                if (d_vecs == 1){
                    
                    for (int j = 0; j < local_max_index; ++j){
                        
                        std::vector<double> temp_vec;
                        
                        for (int i = 0; i < local_max_index; ++i){
                            double w = eigenvector_array[i + j*local_max_index];
                            temp_vec.push_back(w);
                        }
                        
                        eigenvectors.push_back(temp_vec);
                    }
                    
                    results_out.setParam("eigenvectors",eigenvectors);
                    
                }
                
                if (d_cond > 0){
                    
                    double* val_M_xx;
                    val_M_xx = M_xx.getValPtr();
                    
                    double bz_area = 1.0;
                    if (boundary_condition == 1){
                        // need to normalize via BZ area...
                        std::vector< std::vector<double> > supercell = opts.getDoubleMat("supercell");
                        bz_area = Param_tools::computeReciprocalArea(supercell);
                    }
                    
                    for (int i = 0; i < poly_order; ++i){
                        std::vector<double> temp_xx;
                        
                        for (int j = 0; j < poly_order; ++j) {
                            temp_xx.push_back(val_M_xx[i + j*poly_order]);
                        }
                        
                        cond_xx.push_back(temp_xx);
                    }
                    
                    results_out.setParam("M_xx",cond_xx);
                    
                    if (d_cond > 1){
                        
                        double* val_M_yy;
                        val_M_yy = M_yy.getValPtr();
                        
                        double* val_M_xy;
                        val_M_xy = M_xy.getValPtr();
                        
                        for (int i = 0; i < poly_order; ++i){
                            std::vector<double> temp_yy;
                            std::vector<double> temp_xy;
                            
                            for (int j = 0; j < poly_order; ++j) {
                                temp_yy.push_back(val_M_yy[i + j*poly_order]);
                                temp_xy.push_back(val_M_xy[i + j*poly_order]);
                            }
                            
                            cond_yy.push_back(temp_yy);
                            cond_xy.push_back(temp_xy);
                            
                            results_out.setParam("M_yy", cond_yy);
                            results_out.setParam("M_xy", cond_xy);
                            
                        }
                    }
                    
                    // convert to E basis
                    Param_tools::conductivityTransform(results_out);
                    
                }
                
            } else if (complex_matrix == 1){
                
                
                int band_start_idx = 0;
                int band_end_idx = local_max_index;
                int wan_num_bands = opts.getInt("wan_num_bands");
                
                std::vector<double> eigenvalue_array;
                if (d_type == 2){
                    band_start_idx = 0;
                    band_end_idx = 8;
                    eigenvalue_array.resize(8);
                } else if (wan_num_bands > 0){
                    // select wan_num_bands around the center band value
                    int mid_idx = local_max_index/2.0;
                    int width = wan_num_bands/2.0;
                    band_start_idx = mid_idx - (width - 1);
                    band_end_idx = mid_idx + (width + 1);
                    eigenvalue_array.resize(local_max_index);
                } else {
                    eigenvalue_array.resize(local_max_index);
                }
                DMatrix eigvecs;
                
                DMatrix kpm_dos;
                DMatrix M_xx;
                DMatrix M_yy;
                DMatrix M_xy;
                
                computeEigenComplex(eigenvalue_array, eigvecs, kpm_dos, M_xx, M_yy, M_xy, H, dxH, dyH, jobIn, current_index_reduction, local_max_index);
                
                std::vector< std::vector<double> > eigenweights;
                std::vector< std::vector<double> > eigenvectors_r;
                std::vector< std::vector<double> > eigenvectors_c;
                std::vector< std::vector<double> > cond_xx;
                std::vector< std::vector<double> > cond_yy;
                std::vector< std::vector<double> > cond_xy;
                
                std::vector<double> eigenvalue_real_array;
                eigenvalue_real_array.resize(band_end_idx - band_start_idx);
                
                for (int i = band_start_idx; i < band_end_idx; ++i){
                    eigenvalue_real_array[i - band_start_idx] = eigenvalue_array[i];
                    //printf("eigenvalues[%d] = %lf \n",i,eigenvalue_array[i]);
                }
                
                results_out.setParam("eigenvalues",eigenvalue_real_array);
                
                std::complex<double>* eigenvector_array;
                eigenvector_array = eigvecs.getCpxValPtr();
                
                if (d_kpm_dos == 1){
                    
                    std::vector<double> kpm_dos_out;
                    std::vector< std::vector<double> > kpm_dos_in;
                    
                    double* val_kpm_dos;
                    val_kpm_dos = kpm_dos.getValPtr();
                    
                    for (int i = 0; i < poly_order; ++i){
                        std::vector<double> temp_vec;
                        
                        for (int j = 0; j < poly_order; ++j) {
                            temp_vec.push_back(val_kpm_dos[i + j*poly_order]);
                            if (i == 0){
                                printf("temp_vec[%d] = %lf \n",j,temp_vec[j]);
                            }
                        }
                        
                        kpm_dos_in.push_back(temp_vec);
                    }
                    
                    /*
                     Job_params dummy_result(results_out);
                     dummy_result.setParam("kpm_dos_mat",kpm_dos_in);
                     Param_tools::matrixResponseTransform(dummy_result, "kpm_dos_mat");
                     kpm_dos_in = dummy_result.getDoubleMat("kpm_dos_mat");
                     */
                    results_out.setParam("kpm_dos_mat",kpm_dos_in);
                    Param_tools::matrixResponseTransform(results_out, "kpm_dos_mat");
                    kpm_dos_in = results_out.getDoubleMat("kpm_dos_mat");
                    
                    
                    for (int i = 0; i < poly_order; ++i){
                        kpm_dos_out.push_back(kpm_dos_in[i][i]);
                        printf("kpm_dos_out[%d] = %lf \n",i,kpm_dos_out[i]);
                    }
                    
                    results_out.setParam("kpm_dos",kpm_dos_out);
                    
                }
                
                if (d_weights == 1){
                    
                    for (int t = 0; t < num_targets; ++t){
                        
                        int tar_here = target_list[t];
                        std::vector<double> temp_vec;
                        
                        for (int i = band_start_idx; i < band_end_idx; ++i){
                            std::complex<double> w1 = eigenvector_array[i*local_max_index + tar_here];
                            std::complex<double> w2 = w1*std::conj(w1);
                            temp_vec.push_back(w2.real());
                        }
                        
                        eigenweights.push_back(temp_vec);
                    }
                    
                    results_out.setParam("eigenweights",eigenweights);
                    
                }
                
                if (d_vecs == 1){
                    
                    for (int j = band_start_idx; j < band_end_idx; ++j){
                        
                        std::vector<double> temp_vec_r;
                        std::vector<double> temp_vec_c;
                        
                        for (int i = 0; i < local_max_index; ++i){
                            double w = eigenvector_array[i + j*local_max_index].real();
                            temp_vec_r.push_back(w);
                            double v = eigenvector_array[i + j*local_max_index].imag();
                            temp_vec_c.push_back(v);
                        }
                        
                        eigenvectors_r.push_back(temp_vec_r);
                        eigenvectors_c.push_back(temp_vec_c);
                    }
                    
                    results_out.setParam("eigenvectors_r",eigenvectors_r);
                    results_out.setParam("eigenvectors_c",eigenvectors_c);
                    
                }
                
                if (d_cond > 0){
                    
                    std::complex<double>* val_M_xx;
                    val_M_xx = M_xx.getCpxValPtr();
                    
                    double bz_area = 1.0;
                    if (boundary_condition == 1){
                        // need to normalize via BZ area...
                        std::vector< std::vector<double> > supercell = opts.getDoubleMat("supercell");
                        bz_area = Param_tools::computeReciprocalArea(supercell);
                    }
                    
                    for (int i = 0; i < poly_order; ++i){
                        std::vector<double> temp_xx;
                        
                        for (int j = 0; j < poly_order; ++j) {
                            temp_xx.push_back(val_M_xx[i + j*poly_order].real());
                        }
                        
                        cond_xx.push_back(temp_xx);
                        
                    }
                    
                    results_out.setParam("M_xx", cond_xx);
                    
                    if (d_cond > 1){
                        
                        std::complex<double>* val_M_yy;
                        val_M_yy = M_yy.getCpxValPtr();
                        
                        std::complex<double>* val_M_xy;
                        val_M_xy = M_xy.getCpxValPtr();
                        
                        for (int i = 0; i < poly_order; ++i){
                            std::vector<double> temp_yy;
                            std::vector<double> temp_xy;
                            
                            for (int j = 0; j < poly_order; ++j) {
                                temp_yy.push_back(val_M_yy[i + j*poly_order].real());
                                temp_xy.push_back(val_M_xy[i + j*poly_order].real());
                            }
                            
                            cond_yy.push_back(temp_yy);
                            cond_xy.push_back(temp_xy);
                            
                            results_out.setParam("M_yy", cond_yy);
                            results_out.setParam("M_xy", cond_xy);
                        }
                    }
                    
                    
                    // convert to E basis
                    Param_tools::conductivityTransform(results_out);
                    
                }
                
                if (chiral_on == 1){
                    
                    DMatrix J_B_x_cheb;
                    DMatrix J_B_y_cheb;
                    DMatrix J_T_x_cheb;
                    DMatrix J_T_y_cheb;
                    DMatrix J_TB_x_cheb;
                    DMatrix J_TB_y_cheb;
                    
                    // Project results onto chebyshev basis
                    computeMatrixResponse(jobIn, eigenvalue_array, eigvecs, J_B_x, J_B_x_cheb);
                    computeMatrixResponse(jobIn, eigenvalue_array, eigvecs, J_B_y, J_B_y_cheb);
                    computeMatrixResponse(jobIn, eigenvalue_array, eigvecs, J_T_x, J_T_x_cheb);
                    computeMatrixResponse(jobIn, eigenvalue_array, eigvecs, J_T_y, J_T_y_cheb);
                    computeMatrixResponse(jobIn, eigenvalue_array, eigvecs, J_TB_x, J_TB_x_cheb);
                    computeMatrixResponse(jobIn, eigenvalue_array, eigvecs, J_TB_y, J_TB_y_cheb);
                    
                    std::vector< std::vector< std::complex<double> > > J_B_x_mat;
                    std::vector< std::vector< std::complex<double> > > J_B_y_mat;
                    std::vector< std::vector< std::complex<double> > > J_T_x_mat;
                    std::vector< std::vector< std::complex<double> > > J_T_y_mat;
                    std::vector< std::vector< std::complex<double> > > J_TB_x_mat;
                    std::vector< std::vector< std::complex<double> > > J_TB_y_mat;
                    
                    // get chebyshev results as CPP nested vector matrices
                    J_B_x_cheb.getCPPValCopy(J_B_x_mat);
                    J_B_y_cheb.getCPPValCopy(J_B_y_mat);
                    J_T_x_cheb.getCPPValCopy(J_T_x_mat);
                    J_T_y_cheb.getCPPValCopy(J_T_y_mat);
                    J_TB_x_cheb.getCPPValCopy(J_TB_x_mat);
                    J_TB_y_cheb.getCPPValCopy(J_TB_y_mat);
                    
                    
                    // save results
                    results_out.setParam("J_B_x", J_B_x_mat);
                    results_out.setParam("J_B_y", J_B_y_mat);
                    results_out.setParam("J_T_x", J_T_x_mat);
                    results_out.setParam("J_T_y", J_T_y_mat);
                    results_out.setParam("J_TB_x", J_TB_x_mat);
                    results_out.setParam("J_TB_y", J_TB_y_mat);
                    
                    
                    Param_tools::conductivityTransform(results_out);
                }
            }
        }
        
        
        // Save time at which solver finished
        
        time_t tempEnd;
        time(&tempEnd);
        solverTimes.push_back(tempEnd);
        
        clock_t timeSolve;
        timeSolve = clock();
        Param_tools::saveTiming(results_out, ( ((double)timeSolve) - ((double)timeMat) )/((double)CLOCKS_PER_SEC), "SOLVER");
        // send back work to root, with a trash value sent first to get picked up by the recvSpool
        int trash = 1;
        MPI::COMM_WORLD.Send(
                             &trash,             // input buffer
                             1,                    // size of buffer
                             MPI::INT,            // type of buffer
                             root,                // rank to receive
                             0);                    // MPI label
        
        results_out.sendParams(root);
        //if (rank == print_rank)
        //printf("rank %d finished 1 job! \n", rank);
        
        // Cleanup allocated memory
        delete[] alpha_0_x_arr;
        delete[] alpha_0_y_arr;
        
        // End of while(1) means we wait for another instruction from root
    }
    
}

void Locality::setConfigPositions(double* i2pos, double* index_to_pos, int* index_to_grid, std::vector< std::vector<double> >& new_shift_configs, std::vector< std::vector<double> >& shift_configs, std::vector< std::vector< std::vector<double> > > &strain, Job_params jobIn){
    
    int solver_type = jobIn.getInt("solver_type");
    int solver_space = jobIn.getInt("solver_space");
    int strain_type = jobIn.getInt("strain_type");
    int gsfe_z_on = jobIn.getInt("gsfe_z_on");
    int global_shifts_on = jobIn.getInt("global_shifts_on");
    
    std::vector< std::vector<double> > shifts = jobIn.getDoubleMat("shifts");
    
    if (strain_type == 1){
        strainInfo.loadFourierConfigFile(jobIn.getString("strain_file"));
        strainInfo.setOpts(jobIn);
        strain.resize(max_index);
    }
    
    if (strain_type == 2){
        new_shift_configs.resize(max_index);
        strainInfo.loadConfigFile(jobIn.getString("strain_file"));
        strainInfo.setOpts(jobIn);
        strain.resize(max_index);
    }
    
    if (strain_type == 3){
        strainInfo.setOpts(jobIn);
        strain.resize(max_index);
    }
    
    if (strain_type == 5){
        printf("Entering strain_type = PLANEWAVES \n");
        strainInfo.loadFourierConfigFile_interp(jobIn.getString("strain_thetas"),jobIn.getString("strain_x_coeffs"),jobIn.getString("strain_y_coeffs"),jobIn.getString("strain_z_coeffs"));
        strainInfo.setOpts(jobIn);
        strain.resize(max_index);
    }
    
    if (solver_space == 0){
        for (int i = 0; i < max_index; ++i) {
            
            int orbit = index_to_grid[i*4 + 2];
            int s = index_to_grid[i*4 + 3];
            
            double s1_a_base[2][2];
            double s1_a[2][2];
            
            for (int m = 0; m < 2; ++m){
                for (int n = 0; n < 2; ++n){
                    s1_a_base[m][n] = sdata[s].a[m][n];
                }
            }
            
            double theta = angles[s];
            for (int d = 0; d < 2; ++d){
                // d tells us which unit-cell vec we are rotating by theta
                s1_a[d][0] = cos(theta)*s1_a_base[d][0] - sin(theta)*s1_a_base[d][1];
                s1_a[d][1] = sin(theta)*s1_a_base[d][0] + cos(theta)*s1_a_base[d][1];
            }
            
            i2pos[i*3 + 0] = index_to_pos[i*3 + 0] + shifts[s][0]*s1_a[0][0] + shifts[s][1]*s1_a[1][0];
            i2pos[i*3 + 1] = index_to_pos[i*3 + 1] + shifts[s][0]*s1_a[0][1] + shifts[s][1]*s1_a[1][1];
            i2pos[i*3 + 2] = index_to_pos[i*3 + 2];
            
            // The global shift routine has been moved to sheet.cpp
            /*
             if (global_shifts_on == 1){
             std::vector< std::vector<double> > g_shifts = jobIn.getDoubleMat("global_shifts");
             for (int d = 0; d < 3; ++d){
             if (d < 2){
             i2pos[i*3 + d] = i2pos[i*3 + d] + g_shifts[s][0]*s1_a[0][d] + g_shifts[s][1]*s1_a[1][d];
             } else {
             // no z-update here
             }
             }
             }
             */
            if (strain_type == 1){
                
                // sample strain from a supercell grid
                std::vector< std::vector<double> > sc = jobIn.getDoubleMat("supercell");
                
                std::vector< std::vector<double> > sc_inv;
                sc_inv.resize(2);
                sc_inv[0].resize(2);
                sc_inv[1].resize(2);
                
                double det = sc[0][0]*sc[1][1] - sc[0][1]*sc[1][0];
                sc_inv[0][0] =  sc[1][1]/det;
                sc_inv[0][1] = -sc[1][0]/det;
                sc_inv[1][0] = -sc[0][1]/det;
                sc_inv[1][1] =  sc[0][0]/det;
                
                double x = index_to_pos[i*3 + 0] + shifts[s][0]*s1_a[0][0] + shifts[s][1]*s1_a[1][0];
                double y = index_to_pos[i*3 + 1] + shifts[s][0]*s1_a[0][1] + shifts[s][1]*s1_a[1][1];;
                double z = index_to_pos[i*3 + 2];
                
                std::vector<double> sc_pos;
                sc_pos.resize(2);
                sc_pos[0] = sc_inv[0][0]*x + sc_inv[0][1]*y;
                sc_pos[1] = sc_inv[1][0]*x + sc_inv[1][1]*y;
                
                // For fourier strain relaxation method
                double b1[2];
                double b2[2];
                
                std::vector< std::vector<double> > a_sc = opts.getDoubleMat("supercell");
                std::vector< std::vector<double> > b_vec = getReciprocal(a_sc);
                
                b1[0] = b_vec[0][0];
                b1[1] = b_vec[0][1];
                b2[0] = b_vec[1][0];
                b2[1] = b_vec[1][1];
                
                double r[2];
                r[0] = x;
                r[1] = y;
                
                std::vector<double> disp_here = strainInfo.fourierStrainDisp_sc(r, b1, b2, s);
                //printf("disp_here = [%lf, %lf, %lf] \n",disp_here[0], disp_here[1], disp_here[2]);
                
                // old Supercell methods...
                //std::vector<double> disp_here = strainInfo.supercellDisp(sc_pos, s, orbit);
                //strain[i] = strainInfo.supercellStrain(sc_pos, s, orbit);
                
                i2pos[i*3 + 0] = i2pos[i*3 + 0] + disp_here[0];
                i2pos[i*3 + 1] = i2pos[i*3 + 1] + disp_here[1];
                i2pos[i*3 + 2] = i2pos[i*3 + 2] + disp_here[2];
                
            } else if (strain_type == 2){
                // sample strain from the space of configurations
                
                new_shift_configs[i].resize(2);
                
                if (s == 0){
                    new_shift_configs[i][0] = shift_configs[i][0]
                    + (cos(angles[0])*shifts[0][0] - sin(angles[0])*shifts[0][1])
                    - (cos(angles[1])*shifts[1][0] - sin(angles[1])*shifts[1][1]);
                    new_shift_configs[i][1] = shift_configs[i][1]
                    + (sin(angles[0])*shifts[0][0] + cos(angles[0])*shifts[0][1])
                    - (sin(angles[1])*shifts[1][0] + cos(angles[1])*shifts[1][1]);
                } else if (s == 1){
                    new_shift_configs[i][0] = shift_configs[i][0]
                    - (cos(angles[0])*shifts[0][0] - sin(angles[0])*shifts[0][1])
                    + (cos(angles[1])*shifts[1][0] - sin(angles[1])*shifts[1][1]);
                    new_shift_configs[i][1] = shift_configs[i][1]
                    - (sin(angles[0])*shifts[0][0] + cos(angles[0])*shifts[0][1])
                    + (sin(angles[1])*shifts[1][0] + cos(angles[1])*shifts[1][1]);
                }
                
                new_shift_configs[i][0] = fmod(new_shift_configs[i][0],1.0);
                new_shift_configs[i][1] = fmod(new_shift_configs[i][1],1.0);
                
                while (new_shift_configs[i][0] < 0.0){
                    new_shift_configs[i][0] = new_shift_configs[i][0] + 1.0;
                }
                
                while (new_shift_configs[i][1] < 0.0){
                    new_shift_configs[i][1] = new_shift_configs[i][1] + 1.0;
                }
                
                std::vector<double> disp_here;
                disp_here.resize(3);
                
                // Fourier method, not very good for small angles!
                /*
                 //printf("on k=%d, sheet=%d, orbit=%d\n",i,s,orbit);
                 std::vector<double> temp_disp = strainInfo.fourierStrainDisp_old(new_shift_configs[i], s, orbit);
                 double ang = angles[s];
                 ang = 0.0;
                 
                 disp_here[0] = cos(ang)*temp_disp[0] - sin(ang)*temp_disp[1];
                 disp_here[1] = sin(ang)*temp_disp[0] + cos(ang)*temp_disp[1];
                 disp_here[2] = temp_disp[2];
                 */
                
                // Config sampling method, fast but may break symmetries
                // Tries to force rotational symm for numerical data
                // /*
                int rot_max = 2;
                for (int r = 1; r < rot_max; ++r){
                    double x_config = new_shift_configs[i][0] + new_shift_configs[i][1]*0.5;
                    double y_config = new_shift_configs[i][1]*sqrt(3)/2.0;
                    
                    double ang = ((double) r)*2.0*PI/3.0 ;
                    
                    double rot_x = x_config*cos(ang) - y_config*sin(ang);
                    double rot_y = x_config*sin(ang) + y_config*cos(ang);
                    
                    double new_config_x = rot_x - rot_y/sqrt(3);
                    double new_config_y = rot_x*2.0/sqrt(3);
                    
                    while (new_config_x >= 1)
                        new_config_x = new_config_x - 1.0;
                    while (new_config_x < 0)
                        new_config_x = new_config_x + 1.0;
                    while (new_config_y >= 1)
                        new_config_y = new_config_y - 1.0;
                    while (new_config_y < 0)
                        new_config_y = new_config_y + 1.0;
                    
                    std::vector<double> rot_config;
                    rot_config.resize(2);
                    rot_config[0] = new_config_x;
                    rot_config[1] = new_config_y;
                    std::vector<double> new_disp = strainInfo.interpStrainDisp(rot_config, s, orbit);
                    double rot_disp_x = new_disp[0]*cos(-ang) - new_disp[1]*sin(-ang);
                    double rot_disp_y = new_disp[0]*sin(-ang) + new_disp[1]*cos(-ang);
                    
                    disp_here[0] = disp_here[0] + rot_disp_x;
                    disp_here[1] = disp_here[1] + rot_disp_y;
                }
                disp_here[0] = disp_here[0]/3.0;
                disp_here[1] = disp_here[1]/3.0;
                // */
                
                i2pos[i*3 + 0] = i2pos[i*3 + 0] + disp_here[0];
                i2pos[i*3 + 1] = i2pos[i*3 + 1] + disp_here[1];
                i2pos[i*3 + 2] = i2pos[i*3 + 2] + disp_here[2];
                
                if (gsfe_z_on == 1){
                    
                    int s_here = index_to_grid[i*4 + 3];
                    Materials::Mat mat = sdata[s_here].mat;
                    std::array<std::array<double, 2>, 2> lattice = Materials::lattice(mat);
                    
                    std::array<std::array<double, 2>, 2> lattice_inv;
                    double l_det = lattice[0][0]*lattice[1][1] - lattice[0][1]*lattice[1][0];
                    lattice_inv[0][0] =  lattice[1][1]/l_det;
                    lattice_inv[0][1] = -lattice[1][0]/l_det;
                    lattice_inv[1][0] = -lattice[0][1]/l_det;
                    lattice_inv[1][1] =  lattice[0][0]/l_det;
                    
                    new_shift_configs[i][0] = new_shift_configs[i][0] + disp_here[0]*lattice_inv[0][0] + disp_here[1]*lattice_inv[1][0];
                    new_shift_configs[i][1] = new_shift_configs[i][1] + disp_here[0]*lattice_inv[0][1] + disp_here[1]*lattice_inv[1][1];
                    
                    double dz = strainInfo.gsfeHeight(new_shift_configs[i]);
                    if (s_here == 0){
                        i2pos[i*3 + 2] = i2pos[i*3 + 2] - dz/2.0;
                    } else if (s_here == 1){
                        i2pos[i*3 + 2] = i2pos[i*3 + 2] + dz/2.0;
                    }
                    
                    
                }
                
            } else if (strain_type == 3){
                //sample strain by actual realspace position
                
                std::vector<double> pos_in;
                pos_in.resize(3);
                pos_in[0] = i2pos[i*3 + 0];
                pos_in[1] = i2pos[i*3 + 1];
                pos_in[2] = i2pos[i*3 + 2];
                
                std::vector<double> disp_here = strainInfo.realspaceDisp(pos_in, s, orbit);
                strain[i] = strainInfo.realspaceStrain(pos_in, s, orbit);
                
                i2pos[i*3 + 0] = i2pos[i*3 + 0] + disp_here[0];
                i2pos[i*3 + 1] = i2pos[i*3 + 1] + disp_here[1];
                i2pos[i*3 + 2] = i2pos[i*3 + 2] + disp_here[2];
                
            }
            
        }
    } else if (solver_space == 1){
        
        // Here we define new positions as " A = q + K "
        
        for (int i = 0; i < max_index; ++i) {
            
            //int s = index_to_grid[i*4 + 3];
            int s = 0;
            //if (s == 0){
            i2pos[i*3 + 0] = index_to_pos[i*3 + 0] + shifts[s][0];
            i2pos[i*3 + 1] = index_to_pos[i*3 + 1] + shifts[s][1];
            i2pos[i*3 + 2] = index_to_pos[i*3 + 2];
            
            //printf("pos[%d] = [%lf, %lf, %lf] \n",i,i2pos[i*3 + 0],i2pos[i*3 + 1],i2pos[i*3 + 2]);
            /**} if (s == 1){
             i2pos[i*3 + 0] = -index_to_pos[i*3 + 0] + shifts[s*3 + 0]*b1[0][0] + shifts[s*3 + 1]*b1[0][1];
             i2pos[i*3 + 1] = -index_to_pos[i*3 + 1] + shifts[s*3 + 0]*b1[0][1] + shifts[s*3 + 1]*b1[1][1];
             i2pos[i*3 + 2] = -index_to_pos[i*3 + 2];
             }
             **/
            
        }
        
        
    }
    
}

std::vector<double> Locality::getConfigDisp(std::vector<double> config_in, int s){
    
    std::vector<double> disp;
    disp.resize(2);
    disp[0] = 0;
    disp[1] = 0;
    
    return disp;
}

void Locality::generateRealH(SpMatrix &H, SpMatrix &dxH, SpMatrix &dyH, double* alpha_0_x_arr, double* alpha_0_y_arr, Job_params jobIn, int* index_to_grid, double* i2pos,
                             int* inter_pairs, std::vector< std::vector<int> > inter_sc_vecs, int* intra_pairs, double* intra_pairs_t,
                             std::vector< std::vector<int> > intra_sc_vecs, std::vector< std::vector< std::vector<double> > > strain, std::vector<int> current_index_reduction, int local_max_index){
    
    int intra_searchsize = opts.getInt("intra_searchsize");
    
    int solver_type = jobIn.getInt("solver_type");
    int observable_type = jobIn.getInt("observable_type");
    int boundary_condition = jobIn.getInt("boundary_condition");
    int strain_type = jobIn.getInt("strain_type");
    
    int diagonalize = jobIn.getInt("diagonalize");
    int elecOn = jobIn.getInt("elecOn");
    double E = jobIn.getDouble("E");
    double energy_rescale = jobIn.getDouble("energy_rescale");
    double energy_shift = jobIn.getDouble("energy_shift");
    //std::vector<double> k_vec = jobIn.getDoubleVec("k");
    
    //jobIn.printParams();
    
    
    // Indexes how many inter terms we have entered so far
    int inter_counter = 0;
    
    // Indexes how many intra terms we have entered so far
    int intra_counter = 0;
    
    // Total number of expected non-zero matrix elements
    int max_nnz = max_intra_pairs + max_inter_pairs;
    
    // Sparse matrix format is 2 arrays with length = nnz, and 1 array with length = max_index + 1
    
    // col_index tells us the col of element i
    // v tells us the value of element i
    // row_pointer tells us the start of row j (at element i = row_pointer[j]) and the end of row j (at element i = row_pointer[j] - 1)
    
    H.setup(max_nnz, local_max_index, local_max_index);
    dxH.setup(max_nnz, local_max_index, local_max_index);
    dyH.setup(max_nnz, local_max_index, local_max_index);
    
    
    int* col_index = H.allocColIndx();
    int* row_pointer = H.allocRowPtr();
    
    int* col_index_dx = dxH.allocColIndx();
    int* row_pointer_dx = dxH.allocRowPtr();
    
    int* col_index_dy = dyH.allocColIndx();
    int* row_pointer_dy = dyH.allocRowPtr();
    
    double* v;
    double* v_dx;
    double* v_dy;
    
    v = H.allocRealVal();
    v_dx = dxH.allocRealVal();
    v_dy = dyH.allocRealVal();
    
    
    // Count the current element, i.e. "i = input_counter"
    int input_counter = 0;
    
    // Loop through every orbital (i.e. rows of H)
    for (int k_i = 0; k_i < max_index; ++k_i){
        
        double t;
        int skip_here1 = 0;
        
        if (solver_type == 3 || solver_type == 4){
            if (current_index_reduction[k_i] + 1 == current_index_reduction[k_i + 1]){
                skip_here1 = 1;
            }
        }
        
        int k = k_i - current_index_reduction[k_i];
        
        // Save starting point of row k
        row_pointer[k] = input_counter;
        row_pointer_dx[k] = input_counter;
        row_pointer_dy[k] = input_counter;
        
        // While we are still at the correct index in our intra_pairs list:
        bool same_index1 = true;
        while(same_index1) {
            
            int skip_here2 = 0;
            
            // if the first index of intra_pairs changes, we stop
            if ((intra_counter) == max_intra_pairs || intra_pairs[intra_counter*2 + 0] != k_i) {
                same_index1 = false;
                continue; // go to "while(same_index1)" which will end this loop
            }
            
            int new_k = intra_pairs[intra_counter*2 + 1];
            
            // if we are accounting for defects, we check if the other half of the pair has been removed
            if (solver_type == 3 || solver_type == 4){
                if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
                    skip_here2 = 1;
                }
            }
            
            // we save this pair into our sparse matrix format
            if (skip_here1 == 0 && skip_here2 == 0){
                
                int sc_i = 0;
                int sc_j = 0;
                double new_pos_shift_x = 0.0;
                double new_pos_shift_y = 0.0;
                
                if (boundary_condition == 1){
                    
                    int sheet_here = index_to_grid[k_i*4 + 3];
                    std::vector< std::vector<double> > supercell_orig = sdata[sheet_here].supercell;
                    double theta = angles[sheet_here];
                    std::vector< std::vector<double> > supercell;
                    supercell.resize(2);
                    
                    for (int dim = 0; dim < 2; ++dim){
                        supercell[dim].resize(2);
                        
                        supercell[dim][0] = cos(theta)*supercell_orig[dim][0] - sin(theta)*supercell_orig[dim][1];
                        supercell[dim][1] = sin(theta)*supercell_orig[dim][0] + cos(theta)*supercell_orig[dim][1];
                    }
                    
                    
                    sc_i = intra_sc_vecs[intra_counter][0];
                    sc_j = intra_sc_vecs[intra_counter][1];
                    new_pos_shift_x = sc_i*supercell[0][0] + sc_j*supercell[1][0];
                    new_pos_shift_y = sc_i*supercell[0][1] + sc_j*supercell[1][1];
                }
                
                // get positions of orbitals
                double x1 = i2pos[k_i*3 + 0] + new_pos_shift_x;
                double y1 = i2pos[k_i*3 + 1] + new_pos_shift_y;
                double z1 = i2pos[k_i*3 + 2];
                double x2 = i2pos[new_k*3 + 0];
                double y2 = i2pos[new_k*3 + 1];
                double z2 = i2pos[new_k*3 + 2];
                
                double raw_t;
                
                if (strain_type == 1){
                    
                    int i0 = index_to_grid[k_i*4 + 0];
                    int j0 = index_to_grid[k_i*4 + 1];
                    int l0 = index_to_grid[k_i*4 + 2];
                    int s0 = index_to_grid[k_i*4 + 3];
                    
                    int ih = index_to_grid[new_k*4 + 0];
                    int jh = index_to_grid[new_k*4 + 1];
                    int lh = index_to_grid[new_k*4 + 2];
                    int sh = index_to_grid[new_k*4 + 3];
                    
                    double xh = i2pos[new_k*3 + 0];
                    double yh = i2pos[new_k*3 + 1];
                    double zh = i2pos[new_k*3 + 2];
                    
                    std::array<int,2> grid_disp = {{
                        ih-i0,
                        jh-j0 }};
                    
                    // we correct the grid values by the supercell_stride when there are periodic BCs
                    if (boundary_condition == 1){
                        grid_disp[0] = grid_disp[0] - intra_sc_vecs[intra_counter][0]*sdata[s0].supercell_stride[0][0] - intra_sc_vecs[intra_counter][1]*sdata[s0].supercell_stride[1][0];
                        grid_disp[1] = grid_disp[1] - intra_sc_vecs[intra_counter][0]*sdata[s0].supercell_stride[0][1] - intra_sc_vecs[intra_counter][1]*sdata[s0].supercell_stride[1][1];
                    }
                    
                    // take strain as avg of both orbital strains
                    std::vector< std::vector<double> > strain_here;
                    strain_here.resize(2);
                    
                    for (int i = 0; i < 2; ++i){
                        strain_here[i].resize(2);
                        for (int j = 0; j < 2; ++j){
                            // for now we turn off the strain term... it is not implemented for Fourier samp yet.
                            //strain_here[i][j] = (strain[k_i][i][j]  + strain[new_k][i][j])/2.0;
                            strain_here[i][j] = 0.0;
                        }
                    }
                    
                    // we need to find the bonding angle relative the u_ij definitions:
                    std::array<double, 2> strain_dir;
                    
                    strain_dir[0] = xh - x1;
                    strain_dir[1] = yh - y1;
                    
                    double strain_dir_norm = sqrt(strain_dir[0]*strain_dir[0] + strain_dir[1]*strain_dir[1]);
                    
                    // angle (counter-clockwise) from a bonding direction of +x
                    double strain_theta = 0;
                    
                    if (strain_dir_norm != 0){
                        strain_dir[0] = strain_dir[0]/strain_dir_norm;
                        strain_dir[1] = strain_dir[1]/strain_dir_norm;
                        
                        if (strain_dir[0] == 0.0){
                            if (strain_dir[1] > 0){
                                strain_theta = PI_2;
                            } else {
                                strain_theta = -PI_2;
                            }
                        } else {
                            double ratio = strain_dir[1]/strain_dir[0];
                            strain_theta = atan(ratio);
                            if (strain_dir[0] < 0){
                                strain_theta = strain_theta + PI;
                            }
                        }
                    }
                    
                    // now rotate;
                    std::vector< std::vector<double> > strain_rot;
                    strain_rot.resize(2);
                    strain_rot[0].resize(2);
                    strain_rot[1].resize(2);
                    
                    strain_rot[0][0] =  strain_here[0][0]*cos(strain_theta)*cos(strain_theta) +
                    strain_here[1][1]*sin(strain_theta)*sin(strain_theta) +
                    strain_here[0][1]*sin(strain_theta)*cos(strain_theta);
                    strain_rot[1][1] =  strain_here[0][0]*sin(strain_theta)*sin(strain_theta) +
                    strain_here[1][1]*cos(strain_theta)*cos(strain_theta) +
                    -2*strain_here[0][1]*sin(strain_theta)*cos(strain_theta);
                    strain_rot[0][1] =  (strain_here[1][1] - strain_here[0][0])*sin(strain_theta)*cos(strain_theta) +
                    strain_here[0][1]*(cos(strain_theta)*cos(strain_theta) - sin(strain_theta)*sin(strain_theta));
                    strain_rot[1][0] = strain_rot[0][1];
                    
                    Materials::Mat mat = sdata[s0].mat;
                    raw_t = Materials::intralayer_term(l0, lh, grid_disp, strain_rot, mat);
                    
                } else {
                    
                    raw_t = intra_pairs_t[intra_counter];
                    
                }
                
                // if it is the diagonal element, we "shift" the matrix up or down in energy scale (to make sure the spectrum fits in [-1,1] for the Chebyshev method)
                // Also, if electric field is included (elecOn == 1) we add in an on-site offset due to this gate voltage.
                if (new_k == k_i){
                    if (elecOn == 1){
                        t = (raw_t + energy_shift + onSiteE(x1,y1,z1,E))/energy_rescale;
                    } else if (elecOn == 0){
                        t = (raw_t + energy_shift)/energy_rescale;
                    }
                    // Otherwise we enter the value just with rescaling
                }
                else {
                    t = raw_t/energy_rescale;
                }
                
                if (t != 0){
                    //printf("intra coupling (%d, %d) [%lf, %lf, %lf] -> [%lf, %lf, %lf] (%lf, %lf) = %lf \n",k_i,new_k,x1,y1,z1,x2,y2,z2,new_pos_shift_x,new_pos_shift_y,t)
                    //printf("intra_t [%d,%d] = %lf \n",k_i,new_k,t);
                    v[input_counter] += t;
                    v_dx[input_counter] += (x1 - x2)*t; // (delta_x)*t
                    v_dy[input_counter] += (y1 - y2)*t; // (delta_y)*t
                }
                
                // check if next pair is identical (possible with periodic wrapping), and save
                if ((intra_counter+1) == max_intra_pairs || k_i != intra_pairs[(intra_counter+1)*2 + 0] || new_k != intra_pairs[(intra_counter+1)*2 + 1]){
                    col_index[input_counter] = new_k - current_index_reduction[new_k];
                    col_index_dx[input_counter] =  new_k - current_index_reduction[new_k];
                    col_index_dy[input_counter] =  new_k - current_index_reduction[new_k];
                    ++input_counter;
                }
                
            }
            ++intra_counter;
        }
        
        // reset t, shouldn't be necessary but just in case!
        t = 0;
        
        // While we are still at the correct index in our inter_pairs list:
        bool same_index2 = true;
        while(same_index2) {
            
            double t;
            int skip_here2 = 0;
            
            // if the first index of intra_pairs changes, we stop
            if ((inter_counter) == max_inter_pairs ||  inter_pairs[inter_counter*2 + 0] != k_i) {
                same_index2 = false;
                continue; // go to "while(same_index2)" which will end this loop
            }
            
            int new_k = inter_pairs[inter_counter*2 + 1];
            
            // if we are accounting for defects, we check if the other half of the pair has been removed
            if (solver_type == 3 || solver_type == 4){
                if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
                    skip_here2 = 1;
                }
            }
            
            // we save this pair into our sparse matrix format
            if (skip_here1 == 0 && skip_here2 == 0){
                // get the index of the other orbital in this term
                
                int sc_i = 0;
                int sc_j = 0;
                double new_pos_shift_x = 0.0;
                double new_pos_shift_y = 0.0;
                
                if (boundary_condition == 1){
                    
                    std::vector< std::vector<double> > supercell = opts.getDoubleMat("supercell");
                    
                    sc_i = inter_sc_vecs[inter_counter][0];
                    sc_j = inter_sc_vecs[inter_counter][1];
                    
                    new_pos_shift_x = sc_i*supercell[0][0] + sc_j*supercell[1][0];
                    new_pos_shift_y = sc_i*supercell[0][1] + sc_j*supercell[1][1];
                }
                
                // get the position of both orbitals
                double x1 = i2pos[k_i*3 + 0] + new_pos_shift_x;
                double y1 = i2pos[k_i*3 + 1] + new_pos_shift_y;
                double z1 = i2pos[k_i*3 + 2];
                double x2 = i2pos[new_k*3 + 0];
                double y2 = i2pos[new_k*3 + 1];
                double z2 = i2pos[new_k*3 + 2];
                
                // and the orbit tag in their respective unit-cell
                int orbit1 = index_to_grid[k_i*4 + 2];
                int orbit2 = index_to_grid[new_k*4 + 2];
                
                Materials::Mat mat1 = sdata[index_to_grid[k_i*4 + 3]].mat;
                Materials::Mat mat2 = sdata[index_to_grid[new_k*4 + 3]].mat;
                
                // and the angle of the sheet each orbital is on
                double theta1;
                double theta2;
                if (strain_type == 0){
                    theta1 = angles[index_to_grid[k_i*4 + 3]];
                    theta2 = angles[index_to_grid[new_k*4 + 3]];
                } else {
                    // with strain, we need to compute the local theta by nearest neighbor bonding
                    // HARDCODED FOR GRAPHENE CURRENTLY
                    theta1 = getLocalTheta(k_i, index_to_grid, i2pos, max_index);
                    theta2 = getLocalTheta(new_k, index_to_grid, i2pos, max_index);
                }
                
                // use all this information to determine coupling energy
                // !!! Currently NOT generalized for materials other than graphene, need to do a material index call for both sheets and pass to a general "inter_coupling" method !!!
                
                std::array<double, 3> disp = {{
                    x2 - x1,
                    y2 - y1,
                    z2 - z1}};
                
                t = Materials::interlayer_term(orbit1, orbit2,disp, theta1, theta2,mat1, mat2)/energy_rescale;
                
                //if ( (new_k == 568 && k_i == 3941) || (new_k == 3941 && k_i == 568) ){
                //printf("inter[%d,%d] = %lf, [%lf, %lf] -> [%lf, %lf] \n", k_i,new_k,t, x1,y1, x2,y2);
                //}
                
                //t += interlayer_term(x1, y1, z1, x2, y2, z2, orbit1, orbit2, theta1, theta2, mat1, mat2)/energy_rescale;
                if (t != 0 ){
                    v[input_counter] += t;
                    v_dx[input_counter] += (x1 - x2 )*t;
                    v_dy[input_counter] += (y1 - y2)*t;
                }
                
                // check if next pair is identical (possible with periodic wrapping), or if we are at last element, to decide whether to save or not
                if ( (inter_counter+1) == max_inter_pairs || k_i != inter_pairs[(inter_counter+1)*2 + 0] || new_k != inter_pairs[(inter_counter+1)*2 + 1]){
                    col_index[input_counter] = new_k - current_index_reduction[new_k];
                    col_index_dx[input_counter] =  new_k - current_index_reduction[new_k];
                    col_index_dy[input_counter] =  new_k - current_index_reduction[new_k];
                    ++input_counter;
                }
                
            }
            ++inter_counter;
            
        }
        
    }
    
    // Save the end point + 1 of the last row
    row_pointer[local_max_index] = input_counter;
    row_pointer_dx[local_max_index] = input_counter;
    row_pointer_dy[local_max_index] = input_counter;
    
    int num_targets = jobIn.getInt("num_targets");
    std::vector<int> target_list = jobIn.getIntVec("target_list");
    
    if (diagonalize != 1 && observable_type == 1){
        for (int t = 0; t < num_targets; ++t){
            //printf("target = %d, local_max_index = %d \n", target_list[t], local_max_index);
            
            int target_here = target_list[t] - current_index_reduction[target_list[t]];
            //printf("target_here = %d \n",target_here);
            
            double* target_vec = new double[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                target_vec[i] = 0;
            }
            
            target_vec[target_here] = 1;
            
            double temp_vec_x[local_max_index];
            double temp_vec_y[local_max_index];
            
            for (int i = 0; i < local_max_index; ++i){
                temp_vec_x[i] = 0;
                temp_vec_y[i] = 0;
            }
            
            dxH.vectorMultiply(target_vec,temp_vec_x,1,0);
            dyH.vectorMultiply(target_vec,temp_vec_y,1,0);
            
            // want < 0 | dxH, not dxH | 0 >, so need to use: Transpose( < 0 | dxH ) = - dxH | 0 >
            for (int i = 0; i < local_max_index; ++ i){
                alpha_0_x_arr[t*local_max_index + i] = -temp_vec_x[i];
                alpha_0_y_arr[t*local_max_index + i] = -temp_vec_y[i];
            }
            
            delete[] target_vec;
        }
    }
    
    // ------------------------------
    // Following saves Matrix to file
    //
    // Should only used for for 1-job processes, otherwise they will overwrite each other!
    
    int matrix_save = opts.getInt("matrix_save");
    if (matrix_save > 0){
        
        std::ofstream outFile;
        const char* extension = "_matrix.dat";
        outFile.open ((job_name + extension).c_str());
        
        for(int i = 0; i < local_max_index; ++i){
            int start_index = row_pointer[i];
            int stop_index = row_pointer[i+1];
            for(int j = start_index; j < stop_index; ++j){
                outFile << col_index[j] + 1 << ", " << i + 1 << ", " << v[j] << "\n";
            }
        }
        
        if (matrix_save > 1){
            
            outFile.close();
            
            std::ofstream outFile2;
            const char* extension2 = "_dxH_matrix.dat";
            outFile2.open ((job_name + extension2).c_str());
            
            for(int i = 0; i < local_max_index; ++i){
                int start_index = row_pointer_dx[i];
                int stop_index = row_pointer_dx[i+1];
                for(int j = start_index; j < stop_index; ++j){
                    outFile2 << col_index_dx[j] + 1 << ", " << i + 1 << ", " << v_dx[j] << "\n";
                }
            }
            
            outFile2.close();
        }
        //
        
        // End Matrix Save
        // ---------------
    }
    
    // ------------------------------
    // Following saves positions and pairings to file
    
    int matrix_pos_save = opts.getInt("matrix_pos_save");
    if (matrix_pos_save > 0){
        
        std::ofstream outFile3;
        const char* extension3 = "_pos.dat";
        outFile3.open ((job_name + extension3).c_str());
        outFile3 << fixed;
        outFile3.precision(12);
        
        for(int i = 0; i < local_max_index; ++i){
            
            double x = i2pos[i*3 + 0];
            double y = i2pos[i*3 + 1];
            double z = i2pos[i*3 + 2];
            int s =  index_to_grid[i*4 + 3];
            int o =  index_to_grid[i*4 + 2];
            outFile3 << x << "    " << y << "    " << z << "    " << s << "    " << o << "\n";
        }
        
        outFile3.close();
    }
    // ---------------
    if (matrix_pos_save > 1){
        
        std::ofstream outFile4;
        const char* extension4 = "_intra_pos.dat";
        outFile4.open ((job_name + extension4).c_str());
        
        for(int i = 0; i < max_intra_pairs; ++i){
            outFile4 <<
            i2pos[intra_pairs[2*i + 0]*3 + 0] << ", " << i2pos[intra_pairs[2*i + 0]*3 + 1] << ", " << i2pos[intra_pairs[2*i + 0]*3 + 2] << ", " <<
            i2pos[intra_pairs[2*i + 1]*3 + 0] << ", " << i2pos[intra_pairs[2*i + 1]*3 + 1] << ", " << i2pos[intra_pairs[2*i + 1]*3 + 2] << ", " <<
            intra_pairs[2*i + 0] << ", " << intra_pairs[2*i + 1] <<"\n";
        }
        
        outFile4.close();
        // ---------------
        std::ofstream outFile5;
        const char* extension5 = "_inter_pos.dat";
        outFile5.open ((job_name + extension5).c_str());
        
        for(int i = 0; i < max_inter_pairs; ++i){
            outFile5 <<
            i2pos[inter_pairs[2*i + 0]*3 + 0] << ", " << i2pos[inter_pairs[2*i + 0]*3 + 1] << ", " << i2pos[inter_pairs[2*i + 0]*3 + 2] << ", " <<
            i2pos[inter_pairs[2*i + 1]*3 + 0] << ", " << i2pos[inter_pairs[2*i + 1]*3 + 1] << ", " << i2pos[inter_pairs[2*i + 1]*3 + 2] << ", " <<
            inter_pairs[2*i + 0] << ", " << inter_pairs[2*i + 1] << "\n";
        }
        
        outFile5.close();
    }
    
    // End Matrix Save
    // ---------------
    
    
}

void Locality::generateCpxH(SpMatrix &H, SpMatrix &dxH, SpMatrix &dyH,
                            SpMatrix &J_B_x, SpMatrix &J_B_y, SpMatrix &J_T_x, SpMatrix &J_T_y, SpMatrix &J_TB_x, SpMatrix &J_TB_y,
                            double* alpha_0_x_arr, double* alpha_0_y_arr, Job_params jobIn, int* index_to_grid, double* i2pos,
                            int* inter_pairs, std::vector< std::vector<int> > inter_sc_vecs, int* intra_pairs, double* intra_pairs_t,
                            std::vector< std::vector<int> > intra_sc_vecs, std::vector< std::vector< std::vector<double> > > strain, std::vector<int> current_index_reduction, int local_max_index){
    
    int intra_searchsize = opts.getInt("intra_searchsize");
    
    int mat_from_file = opts.getInt("mat_from_file");
    
    int solver_type = jobIn.getInt("solver_type");
    int observable_type = jobIn.getInt("observable_type");
    int boundary_condition = jobIn.getInt("boundary_condition");
    int strain_type = opts.getInt("strain_type");
    
    int diagonalize = jobIn.getInt("diagonalize");
    int chiral_on = jobIn.getInt("chiral_on");
    int magOn = jobIn.getInt("magOn");
    int elecOn = jobIn.getInt("elecOn");
    double B = jobIn.getDouble("B");
    double E = jobIn.getDouble("E");
    double energy_rescale = jobIn.getDouble("energy_rescale");
    double energy_shift = jobIn.getDouble("energy_shift");
    
    // Check to see if we will do k sampling
    int k_sampling = jobIn.getInt("k_sampling");
    
    // Initialize "k" (k_vec) to 0, but update if k_sampling is turned on
    std::vector<double> k_vec;
    k_vec.resize(2);
    k_vec[0] = 0.0;
    k_vec[1] = 0.0;
    if (k_sampling == 1){
        k_vec = jobIn.getDoubleVec("k_vec");
        //printf("k_vec [%lf, %lf]\n",k_vec[0],k_vec[1]);
    }
    
    //
    /*
     double phi = 0.0;
     double r = 2.2;
     
     while(phi < 2*PI){
     Materials::Mat mat1 = sdata[0].mat;
     std::array<double, 3> disp;
     disp[0] = r*cos(phi);
     disp[1] = r*sin(phi);
     disp[2] = 3.015;
     int orbit1 = 5;
     int orbit2 = 9;
     double theta1 = 0.0;
     double theta2 = 0.0;
     double t = Materials::interlayer_term(orbit1, orbit2, disp, theta1, theta2, mat1, mat1);
     printf("phi = %lf, t = %lf \n",phi*180.0/PI,t);
     phi = phi + PI/12;
     
     }
     */
    /*
     std::array<double, 3> disp;
     disp[0] = 0.0000;
     disp[1] = 1.8373;
     disp[2] = 3.0128;
     for (int o1 = 8; o1 < 11; ++o1){
     printf("[");
     for (int o2 = 5; o2 < 8; ++o2){
     Materials::Mat mat1 = sdata[0].mat;
     double theta1 = 0.0;
     double theta2 = 0.0;
     double t = Materials::interlayer_term(o1, o2, disp, theta1, theta2, mat1, mat1);
     printf("%lf ,",t);
     }
     printf("]\n");
     }
     */
    //
    
    //jobIn.printParams();
    
    
    // Indexes how many inter terms we have entered so far
    int inter_counter = 0;
    
    // Indexes how many intra terms we have entered so far
    int intra_counter = 0;
    
    // Total number of expected non-zero matrix elements
    int max_nnz = max_intra_pairs + max_inter_pairs;
    
    // Sparse matrix format is 2 arrays with length = nnz, and 1 array with length = max_index + 1
    
    // col_index tells us the col of element i
    // v tells us the value of element i
    // row_pointer tells us the start of row j (at element i = row_pointer[j]) and the end of row j (at element i = row_pointer[j] - 1)
    
    H.setup(max_nnz, local_max_index, local_max_index);
    dxH.setup(max_nnz, local_max_index, local_max_index);
    dyH.setup(max_nnz, local_max_index, local_max_index);
    
    if (chiral_on == 1){
        J_B_x.setup(max_nnz, local_max_index, local_max_index);
        J_B_y.setup(max_nnz, local_max_index, local_max_index);
        J_T_x.setup(max_nnz, local_max_index, local_max_index);
        J_T_y.setup(max_nnz, local_max_index, local_max_index);
        J_TB_x.setup(max_nnz, local_max_index, local_max_index);
        J_TB_y.setup(max_nnz, local_max_index, local_max_index);
    }
    
    int* col_index = H.allocColIndx();
    int* row_pointer = H.allocRowPtr();
    
    int* col_index_dx = dxH.allocColIndx();
    int* row_pointer_dx = dxH.allocRowPtr();
    
    int* col_index_dy = dyH.allocColIndx();
    int* row_pointer_dy = dyH.allocRowPtr();
    
    int* col_index_J_B_x;
    int* col_index_J_B_y;
    int* col_index_J_T_x;
    int* col_index_J_T_y;
    int* col_index_J_TB_x;
    int* col_index_J_TB_y;
    
    int* row_pointer_J_B_x;
    int* row_pointer_J_B_y;
    int* row_pointer_J_T_x;
    int* row_pointer_J_T_y;
    int* row_pointer_J_TB_x;
    int* row_pointer_J_TB_y;
    
    std::complex<double>* v_c;
    std::complex<double>* v_c_dx;
    std::complex<double>* v_c_dy;
    std::complex<double>* v_c_J_B_x;
    std::complex<double>* v_c_J_B_y;
    std::complex<double>* v_c_J_T_x;
    std::complex<double>* v_c_J_T_y;
    std::complex<double>* v_c_J_TB_x;
    std::complex<double>* v_c_J_TB_y;
    
    
    v_c = H.allocCpxVal();
    v_c_dx = dxH.allocCpxVal();
    v_c_dy = dyH.allocCpxVal();
    
    if (chiral_on == 1){
        
        col_index_J_B_x = J_B_x.allocColIndx();
        col_index_J_B_y = J_B_y.allocColIndx();
        col_index_J_T_x = J_T_x.allocColIndx();
        col_index_J_T_y = J_T_y.allocColIndx();
        col_index_J_TB_x = J_TB_x.allocColIndx();
        col_index_J_TB_y = J_TB_y.allocColIndx();
        
        row_pointer_J_B_x = J_B_x.allocRowPtr();
        row_pointer_J_B_y = J_B_y.allocRowPtr();
        row_pointer_J_T_x = J_T_x.allocRowPtr();
        row_pointer_J_T_y = J_T_y.allocRowPtr();
        row_pointer_J_TB_x = J_TB_x.allocRowPtr();
        row_pointer_J_TB_y = J_TB_y.allocRowPtr();
        
        v_c_J_B_x = J_B_x.allocCpxVal();
        v_c_J_B_y = J_B_y.allocCpxVal();
        v_c_J_T_x = J_T_x.allocCpxVal();
        v_c_J_T_y = J_T_y.allocCpxVal();
        v_c_J_TB_x = J_TB_x.allocCpxVal();
        v_c_J_TB_y = J_TB_y.allocCpxVal();
    }
    
    // Count the current element, i.e. "i = input_counter"
    int input_counter = 0;
    
    // Loop through every orbital (i.e. rows of H)
    for (int k_i = 0; k_i < max_index; ++k_i){
        
        std::complex<double> t_cpx;
        t_cpx = 0;
        int skip_here1 = 0;
        
        if (solver_type == 3 || solver_type == 4){
            if (current_index_reduction[k_i] + 1 == current_index_reduction[k_i + 1]){
                skip_here1 = 1;
            }
        }
        
        int k = k_i - current_index_reduction[k_i];
        
        // Save starting point of row k
        row_pointer[k] = input_counter;
        row_pointer_dx[k] = input_counter;
        row_pointer_dy[k] = input_counter;
        
        if (chiral_on == 1){
            row_pointer_J_B_x[k] = input_counter;
            row_pointer_J_B_y[k] = input_counter;
            row_pointer_J_T_x[k] = input_counter;
            row_pointer_J_T_y[k] = input_counter;
            row_pointer_J_TB_x[k] = input_counter;
            row_pointer_J_TB_y[k] = input_counter;
        }
        
        // While we are still at the correct index in our intra_pairs list:
        bool same_index1 = true;
        while(same_index1) {
            
            int skip_here2 = 0;
            
            // if the first index of intra_pairs changes, we stop
            if (intra_counter == max_intra_pairs || intra_pairs[intra_counter*2 + 0] != k_i) {
                same_index1 = false;
                continue; // go to "while(same_index1)" which will end this loop
            }
            
            int new_k = intra_pairs[intra_counter*2 + 1];
            
            // if we are accounting for defects, we check if the other half of the pair has been removed
            if (solver_type == 3 || solver_type == 4){
                if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
                    skip_here2 = 1;
                }
            }
            
            // we save this pair into our sparse matrix format
            if (skip_here1 == 0 && skip_here2 == 0){
                
                int sc_i = 0;
                int sc_j = 0;
                double new_pos_shift_x = 0.0;
                double new_pos_shift_y = 0.0;
                
                // Done
                if (boundary_condition == 1){
                    
                    int sheet_here = index_to_grid[k_i*4 + 3];
                    std::vector< std::vector<double> > supercell_orig = sdata[sheet_here].supercell;
                    double theta = angles[sheet_here];
                    std::vector< std::vector<double> > supercell;
                    supercell.resize(2);
                    
                    for (int dim = 0; dim < 2; ++dim){
                        supercell[dim].resize(2);
                        
                        supercell[dim][0] = cos(theta)*supercell_orig[dim][0] - sin(theta)*supercell_orig[dim][1];
                        supercell[dim][1] = sin(theta)*supercell_orig[dim][0] + cos(theta)*supercell_orig[dim][1];
                    }
                    
                    
                    sc_i = intra_sc_vecs[intra_counter][0];
                    sc_j = intra_sc_vecs[intra_counter][1];
                    new_pos_shift_x = sc_i*supercell[0][0] + sc_j*supercell[1][0];
                    new_pos_shift_y = sc_i*supercell[0][1] + sc_j*supercell[1][1];
                }
                
                //printf("rank %d added intra_pair for index %d: [%d,%d] = %f \n", rank, k, intra_pairs[intra_counter*2 + 0], intra_pairs[intra_counter*2 + 1],v[input_counter]);
                
                double x1 = i2pos[k_i*3 + 0] + new_pos_shift_x;
                double y1 = i2pos[k_i*3 + 1] + new_pos_shift_y;
                double z1 = i2pos[k_i*3 + 2];
                double x2 = i2pos[new_k*3 + 0];
                double y2 = i2pos[new_k*3 + 1];
                double z2 = i2pos[new_k*3 + 2];
                
                double t;
                double raw_t;
                
                if (strain_type == 1){
                    
                    int i0 = index_to_grid[k_i*4 + 0];
                    int j0 = index_to_grid[k_i*4 + 1];
                    int l0 = index_to_grid[k_i*4 + 2];
                    int s0 = index_to_grid[k_i*4 + 3];
                    
                    int ih = index_to_grid[new_k*4 + 0];
                    int jh = index_to_grid[new_k*4 + 1];
                    int lh = index_to_grid[new_k*4 + 2];
                    int sh = index_to_grid[new_k*4 + 3];
                    
                    double xh = i2pos[new_k*3 + 0];
                    double yh = i2pos[new_k*3 + 1];
                    double zh = i2pos[new_k*3 + 2];
                    
                    std::array<int,2> grid_disp = {{
                        ih-i0,
                        jh-j0 }};
                    
                    // we correct the grid values by the supercell_stride when there are periodic BCs
                    if (boundary_condition == 1){
                        grid_disp[0] = grid_disp[0] - intra_sc_vecs[intra_counter][0]*sdata[s0].supercell_stride[0][0] - intra_sc_vecs[intra_counter][1]*sdata[s0].supercell_stride[1][0];
                        grid_disp[1] = grid_disp[1] - intra_sc_vecs[intra_counter][0]*sdata[s0].supercell_stride[0][1] - intra_sc_vecs[intra_counter][1]*sdata[s0].supercell_stride[1][1];
                    }
                    
                    // take strain as avg of both orbital strains
                    std::vector< std::vector<double> > strain_here;
                    strain_here.resize(2);
                    
                    for (int i = 0; i < 2; ++i){
                        strain_here[i].resize(2);
                        for (int j = 0; j < 2; ++j){
                            //strain_here[i][j] = (strain[k_i][i][j]  + strain[new_k][i][j])/2.0;
                            //strain_here[i][j] = (strain[k_i][i][j])/2.0;
                            strain_here[i][j] = 0.0; // for now turn of in-plane strain corrections!!
                            
                        }
                    }
                    
                    // we need to find the bonding angle relative the u_ij definitions:
                    std::array<double, 2> strain_dir;
                    
                    strain_dir[0] = xh - x1;
                    strain_dir[1] = yh - y1;
                    
                    double strain_dir_norm = sqrt(strain_dir[0]*strain_dir[0] + strain_dir[1]*strain_dir[1]);
                    
                    // angle (counter-clockwise) from a bonding direction of +x
                    double strain_theta = 0;
                    
                    if (strain_dir_norm != 0){
                        strain_dir[0] = strain_dir[0]/strain_dir_norm;
                        strain_dir[1] = strain_dir[1]/strain_dir_norm;
                        
                        if (strain_dir[0] == 0.0){
                            if (strain_dir[1] > 0){
                                strain_theta = PI_2;
                            } else {
                                strain_theta = -PI_2;
                            }
                        } else {
                            double ratio = strain_dir[1]/strain_dir[0];
                            strain_theta = atan(ratio);
                            if (strain_dir[0] < 0){
                                strain_theta = strain_theta + PI;
                            }
                        }
                    }
                    
                    // now rotate;
                    std::vector< std::vector<double> > strain_rot;
                    strain_rot.resize(2);
                    strain_rot[0].resize(2);
                    strain_rot[1].resize(2);
                    
                    strain_rot[0][0] =  strain_here[0][0]*cos(strain_theta)*cos(strain_theta) +
                    strain_here[1][1]*sin(strain_theta)*sin(strain_theta) +
                    strain_here[0][1]*sin(strain_theta)*cos(strain_theta);
                    strain_rot[1][1] =  strain_here[0][0]*sin(strain_theta)*sin(strain_theta) +
                    strain_here[1][1]*cos(strain_theta)*cos(strain_theta) +
                    -2*strain_here[0][1]*sin(strain_theta)*cos(strain_theta);
                    strain_rot[0][1] =  (strain_here[1][1] - strain_here[0][0])*sin(strain_theta)*cos(strain_theta) +
                    strain_here[0][1]*(cos(strain_theta)*cos(strain_theta) - sin(strain_theta)*sin(strain_theta));
                    strain_rot[1][0] = strain_rot[0][1];
                    
                    Materials::Mat mat = sdata[s0].mat;
                    
                    if (mat_from_file == 0){
                        raw_t = Materials::intralayer_term(l0, lh, grid_disp, strain_rot, mat)/energy_rescale;
                        //raw_t = intra_pairs_t[intra_counter];
                        
                    } else {
                        raw_t = ReadMat::intralayer_term(l0, lh, grid_disp, loadedMatData, s0)/energy_rescale;
                    }
                    
                } else {
                    
                    raw_t = intra_pairs_t[intra_counter];
                    
                }
                
                // if it is the diagonal element, we "shift" the matrix up or down in energy scale (to make sure the spectrum fits in [-1,1] for the Chebyshev method)
                // Also, if electric field is included (elecOn == 1) we add in an on-site offset due to this gate voltage.
                
                if (new_k == k_i){
                    if (elecOn == 1){
                        t = (raw_t + energy_shift + onSiteE(x1,y1,z1,E))/energy_rescale;
                    } else if (elecOn == 0)
                        t = (raw_t + energy_shift)/energy_rescale;
                    // Otherwise we enter the value just with rescaling
                }
                else {
                    t = raw_t/energy_rescale;
                }
                
                // First get Magnetic field phase
                double phase = peierlsPhase(x1, x2, y1, y2, B);
                
                // Then get k dot R phase from the supercell wrapping
                // R is the vector connecting orbital k_i to orbital k_new
                phase = phase + k_vec[0]*(-new_pos_shift_x) + k_vec[1]*(-new_pos_shift_y);
                
                //t_cpx = t_cpx + std::polar(t, phase);
                t_cpx = std::polar(t, phase);
                //printf("sc_vec = [%lf, %lf], phase = %lf\n",-new_pos_shift_x,-new_pos_shift_y,phase);
                if (t_cpx != std::complex<double>(0.0)){
                    //printf("intra coupling (%d, %d) [%lf, %lf, %lf] -> [%lf, %lf, %lf] (%lf, %lf) = %lf \n",k_i,new_k,x1,y1,z1,x2,y2,z2,new_pos_shift_x,new_pos_shift_y,t);
                    
                    // moved to index checking block, as we want column info saved even if t != 0!
                    /*
                     col_index[input_counter] = new_k - current_index_reduction[new_k];
                     col_index_dx[input_counter] =  new_k - current_index_reduction[new_k];
                     col_index_dy[input_counter] =  new_k - current_index_reduction[new_k];
                     */
                    
                    v_c[input_counter] += t_cpx;
                    v_c_dx[input_counter] += (x1 - x2)*t_cpx;
                    v_c_dy[input_counter] += (y1 - y2)*t_cpx;
                    
                    if (chiral_on == 1){
                        
                        col_index_J_B_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_B_y[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_T_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_T_y[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_TB_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_TB_y[input_counter] = new_k - current_index_reduction[new_k];
                        
                        
                        if (index_to_grid[k_i*4 + 3] == 0){
                            v_c_J_B_x[input_counter] += (x2 - x1)*t_cpx;
                            v_c_J_B_y[input_counter] += (y2 - y1)*t_cpx;
                        } else if (index_to_grid[k_i*4 + 3] == 1) {
                            v_c_J_T_x[input_counter] += (x2 - x1)*t_cpx;
                            v_c_J_T_y[input_counter] += (y2 - y1)*t_cpx;
                        }
                    }
                }
                // check if next pair is identical (possible with periodic wrapping), or if we are at last element, to decide whether to save or not
                if ((intra_counter+1) == max_intra_pairs || k_i != intra_pairs[(intra_counter+1)*2 + 0] || new_k != intra_pairs[(intra_counter+1)*2 + 1]){
                    col_index[input_counter] = new_k - current_index_reduction[new_k];
                    col_index_dx[input_counter] =  new_k - current_index_reduction[new_k];
                    col_index_dy[input_counter] =  new_k - current_index_reduction[new_k];
                    if (chiral_on == 1){
                        col_index_J_B_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_B_y[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_T_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_T_y[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_TB_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_TB_y[input_counter] = new_k - current_index_reduction[new_k];
                    }
                    ++input_counter;
                }
            }
            
            ++intra_counter;
            
        }
        
        // set t_cpx to zero again, shouldn't be necessary but just in case!
        t_cpx = 0;
        
        // While we are still at the correct index in our inter_pairs list:
        bool same_index2 = true;
        while(same_index2) {
            
            int skip_here2 = 0;
            
            // if the first index of inter_pairs changes, we stop
            if (inter_counter == max_inter_pairs || inter_pairs[inter_counter*2 + 0] != k_i) {
                same_index2 = false;
                continue; // go back up to "while(same_index2)" which will end this loop
            }
            
            int new_k = inter_pairs[inter_counter*2 + 1];
            
            // if we are accounting for defects, we check if the other half of the pair has been removed
            if (solver_type == 3 || solver_type == 4){
                if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
                    skip_here2 = 1;
                }
            }
            
            // we save this pair into our sparse matrix format
            if (skip_here1 == 0 && skip_here2 == 0){
                
                int sc_i = 0;
                int sc_j = 0;
                double new_pos_shift_x = 0.0;
                double new_pos_shift_y = 0.0;
                
                if (boundary_condition == 1){
                    
                    std::vector< std::vector<double> > supercell = opts.getDoubleMat("supercell");
                    
                    sc_i = inter_sc_vecs[inter_counter][0];
                    sc_j = inter_sc_vecs[inter_counter][1];
                    new_pos_shift_x = sc_i*supercell[0][0] + sc_j*supercell[1][0];
                    new_pos_shift_y = sc_i*supercell[0][1] + sc_j*supercell[1][1];
                }
                
                //printf("rank %d added inter_pair for index %d: [%d,%d] = %f \n", rank, k, inter_pairs[inter_counter*2 + 0], inter_pairs[inter_counter*2 + 1],v[input_counter]);
                
                // get the position of both orbitals
                double x1 = i2pos[k_i*3 + 0] + new_pos_shift_x;
                double y1 = i2pos[k_i*3 + 1] + new_pos_shift_y;
                double z1 = i2pos[k_i*3 + 2];
                double x2 = i2pos[new_k*3 + 0];
                double y2 = i2pos[new_k*3 + 1];
                double z2 = i2pos[new_k*3 + 2];
                
                // and the orbit tag in their respective unit-cell
                int orbit1 = index_to_grid[k_i*4 + 2];
                int orbit2 = index_to_grid[new_k*4 + 2];
                
                // and material information
                Materials::Mat mat1 = sdata[index_to_grid[k_i*4 + 3]].mat;
                Materials::Mat mat2 = sdata[index_to_grid[new_k*4 + 3]].mat;
                
                // and the angle of the sheet each orbital is on
                double theta1;
                double theta2;
                if (strain_type == 0){
                    theta1 = angles[index_to_grid[k_i*4 + 3]];
                    theta2 = angles[index_to_grid[new_k*4 + 3]];
                } else {
                    // with strain, we need to compute the local theta by nearest neighbour bonding
                    theta1 = getLocalTheta(k_i, index_to_grid, i2pos, max_index);
                    theta2 = getLocalTheta(new_k, index_to_grid, i2pos, max_index);
                    //printf("(%d, %d), theta1 = %lf, theta2 = %lf \n",k_i, new_k, theta1*(180.0/M_PI), theta2*(180.0/M_PI));
                }
                
                std::array<double, 3> disp = {{
                    x2 - x1,
                    y2 - y1,
                    z2 - z1}};
                
                double t = Materials::interlayer_term(orbit1, orbit2, disp, theta1, theta2, mat1, mat2)/energy_rescale;
                
                //double t = Materials::interlayer_term(orbit1, orbit2, disp, 0.0, 0.0, mat1, mat2)/energy_rescale;
                //if(abs(t) > 1e-4){
                //if ((orbit1 == 7 && orbit2 == 10) || (orbit1 == 10 && orbit2 == 7)){
                //printf("sheet %d to %d, orbit %d to %d, [%lf, %lf, %lf]: %lf \n",index_to_grid[k_i*4 + 3],index_to_grid[new_k*4 + 3],orbit1, orbit2, disp[0],disp[1],disp[2],t);
                //printf("%lf, %lf, %lf, %d, %d, %d, %d, %lf \n",disp[0],disp[1],disp[2],orbit1,orbit2,sc_i,sc_j,t);
                //}
                //if (t != 0 ){
                
                // First get Magnetic field phase
                double phase = peierlsPhase(x1, x2, y1, y2, B);
                
                // Then get k dot R phase from the supercell wrapping
                phase = phase + k_vec[0]*(-new_pos_shift_x) + k_vec[1]*(-new_pos_shift_y);
                //printf("[%d, %d]: phase = %lf, k_vec = [%lf, %lf], R = [%lf, %lf] \n",k_i,new_k,phase, k_vec[0],k_vec[1], -new_pos_shift_x, -new_pos_shift_y);
                //t_cpx = t_cpx + std::polar(t, phase);
                t_cpx = std::polar(t,phase);
                if (t_cpx != std::complex<double>(0.0)){
                    
                    // moved to index checking block, as we want column info saved even if t != 0!
                    /*
                     col_index[input_counter] = new_k - current_index_reduction[new_k];
                     col_index_dx[input_counter] = new_k - current_index_reduction[new_k];
                     col_index_dy[input_counter] = new_k - current_index_reduction[new_k];
                     */
                    //printf("[%d, %d] @ (%d, %d): v_c[%d] += %lf | ",k_i,new_k,k_i,col_index[input_counter],input_counter,t_cpx.real());
                    //printf("new_k = %d, next_k = %d \n",inter_pairs[(inter_counter+0)*2 + 1], inter_pairs[(inter_counter+1)*2 + 1]);
                    
                    v_c[input_counter] += t_cpx;
                    v_c_dx[input_counter] += (x1 - x2)*t_cpx; // (delta_x)*t
                    v_c_dy[input_counter] += (y1 - y2)*t_cpx; // (delta_y)*t
                    
                    if (chiral_on == 1){
                        
                        col_index_J_B_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_B_y[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_T_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_T_y[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_TB_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_TB_y[input_counter] = new_k - current_index_reduction[new_k];
                        
                        v_c_J_TB_x[input_counter] += (x2 - x1)*t_cpx;
                        v_c_J_TB_y[input_counter] += (y2 - y1)*t_cpx;
                    }
                    
                }
                
                // check if next pair is identical (possible with periodic wrapping), or if we are at last element, to decide whether to save or not
                if ( (inter_counter+1) == max_inter_pairs || k_i != inter_pairs[(inter_counter+1)*2 + 0] || new_k != inter_pairs[(inter_counter+1)*2 + 1]){
                    col_index[input_counter] = new_k - current_index_reduction[new_k];
                    col_index_dx[input_counter] = new_k - current_index_reduction[new_k];
                    col_index_dy[input_counter] = new_k - current_index_reduction[new_k];
                    if (chiral_on == 1){
                        col_index_J_B_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_B_y[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_T_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_T_y[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_TB_x[input_counter] = new_k - current_index_reduction[new_k];
                        col_index_J_TB_y[input_counter] = new_k - current_index_reduction[new_k];
                    }
                    ++input_counter;
                }
                
            }
            
            ++inter_counter;
            
        }
        
    }
    
    // Save the end point + 1 of the last row
    row_pointer[local_max_index] = input_counter;
    row_pointer_dx[local_max_index] = input_counter;
    row_pointer_dy[local_max_index] = input_counter;
    
    if (chiral_on == 1){
        row_pointer_J_B_x[local_max_index] = input_counter;
        row_pointer_J_B_y[local_max_index] = input_counter;
        row_pointer_J_T_x[local_max_index] = input_counter;
        row_pointer_J_T_y[local_max_index] = input_counter;
        row_pointer_J_TB_x[local_max_index] = input_counter;
        row_pointer_J_TB_y[local_max_index] = input_counter;
    }
    
    int num_targets = jobIn.getInt("num_targets");
    std::vector<int> target_list = jobIn.getIntVec("target_list");
    
    // Generates the "alpha" vectors for conductivity calculation via the kernel polynomial method
    if (diagonalize != 1 && observable_type == 1){
        for (int t = 0; t < num_targets; ++t){
            //printf("target = %d, local_max_index = %d \n", target_list[t], local_max_index);
            
            int target_here = target_list[t] - current_index_reduction[target_list[t]];
            //printf("target_here = %d \n",target_here);
            
            double* target_vec = new double[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                target_vec[i] = 0;
            }
            
            target_vec[target_here] = 1;
            
            double temp_vec_x[local_max_index];
            double temp_vec_y[local_max_index];
            
            for (int i = 0; i < local_max_index; ++i){
                temp_vec_x[i] = 0;
                temp_vec_y[i] = 0;
            }
            
            dxH.vectorMultiply(target_vec,temp_vec_x,1,0);
            dyH.vectorMultiply(target_vec,temp_vec_y,1,0);
            
            // want < 0 | dxH, not dxH | 0 >, so need to use: Transpose( < 0 | dxH ) = - dxH | 0 >
            for (int i = 0; i < local_max_index; ++ i){
                alpha_0_x_arr[t*local_max_index + i] = -temp_vec_x[i];
                alpha_0_y_arr[t*local_max_index + i] = -temp_vec_y[i];
            }
            
            delete[] target_vec;
        }
    }
    
    // ------------------------------
    // Following saves Matrix to file
    //
    // Should only be used for 1-job processes, otherwise they will overwrite each other!
    
    int matrix_save = opts.getInt("matrix_save");
    if (matrix_save > 0){
        
        std::ofstream outFile;
        const char* extension = "_matrix.dat";
        outFile.open ((job_name + extension).c_str());
        
        for(int i = 0; i < local_max_index; ++i){
            int start_index = row_pointer[i];
            int stop_index = row_pointer[i+1];
            for(int j = start_index; j < stop_index; ++j){
                outFile << col_index[j] + 1 << ", " << i + 1 << ", " << v_c[j].real() << ", " << v_c[j].imag() << "\n";
            }
        }
        
        outFile.close();
        
        if (matrix_save > 1){
            
            std::ofstream outFile2;
            const char* extension2 = "_dxH_matrix.dat";
            outFile2.open ((job_name + extension2).c_str());
            
            for(int i = 0; i < local_max_index; ++i){
                int start_index = row_pointer_dx[i];
                int stop_index = row_pointer_dx[i+1];
                for(int j = start_index; j < stop_index; ++j){
                    outFile2 << col_index_dx[j] + 1 << ", " << i + 1 << ", " << v_c_dx[j].real() << ", " << v_c_dx[j].imag() << "\n";
                }
            }
            
            outFile2.close();
            //
            
            // End Matrix Save
            // ---------------
        }
    }
    
    
    
    // ------------------------------
    // Following saves positions and pairings to file
    
    int matrix_pos_save = opts.getInt("matrix_pos_save");
    if (matrix_pos_save > 0){
        
        std::ofstream outFile3;
        const char* extension3 = "_pos.dat";
        outFile3.open ((job_name + extension3).c_str());
        outFile3 << fixed;
        outFile3.precision(12);
        
        for(int i = 0; i < local_max_index; ++i){
            
            double x = i2pos[i*3 + 0];
            double y = i2pos[i*3 + 1];
            double z = i2pos[i*3 + 2];
            int s =  index_to_grid[i*4 + 3];
            int o =  index_to_grid[i*4 + 2];
            outFile3 << x << "    " << y << "    " << z << "    " << s << "    " << o << "\n";
        }
        
        outFile3.close();
    }
    
    if (matrix_pos_save > 1){
        std::ofstream outFile4;
        const char* extension4 = "_intra_pos.dat";
        outFile4.open ((job_name + extension4).c_str());
        
        for(int i = 0; i < max_intra_pairs; ++i){
            outFile4 <<
            i2pos[intra_pairs[2*i + 0]*3 + 0] << ", " << i2pos[intra_pairs[2*i + 0]*3 + 1] << ", " << i2pos[intra_pairs[2*i + 0]*3 + 2] << ", " <<
            i2pos[intra_pairs[2*i + 1]*3 + 0] << ", " << i2pos[intra_pairs[2*i + 1]*3 + 1] << ", " << i2pos[intra_pairs[2*i + 1]*3 + 2] << ", " <<
            intra_pairs[2*i + 0] << ", " << intra_pairs[2*i + 1] <<"\n";
        }
        
        outFile4.close();
        // ---------------
        std::ofstream outFile5;
        const char* extension5 = "_inter_pos.dat";
        outFile5.open ((job_name + extension5).c_str());
        
        for(int i = 0; i < max_inter_pairs; ++i){
            outFile5 <<
            i2pos[inter_pairs[2*i + 0]*3 + 0] << ", " << i2pos[inter_pairs[2*i + 0]*3 + 1] << ", " << i2pos[inter_pairs[2*i + 0]*3 + 2] << ", " <<
            i2pos[inter_pairs[2*i + 1]*3 + 0] << ", " << i2pos[inter_pairs[2*i + 1]*3 + 1] << ", " << i2pos[inter_pairs[2*i + 1]*3 + 2] << ", " <<
            inter_pairs[2*i + 0] << ", " << inter_pairs[2*i + 1] << "\n";
        }
        
        outFile5.close();
    }
    
    // End Matrix Save
    // ---------------
    
}

void Locality::generateMomH(SpMatrix &H, Job_params jobIn, int* index_to_grid, double* i2pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t, std::vector<int> current_index_reduction, int local_max_index){
    
    int solver_type = jobIn.getInt("solver_type");
    int jobID = jobIn.getInt("jobID");
    int magOn = jobIn.getInt("magOn");
    int elecOn = jobIn.getInt("elecOn");
    double B = jobIn.getDouble("B");
    double E = jobIn.getDouble("E");
    double energy_rescale = jobIn.getDouble("energy_rescale");
    double energy_shift = jobIn.getDouble("energy_shift");
    
    int mat_from_file = jobIn.getInt("mat_from_file");
    
    int num_mom_groups = opts.getInt("num_mom_groups");
    std::vector< std::vector<int> > mom_groups = opts.getIntMat("mom_groups");
    
    std::vector< std::vector<double> > shifts = jobIn.getDoubleMat("shifts");
    
    // Here we assume every layer has the same k-shift
    double shift_x = shifts[0][0];
    double shift_y = shifts[0][1];
    
    // Indexes how many inter terms we have entered so far
    int inter_counter = 0;
    
    // Indexes how many intra terms we have entered so far
    int intra_counter = 0;
    
    // Total number of expected non-zero matrix elements
    int max_nnz = max_intra_pairs + max_inter_pairs;
    
    // Do monolayer bloch theory for each sheet:
    
    int num_sheets = jobIn.getInt("num_sheets");
    
    int L = 5; // L is the half-side length of a 2L+1 x 2L+1 search grid for the monolayer bloch terms
    
    std::vector<int> N_R_array; // num_mom_groups vector of the number of real-space positions in each group
    std::vector< std::vector< std::vector<double> > > R_array; // (num_mom_groups) x (N_R) x (3) array of real-space positions
    std::vector< std::vector< std::vector< std::vector<double> > > > bloch_t_array; // (num_mom_groups) x (N_R) x (num_orbitals) x (num_orbitals) of real-space t_ij coupling values
    std::vector< std::vector<int> > orb_sheet_array; // (num_mom_groups) x (num_orbitals) of sheet indices of orbitals
    std::vector< std::vector<int> > orb_index_array; // (num_mom_groups) x (num_orbitals) of local indices of orbitals
    
    // Start loop over groups to create bloch-theory monolayer bands
    // For now, we need each layer to never have an orbital outside the unit-cell!
    // Will eventually want to add some controls for when the fixed-shift methods
    // push the individual layer atoms out of the original unit-cell (need to make
    // sure all orbitals are in the unit-cell when doing Bloch theory!)
    
    for (int g = 0; g < num_mom_groups; ++g){
        
        std::vector<int> sheets_here = mom_groups[g];
        int num_sheets_here = sheets_here.size();
        int total_num_orbs = 0;
        std::vector<std::vector<double> > orb_pos;
        double theta = angles[sheets_here[0]]; // we assume the lowest sheet has the angle info
        
        for (int s = 0; s < num_sheets_here; ++s){
            
            Materials::Mat local_mat = sdata[sheets_here[s]].mat;
            total_num_orbs += Materials::n_orbitals(sdata[s].mat);
            
        }
        
        
        orb_pos.resize(total_num_orbs);
        
        // Keeps track of if this is an intra or interlayer coupling
        std::vector<int> orb_sheet_indices;
        std::vector<int> orb_indices;
        std::vector<Materials::Mat> orb_materials;
        
        orb_sheet_indices.resize(total_num_orbs);
        orb_indices.resize(total_num_orbs);
        orb_materials.resize(total_num_orbs);
        int total_orb_index = 0;
        
        // Now again for orb positions
        for (int s = 0; s < num_sheets_here; ++s){
            
            Materials::Mat local_mat = sdata[sheets_here[s]].mat;
            int local_num_orbs = Materials::n_orbitals(local_mat);
            
            for (int o = 0; o < local_num_orbs; ++o){
                
                orb_pos[o].resize(3);
                for (int d = 0; d < 3; ++d){
                    orb_pos[o][d] = Materials::orbital_pos(local_mat, o, d);
                    if (d == 3){
                        orb_pos[o][d] += heights[sheets_here[s]];
                    }
                }
                
                orb_sheet_indices[total_orb_index] = sheets_here[s];
                orb_indices[total_orb_index] = o;
                orb_materials[total_orb_index] = local_mat;
                total_orb_index++;
            }
        }
        
        orb_sheet_array.push_back(orb_sheet_indices);
        orb_index_array.push_back(orb_indices);
        
        // compute bloch terms, we assume all sheets in the same group have a single unit-cell
        
        int N_R = 0;
        
        // Array of real-space lattice positions R
        std::vector< std::vector<double> > temp_R_array;
        // Array of bloch hopping parameters t (i.e. 2x2 blocks for 2 band graphene model)
        std::vector< std::vector< std::vector<double> > > temp_bloch_t_array;
        
        std::vector<std::vector<double> > a = sdata[sheets_here[0]].a;
        
        for (int i = -L; i < L+1; ++i){
            for (int j = -L; j < L+1; ++j){
                
                
                // unrotated variables
                
                double temp_x = i*a[0][0] + j*a[1][0];
                double temp_y = i*a[0][1] + j*a[1][1];
                double temp_z = 0;
                
                // rotated variables
                double x = temp_x*cos(theta) - temp_y*sin(theta);
                double y = temp_x*sin(theta) + temp_y*cos(theta);
                double z = temp_z;
                
                std::array<int,2> grid_disp = {{ i, j }};
                
                std::vector< std::vector<double> > temp_o1_bloch_t_array;
                for (int o1 = 0; o1 < total_num_orbs; ++o1){
                    std::vector<double> temp_o2_bloch_t_array;
                    for (int o2 = 0; o2 < total_num_orbs; ++o2){
                        
                        // compute R(i,j) based on <a> and <theta>
                        // find monolayer t_ij(R(i,j),o1,o2)
                        // if t_ij != 0, add R and t_ij to respective vectors, and ++N_R;
                        double t;
                        
                        // same layer
                        if (orb_sheet_indices[o1] == orb_sheet_indices[o2]){
                            t = Materials::intralayer_term(orb_indices[o1], orb_indices[o2], grid_disp, orb_materials[o1]);
                        } else { // different layers
                            
                            int global_shifts_on = jobIn.getInt("global_shifts_on");
                            
                            std::vector<double> o1_offset;
                            std::vector<double> o2_offset;
                            o1_offset.resize(3);
                            o2_offset.resize(3);
                            
                            for (int d = 0; d < 3; ++d){
                                o1_offset[d] = Materials::orbital_pos(orb_materials[o1], orb_indices[o1], d);
                                o2_offset[d] = Materials::orbital_pos(orb_materials[o2], orb_indices[o2], d);
                                if (global_shifts_on == 1){
                                    std::vector< std::vector<double> > g_shifts = jobIn.getDoubleMat("global_shifts");
                                    int s1 = orb_sheet_indices[o1];
                                    int s2 = orb_sheet_indices[o2];
                                    if (d < 2){
                                        o1_offset[d] = o1_offset[d] + g_shifts[s1][0]*a[0][d] + g_shifts[s1][1]*a[1][d];
                                        o2_offset[d] = o2_offset[d] + g_shifts[s2][0]*a[0][d] + g_shifts[s2][1]*a[1][d];
                                    } else {
                                        // no z-update here
                                    }
                                }
                            }
                            
                            double dx = temp_x + o2_offset[0] - o1_offset[0];
                            double dy = temp_y + o2_offset[1] - o1_offset[1];
                            double dz = (o2_offset[2] + heights[orb_sheet_indices[o2]]) - ( o1_offset[2] + heights[orb_sheet_indices[o1]]);
                            
                            std::array<double,3> pos_disp = {{dx, dy, dz}};
                            t = Materials::interlayer_term(orb_indices[o1], orb_indices[o2], pos_disp, theta, theta, orb_materials[o1], orb_materials[o2]);
                        }
                        
                        temp_o2_bloch_t_array.push_back(t);
                        
                    } // end of o2 loop
                    temp_o1_bloch_t_array.push_back(temp_o2_bloch_t_array);
                } // end of o1 loop
                
                std::vector<double> temp_r_vec;
                temp_r_vec.push_back(x);
                temp_r_vec.push_back(y);
                temp_r_vec.push_back(z);
                
                temp_R_array.push_back(temp_r_vec);
                temp_bloch_t_array.push_back(temp_o1_bloch_t_array);
            } // end of j loop
        } // end of i loop
        
        R_array.push_back(temp_R_array);
        bloch_t_array.push_back(temp_bloch_t_array);
        N_R_array.push_back( (int) temp_R_array.size());
        
    } // end of g loop
    
    
    // Sparse matrix format is 2 arrays with length = nnz, and 1 array with length = max_index + 1
    
    // col_index tells us the col of element i
    // v tells us the value of element i
    // row_pointer tells us the start of row j (at element i = row_pointer[j]) and the end of row j (at element i = row_pointer[j] - 1)
    
    H.setup(max_nnz, local_max_index, local_max_index);
    
    
    int* col_index = H.allocColIndx();
    int* row_pointer = H.allocRowPtr();
    
    std::complex<double>* v_c;
    
    v_c = H.allocCpxVal();
    
    // Count the current element, i.e. "i = input_counter"
    int input_counter = 0;
    
    // We keep the vacancy methods here, although momentum space in general has no vacancies...
    
    // Loop through every orbital (i.e. rows of H)
    for (int k_i = 0; k_i < max_index; ++k_i){
        
        int skip_here1 = 0;
        
        if (solver_type == 3 || solver_type == 4){
            if (current_index_reduction[k_i] + 1 == current_index_reduction[k_i + 1]){
                skip_here1 = 1;
            }
        }
        
        int k = k_i - current_index_reduction[k_i];
        
        // Save starting point of row k
        row_pointer[k] = input_counter;
        
        // While we are still at the correct index in our intra_pairs list:
        bool same_index1 = true;
        while(same_index1) {
            
            int skip_here2 = 0;
            
            // if the first index of intra_pairs changes, we stop
            if (intra_pairs[intra_counter*2 + 0] != k_i) {
                same_index1 = false;
                continue; // go to "while(same_index1)" which will end this loop
            }
            
            int new_k = intra_pairs[intra_counter*2 + 1];
            
            
            
            // if we are accounting for defects, we check if the other half of the pair has been removed
            if (solver_type == 3 || solver_type == 4){
                if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
                    skip_here2 = 1;
                }
            }
            
            // we save this pair into our sparse matrix format
            if (skip_here1 == 0 && skip_here2 == 0){
                col_index[input_counter] = new_k - current_index_reduction[new_k];
                
                std::complex<double> t(0.0,0.0);
                
                // Need to do monolayer bloch theory here!
                
                // finding pairing between k_i and new_k
                // and the orbit tag in their respective unit-cell
                int orbit1 = index_to_grid[k_i*4 + 2];
                int orbit2 = index_to_grid[new_k*4 + 2];
                int s1 = index_to_grid[k_i*4 + 3];
                int s2 = index_to_grid[new_k*4 + 3];
                
                // Find which group we are in
                int group_idx = -1;
                
                for (int g = 0; g < num_mom_groups; ++g){
                    std::vector<int> sheets_here = mom_groups[g];
                    int num_sheets_here = sheets_here.size();
                    
                    for (int s = 0; s < num_sheets_here; ++s){
                        if (sheets_here[s] == s1){
                            group_idx = g;
                        }
                    }
                    
                }
                
                int N_R = N_R_array[group_idx];
                
                // Now find what global orb index we want for this specific orb and sheet choice
                
                std::vector<int> local_sheet_indices = orb_sheet_array[group_idx];
                std::vector<int> local_orb_indices = orb_index_array[group_idx];
                int local_num_orbs = local_orb_indices.size();
                
                int local_o1 = -1;
                int local_o2 = -1;
                
                for (int o_search = 0; o_search < local_num_orbs; ++o_search){
                    if (s1 == local_sheet_indices[o_search] && orbit1 == local_orb_indices[o_search]){
                        local_o1 = o_search;
                    }
                    if (s2 == local_sheet_indices[o_search] && orbit2 == local_orb_indices[o_search]){
                        local_o2 = o_search;
                    }
                }
                
                for (int n = 0; n < N_R; ++n){
                    
                    double q1 = i2pos[k_i*3 + 0];
                    double q2 = i2pos[k_i*3 + 1];
                    double q3 = i2pos[k_i*3 + 2];
                    
                    double r1 = R_array[group_idx][n][0];
                    double r2 = R_array[group_idx][n][1];
                    double r3 = R_array[group_idx][n][2];
                    
                    double phase = q1*r1 + q2*r2 + q3*r3;
                    
                    double t_here = bloch_t_array[group_idx][n][local_o1][local_o2];
                    
                    
                    t = t + std::polar(t_here,phase);
                }
                
                // if it is the diagonal element, we "shift" the matrix up or down in energy scale (to make sure the spectrum fits in [-1,1] for the Chebyshev method)
                // Also, if electric field is included (elecOn == 1) we add in an on-site offset due to this gate voltage.
                if (new_k == k_i){
                    if (elecOn == 1){
                        double q1 = i2pos[k_i*3 + 0];
                        double q2 = i2pos[k_i*3 + 1];
                        double q3 = i2pos[k_i*3 + 2];
                        t = t + energy_shift + onSiteE(q1,q2,q3,E);
                    } else if (elecOn == 0)
                        t = t + energy_shift;
                    // Otherwise we enter the value just with rescaling
                }
                
                v_c[input_counter] = t/energy_rescale;
                
                ++input_counter;
            }
            
            ++intra_counter;
            
        }
        
        // While we are still at the correct index in our inter_pairs list:
        bool same_index2 = true;
        while(same_index2) {
            
            int skip_here2 = 0;
            
            // if the first index of inter_pairs changes, we stop
            if (inter_pairs[inter_counter*2 + 0] != k_i) {
                same_index2 = false;
                continue; // go to "while(same_index2)" which will end this loop
            }
            
            int new_k = inter_pairs[inter_counter*2 + 1];
            
            // if we are accounting for defects, we check if the other half of the pair has been removed
            if (solver_type == 3 || solver_type == 4){
                if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
                    skip_here2 = 1;
                }
            }
            
            // we save this pair into our sparse matrix format
            if (skip_here1 == 0 && skip_here2 == 0){
                
                // get the index of the other orbital in this term
                col_index[input_counter] = new_k - current_index_reduction[new_k];
                
                // get the position of both orbitals
                double x1 = i2pos[k_i*3 + 0];
                double y1 = i2pos[k_i*3 + 1];
                double z1 = i2pos[k_i*3 + 2];
                double x2 = i2pos[new_k*3 + 0];
                double y2 = i2pos[new_k*3 + 1];
                double z2 = i2pos[new_k*3 + 2];
                
                // get the sheet indices
                int s1 = index_to_grid[k_i*4 + 3];
                int s2 = index_to_grid[new_k*4 + 3];
                
                // and the orbit tag in their respective unit-cell
                int orbit1 = index_to_grid[k_i*4 + 2];
                int orbit2 = index_to_grid[new_k*4 + 2];
                
                // Momentum-space Interlayer coupling is evaluated at k = k2 + k1 - q
                //printf("[x1,y1] = [%lf, %lf], [x2,y2] = [%lf, %lf] \n",x1,y1,x2,y2);
                double dx = x2 + x1 - shift_x;
                double dy = y2 + y1 - shift_y;
                
                // The graphene interlayer functions are always 2pi/3 symmetric
                // we will impose this by averaging over rotated momentum values
                /*
                 double rot_param = (sqrt(3.0)/2.0);
                 double rot_dx = (-0.5)*dx - (rot_param)*dy;
                 double rot_dy = (rot_param)*dx + (-0.5)*dy;
                 double rot2_dx = (-0.5)*rot_dx - (rot_param)*rot_dy;
                 double rot2_dy = (rot_param)*rot_dx + (-0.5)*rot_dy;
                 */
                
                std::complex<double> t;
                
                // FFT file is based on sheet1 -> sheet2, so we need to keep track of which sheet k_i and new_k are on (i.e. k_i < new_k -> k_i on 1st sheet, k_i > new_k -> k_i on 2nd sheet)
                // Not checking for this is the equivalent of twisting the layer "theta" for 1->2 coupling and "-theta" for 2->1 coupling, which makes H non-Hermitian.
                
                double orb_disp_x = fftw_inter.get_fft_orb_disps(s1,s2,orbit1,orbit2,0);
                double orb_disp_y = fftw_inter.get_fft_orb_disps(s1,s2,orbit1,orbit2,1);
                //printf("[dx, dy] = [%lf, %lf], orb_disp_x = %lf, orb_disp_y = %lf \n",dx,dy, orb_disp_x, orb_disp_y);
                
                std::complex<double> phase_here;
                phase_here = std::polar(1.0, -dx*orb_disp_x - dy*orb_disp_y);
                
                std::complex<double> t_prephase = std::complex<double>(0.0,0.0);
                
                // impose rotation and mirror symm on the FFT of the interlayer Coupling
                // our graphene model has 3-fold rotational symm and a mirror-plane symmetry
                int rot_max = 3;
                // Dont mirror, seems to break electron-hole assymmetry (i.e. angular dependence of interaction)
                int mirror_max = 1;
                for (int rot_index = 0; rot_index < rot_max; ++rot_index){
                    for (int mirror_index = 0; mirror_index < mirror_max; ++mirror_index){
                        
                        double new_dx, new_dy;
                        double rot_dx, rot_dy;
                        
                        double theta = (2.0*M_PI/3.0)*((double)rot_index);
                        
                        rot_dx = cos(theta)*dx - sin(theta)*dy;
                        rot_dy = sin(theta)*dx + cos(theta)*dy;
                        
                        // when mirroring, new_vec = -rot_vec + 2*v*(v dot rot_vec)
                        // where v is the vector [mirror_plane_x mirror_plane_y]
                        
                        if (mirror_index == 0){
                            new_dx = rot_dx;
                            new_dy = rot_dy;
                        } else {
                            double mirror_theta = M_PI/2.0 + (angles[s1] + angles[s2])/2.0;
                            double mirror_plane_x = cos(mirror_theta);
                            double mirror_plane_y = sin(mirror_theta);
                            double mirror_dot_prod = mirror_plane_x*rot_dx + mirror_plane_y*rot_dy;
                            new_dx = -rot_dx + 2*mirror_plane_x*mirror_dot_prod;
                            new_dy = -rot_dy + 2*mirror_plane_x*mirror_dot_prod;
                        }
                        
                        if (1){
                            //if (k_i <= new_k){
                            // factor of 1/6 for the 3*2=6 symmetric points
                            t_prephase += (1.0 / ((double)(rot_max*mirror_max)) )*
                            (std::complex<double>(fftw_inter.interp_fft(new_dx,new_dy,s1,s2,orbit1,orbit2,0),
                                                  fftw_inter.interp_fft(new_dx,new_dy,s1,s2,orbit1,orbit2,1))
                             );
                        } else if (k_i > new_k) {
                            
                            // factor of 1/6 for the 3*2=6 symmetric points
                            t_prephase += (1.0 / ((double)(rot_max*mirror_max)) )*
                            (std::complex<double>(fftw_inter.interp_fft(new_dx,new_dy,s2,s1,orbit2,orbit1,0),
                                                  -fftw_inter.interp_fft(new_dx,new_dy,s2,s1,orbit2,orbit1,1))
                             );
                            
                        }
                    }
                }
                
                //printf("phase_here = [%lf,%lf] \n",phase_here.real(), phase_here.imag());
                //printf("t_prephase = [%lf,%lf] \n",t_prephase.real(), t_prephase.imag());
                
                
                t = phase_here*t_prephase;
                //printf("[%d, %d] = [%lf, %lf] \n",k_i,new_k,t.real(),t.imag());
                
                /*
                 if (k_i <= new_k){
                 t = (1.0/3.0)*(  std::complex<double>(fftw_inter.interp_fft(dx,dy,orbit1,orbit2,0),fftw_inter.interp_fft(dx,dy,orbit1,orbit2,1))
                 +std::complex<double>(fftw_inter.interp_fft(rot_dx,rot_dy,orbit1,orbit2,0),fftw_inter.interp_fft(rot_dx,rot_dy,orbit1,orbit2,1))
                 +std::complex<double>(fftw_inter.interp_fft(rot2_dx,rot2_dy,orbit1,orbit2,0),fftw_inter.interp_fft(rot2_dx,rot2_dy,orbit1,orbit2,1))
                 );
                 } else if (k_i > new_k) {
                 t = (1.0/3.0)*(  std::complex<double>(fftw_inter.interp_fft(dx,dy,orbit2,orbit1,0),-fftw_inter.interp_fft(dx,dy,orbit2,orbit1,1))
                 +std::complex<double>(fftw_inter.interp_fft(rot_dx,rot_dy,orbit2,orbit1,0),-fftw_inter.interp_fft(rot_dx,rot_dy,orbit2,orbit1,1))
                 +std::complex<double>(fftw_inter.interp_fft(rot2_dx,rot2_dy,orbit2,orbit1,0),-fftw_inter.interp_fft(rot2_dx,rot2_dy,orbit2,orbit1,1))
                 );
                 //t = std::complex<double>(fftw_inter.interp_fft(dx,dy,orbit1,orbit2,0),fftw_inter.interp_fft(dx,dy,orbit1,orbit2,1));
                 }
                 */
                
                // following used for debugging specific elements of H
                
                /*
                 if ((new_k == 120  && k_i == 243) || (new_k == 243 && k_i == 120)){
                 printf("[%d, %d]: t = [%lf, %lf] \n",k_i, new_k, t.real(),t.imag());
                 printf("dx = %lf, dy = %lf, o1 = %d, o2 = %d \n",dx,dy,orbit1,orbit2);
                 
                 }
                 */
                
                /*
                 if ( (k_i == 555 && new_k == 931) ) {
                 printf("interpair: [%d, %d] \n",k_i, new_k);
                 printf("orbits = [%d, %d] \n",orbit1, orbit2);
                 printf("pos1 = [%lf, %lf, %lf] \n", x1, y1, z1);
                 printf("pos2 = [%lf, %lf, %lf] \n", x2, y2, z2);
                 printf("dr = [%lf, %lf] \n", dx, dy);
                 printf("abs(t) = %lf \n",std::abs(t));
                 double temp_verbose;
                 temp_verbose = fftw_inter.interp_fft_v(dx,dy,orbit1,orbit2,0);
                 temp_verbose = fftw_inter.interp_fft_v(dx,dy,orbit1,orbit2,1);
                 printf("----------------------------- \n");
                 }
                 */
                
                
                if (std::abs(t) != 0){
                    //if (0){ // swapping this with above line turns off interlayer coupling
                    
                    //printf("[%d,%d] = %lf + %lf i \n",k_i,new_k,t.real(),t.imag());
                    
                    // Controls Momentum Interlayer coupling (set to 0 to turn off interlayer interaction!)
                    /*
                     if (k_i <= new_k){
                     v_c[input_counter] = t/energy_rescale;
                     } else if (k_i > new_k) {
                     v_c[input_counter] = std::conj(t)/energy_rescale;
                     }
                     */
                    v_c[input_counter] = t/energy_rescale;
                    ++input_counter;
                }
            }
            ++inter_counter;
            
        }
    }
    
    // Save the end point + 1 of the last row
    row_pointer[local_max_index] = input_counter;
    
    // ------------------------------
    // Following saves Matrix to file
    //
    // Should only be used for 1-job processes, otherwise they will overwrite each other!
    
    int matrix_save = opts.getInt("matrix_save");
    if (matrix_save > 0){
        
        std::ofstream outFile;
        const char* extension = "_matrix.dat";
        outFile.open ((job_name + extension).c_str());
        
        for(int i = 0; i < local_max_index; ++i){
            int start_index = row_pointer[i];
            int stop_index = row_pointer[i+1];
            for(int j = start_index; j < stop_index; ++j){
                outFile << col_index[j] + 1 << ", " << i + 1 << ", " << v_c[j].real() << ", " << v_c[j].imag() << "\n";
            }
        }
        
        outFile.close();
        
        if (matrix_save > 1){
            
            //
            // End Matrix Save
            // ---------------
        }
    }
    
    
    
    // ------------------------------
    // Following saves positions and pairings to file
    
    int matrix_pos_save = opts.getInt("matrix_pos_save");
    if (matrix_pos_save == 1){
        
        std::ofstream outFile3;
        const char* extension3 = "_pos.dat";
        outFile3.open ((job_name + extension3).c_str());
        outFile3 << fixed;
        outFile3.precision(12);
        
        for(int i = 0; i < local_max_index; ++i){
            
            double x = i2pos[i*3 + 0];
            double y = i2pos[i*3 + 1];
            double z = i2pos[i*3 + 2];
            int s =  index_to_grid[i*4 + 3];
            int o =  index_to_grid[i*4 + 2];
            outFile3 << x << "    " << y << "    " << z << "    " << s << "    " << o << "\n";
        }
        
        outFile3.close();
        // ---------------
        std::ofstream outFile4;
        const char* extension4 = "_intra_pos.dat";
        outFile4.open ((job_name + extension4).c_str());
        
        for(int i = 0; i < max_intra_pairs; ++i){
            outFile4 <<
            i2pos[intra_pairs[2*i + 0]*3 + 0] << ", " << i2pos[intra_pairs[2*i + 0]*3 + 1] << ", " << i2pos[intra_pairs[2*i + 0]*3 + 2] << ", " <<
            i2pos[intra_pairs[2*i + 1]*3 + 0] << ", " << i2pos[intra_pairs[2*i + 1]*3 + 1] << ", " << i2pos[intra_pairs[2*i + 1]*3 + 2] << ", " <<
            intra_pairs[2*i + 0] << ", " << intra_pairs[2*i + 1] <<"\n";
        }
        
        outFile4.close();
        // ---------------
        std::ofstream outFile5;
        const char* extension5 = "_inter_pos.dat";
        outFile5.open ((job_name + extension5).c_str());
        
        for(int i = 0; i < max_inter_pairs; ++i){
            outFile5 <<
            i2pos[inter_pairs[2*i + 0]*3 + 0] << ", " << i2pos[inter_pairs[2*i + 0]*3 + 1] << ", " << i2pos[inter_pairs[2*i + 0]*3 + 2] << ", " <<
            i2pos[inter_pairs[2*i + 1]*3 + 0] << ", " << i2pos[inter_pairs[2*i + 1]*3 + 1] << ", " << i2pos[inter_pairs[2*i + 1]*3 + 2] << ", " <<
            inter_pairs[2*i + 0] << ", " << inter_pairs[2*i + 1] << "\n";
        }
        
        outFile5.close();
    }
    
    // End Matrix Save
    // ---------------
    
    
}

void Locality::computeDosKPM(std::vector< std::vector<double> > &cheb_coeffs, SpMatrix &H, Job_params jobIn, std::vector<int> current_index_reduction, int local_max_index){
    
    int magOn = jobIn.getInt("magOn");
    int poly_order = jobIn.getInt("poly_order");
    int solver_space = jobIn.getInt("solver_space");
    
    int num_targets = jobIn.getInt("num_targets");
    std::vector<int> target_list = jobIn.getIntVec("target_list");
    
    // 0 for Real, 1 for Complex
    int complex_matrix = H.getType();
    
    if (complex_matrix == 0){
        
        // Starting vector for Chebyshev method is a unit-vector at the target orbital
        
        for (int t_count = 0; t_count < num_targets; ++t_count){
            
            int target_index = target_list[t_count] - current_index_reduction[target_list[t_count]];
            
            double T_prev[local_max_index];
            
            for (int i = 0; i < local_max_index; ++i){
                T_prev[i] = 0.0;
            }
            
            T_prev[target_index] = 1.0;
            
            // Temporary vector for algorithm ("current" vector T_j)
            double T_j[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                T_j[i] = 0.0;
            }
            
            H.vectorMultiply(T_prev, T_j, 1, 0);
            
            // Temporary vector for algorithm ("next" vector T_j+1)
            double T_next[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                T_next[i] = 0.0;
            }
            
            // first T value is always 1
            cheb_coeffs[t_count][0] = 1;
            
            // Next one is calculated simply
            cheb_coeffs[t_count][1] = T_j[target_index];
            
            // Now loop algorithm up to poly_order to find all T values
            // double alpha2 = 2;
            
            // want to do: T_next = 2*H*T_j - T_prev;
            H.vectorMultiply(T_j, T_next, 2, 0);
            for (int c = 0; c < local_max_index; ++c){
                T_next[c] = T_next[c] - T_prev[c];
            }
            
            for (int j = 2; j < poly_order/2; ++j){
                
                // reassign values from previous iteration
                for (int c = 0; c < local_max_index; ++c){
                    T_prev[c] = T_j[c];    //T_prev = T_j;
                    T_j[c] = T_next[c];    //T_j = T_next;
                }
                
                // get the jth entry
                cheb_coeffs[t_count][j] = T_j[target_index];
                
                // compute the next vector
                H.vectorMultiply(T_j, T_next, 2, 0);
                for (int c = 0; c < local_max_index; ++c){
                    T_next[c] = T_next[c] - T_prev[c];
                }
                
                // use Chebyshev recursion relations to populate the {2j,2j+1} entries.
                if (j >= poly_order/4){
                    
                    double an_an = 0;
                    double anp_an = 0;
                    for (int c = 0; c < local_max_index; ++c){
                        an_an += T_j[c]*T_j[c];
                        anp_an += T_next[c]*T_j[c];
                    }
                    
                    // u_{2n}     = 2*<a_n|a_n>        - u_0;
                    // u_{2n+1} = 2*<a_{n+1}|a_n>     - u_1;
                    cheb_coeffs[t_count][2*j]     = 2*an_an  - cheb_coeffs[t_count][0];
                    cheb_coeffs[t_count][2*j + 1] = 2*anp_an - cheb_coeffs[t_count][1];
                    
                }
                
                
                // print every 100 steps on print rank
                //if (rank == print_rank)
                //if (j%100 == 0)
                //printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
            }
            
        }
        
    }    // end real matrix block
    else if (complex_matrix == 1) {
        // Starting vector for Chebyshev method is a unit-vector at the target rbital
        for (int t_count = 0; t_count < num_targets; ++t_count){
            
            int target_index = target_list[t_count] - current_index_reduction[target_list[t_count]];
            
            std::complex<double>* T_prev = new std::complex<double>[local_max_index];
            
            for (int i = 0; i < local_max_index; ++i){
                T_prev[i] = 0.0;
            }
            
            T_prev[target_index] = 1.0;
            
            // Temporary vector for algorithm ("current" vector T_j)
            std::complex<double>* T_j = new std::complex<double>[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                T_j[i] = 0.0;
            }
            
            H.vectorMultiply(T_prev, T_j, 1, 0);
            
            // Temporary vector for algorithm ("next" vector T_j+1)
            std::complex<double>* T_next = new std::complex<double>[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                T_next[i] = 0.0;
            }
            
            // first T value is always 1
            cheb_coeffs[t_count][0] = 1;
            
            // Next one is calculated simply
            cheb_coeffs[t_count][1] = T_j[target_index].real();
            
            // Now loop algorithm up to poly_order to find all T values
            // double alpha2 = 2;
            
            // want to do: T_next = 2*H*T_j - T_prev;
            H.vectorMultiply(T_j, T_next, 2, 0);
            for (int c = 0; c < local_max_index; ++c){
                T_next[c] = T_next[c] - T_prev[c];
            }
            
            for (int j = 2; j < poly_order/2; ++j){
                
                // reassign values from previous iteration
                for (int c = 0; c < local_max_index; ++c){
                    T_prev[c] = T_j[c];    //T_prev = T_j;
                    T_j[c] = T_next[c];    //T_j = T_next;
                }
                
                // get the jth entry
                cheb_coeffs[t_count][j] = T_j[target_index].real();
                
                // compute the next vector
                H.vectorMultiply(T_j, T_next, 2, 0);
                for (int c = 0; c < local_max_index; ++c){
                    T_next[c] = T_next[c] - T_prev[c];
                }
                
                // use Chebyshev recursion relations to populate the {2j,2j+1} entries.
                if (j >= poly_order/4){
                    
                    double an_an = 0;
                    double anp_an = 0;
                    for (int c = 0; c < local_max_index; ++c){
                        an_an += (std::conj(T_j[c])*T_j[c]).real();
                        anp_an += (std::conj(T_next[c])*T_j[c]).real();
                    }
                    
                    // u_{2n}     = 2*<a_n|a_n>        - u_0;
                    // u_{2n+1} = 2*<a_{n+1}|a_n>     - u_1;
                    cheb_coeffs[t_count][2*j]     = 2*an_an  - cheb_coeffs[t_count][0];
                    cheb_coeffs[t_count][2*j + 1] = 2*anp_an - cheb_coeffs[t_count][1];
                    
                }
                
                // print every 100 steps on print rank
                //if (rank == print_rank)
                //if (j%100 == 0)
                //printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
            }
            
            delete[] T_prev;
            delete[] T_j;
            delete[] T_next;
            
        }
        
    } // end complex matrix block
}

void Locality::computeDosTraceKPM(std::vector< std::vector<double> > &cheb_coeffs, SpMatrix &H, Job_params jobIn, std::vector<int> current_index_reduction, int local_max_index){
    
    int magOn = jobIn.getInt("magOn");
    int poly_order = jobIn.getInt("poly_order");
    int solver_space = jobIn.getInt("solver_space");
    
    
    int num_targets = jobIn.getInt("num_targets");
    std::vector< std::vector<double> > kpm_trace_vectors = jobIn.getDoubleMat("kpm_trace_vectors");
    
    // 0 for Real, 1 for Complex
    int complex_matrix = H.getType();
    
    if (complex_matrix == 0){
        
        // Starting vector for Chebyshev method is a unit-vector at the target orbital
        
        for (int t_count = 0; t_count < num_targets; ++t_count){
            
            double T_prev[local_max_index];
            
            for (int i = 0; i < local_max_index; ++i){
                T_prev[i] = kpm_trace_vectors[t_count][i];
            }
            
            // Temporary vector for algorithm ("current" vector T_j)
            double T_j[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                T_j[i] = 0.0;
            }
            
            H.vectorMultiply(T_prev, T_j, 1, 0);
            
            // Temporary vector for algorithm ("next" vector T_j+1)
            double T_next[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                T_next[i] = 0.0;
            }
            
            // first T value
            double temp_sum = 0.0;
            for (int i = 0; i < local_max_index; ++i){
                temp_sum += kpm_trace_vectors[t_count][i]*kpm_trace_vectors[t_count][i];
            }
            cheb_coeffs[t_count][0] = temp_sum;
            
            // Next one
            temp_sum = 0.0;
            for (int i = 0; i < local_max_index; ++i){
                temp_sum += T_j[i]*kpm_trace_vectors[t_count][i];
            }
            cheb_coeffs[t_count][1] = temp_sum;
            
            // Now loop algorithm up to poly_order to find all T values
            // double alpha2 = 2;
            
            // want to do: T_next = 2*H*T_j - T_prev;
            H.vectorMultiply(T_j, T_next, 2, 0);
            for (int c = 0; c < local_max_index; ++c){
                T_next[c] = T_next[c] - T_prev[c];
            }
            
            for (int j = 2; j < poly_order/2; ++j){
                
                // reassign values from previous iteration
                for (int c = 0; c < local_max_index; ++c){
                    T_prev[c] = T_j[c];    //T_prev = T_j;
                    T_j[c] = T_next[c];    //T_j = T_next;
                }
                
                // get the jth entry
                temp_sum = 0.0;
                for (int i = 0; i < local_max_index; ++i){
                    temp_sum += T_j[i]*kpm_trace_vectors[t_count][i];
                }
                cheb_coeffs[t_count][j] = temp_sum;
                
                // compute the next vector
                H.vectorMultiply(T_j, T_next, 2, 0);
                for (int c = 0; c < local_max_index; ++c){
                    T_next[c] = T_next[c] - T_prev[c];
                }
                
                // use Chebyshev recursion relations to populate the {2j,2j+1} entries.
                if (j >= poly_order/4){
                    
                    double an_an = 0;
                    double anp_an = 0;
                    for (int c = 0; c < local_max_index; ++c){
                        an_an += T_j[c]*T_j[c];
                        anp_an += T_next[c]*T_j[c];
                    }
                    
                    // u_{2n}     = 2*<a_n|a_n>        - u_0;
                    // u_{2n+1} = 2*<a_{n+1}|a_n>     - u_1;
                    cheb_coeffs[t_count][2*j]     = 2*an_an  - cheb_coeffs[t_count][0];
                    cheb_coeffs[t_count][2*j + 1] = 2*anp_an - cheb_coeffs[t_count][1];
                    
                }
                
                
                // print every 100 steps on print rank
                //if (rank == print_rank)
                //if (j%100 == 0)
                //printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
            }
            
        }
        
    }    // end real matrix block
    else if (complex_matrix == 1) {
        // Starting vector for Chebyshev method is a unit-vector at the target rbital
        for (int t_count = 0; t_count < num_targets; ++t_count){
            
            std::complex<double>* T_prev = new std::complex<double>[local_max_index];
            
            for (int i = 0; i < local_max_index; ++i){
                T_prev[i] = kpm_trace_vectors[t_count][i];
            }
            
            // Temporary vector for algorithm ("current" vector T_j)
            std::complex<double>* T_j = new std::complex<double>[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                T_j[i] = 0.0;
            }
            
            H.vectorMultiply(T_prev, T_j, 1, 0);
            
            // Temporary vector for algorithm ("next" vector T_j+1)
            std::complex<double>* T_next = new std::complex<double>[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                T_next[i] = 0.0;
            }
            
            // first T value
            std::complex<double> temp_sum = 0.0;
            for (int i = 0; i < local_max_index; ++i){
                temp_sum += kpm_trace_vectors[t_count][i]*kpm_trace_vectors[t_count][i];
            }
            cheb_coeffs[t_count][0] = temp_sum.real();
            
            // Next one is calculated simply
            temp_sum = 0.0;
            for (int i = 0; i < local_max_index; ++i){
                temp_sum += T_j[i]*kpm_trace_vectors[t_count][i];
            }
            cheb_coeffs[t_count][1] = temp_sum.real();
            
            // Now loop algorithm up to poly_order to find all T values
            // double alpha2 = 2;
            
            // want to do: T_next = 2*H*T_j - T_prev;
            H.vectorMultiply(T_j, T_next, 2, 0);
            for (int c = 0; c < local_max_index; ++c){
                T_next[c] = T_next[c] - T_prev[c];
            }
            
            for (int j = 2; j < poly_order/2; ++j){
                
                // reassign values from previous iteration
                for (int c = 0; c < local_max_index; ++c){
                    T_prev[c] = T_j[c];    //T_prev = T_j;
                    T_j[c] = T_next[c];    //T_j = T_next;
                }
                
                // get the jth entry
                temp_sum = 0.0;
                for (int i = 0; i < local_max_index; ++i){
                    temp_sum += T_j[i]*kpm_trace_vectors[t_count][i];
                }
                cheb_coeffs[t_count][j] = temp_sum.real();
                
                // compute the next vector
                H.vectorMultiply(T_j, T_next, 2, 0);
                for (int c = 0; c < local_max_index; ++c){
                    T_next[c] = T_next[c] - T_prev[c];
                }
                
                // use Chebyshev recursion relations to populate the {2j,2j+1} entries.
                if (j >= poly_order/4){
                    
                    double an_an = 0;
                    double anp_an = 0;
                    for (int c = 0; c < local_max_index; ++c){
                        an_an += (std::conj(T_j[c])*T_j[c]).real();
                        anp_an += (std::conj(T_next[c])*T_j[c]).real();
                    }
                    
                    // u_{2n}     = 2*<a_n|a_n>        - u_0;
                    // u_{2n+1} = 2*<a_{n+1}|a_n>     - u_1;
                    cheb_coeffs[t_count][2*j]     = 2*an_an  - cheb_coeffs[t_count][0];
                    cheb_coeffs[t_count][2*j + 1] = 2*anp_an - cheb_coeffs[t_count][1];
                    
                }
                
                // print every 100 steps on print rank
                //if (rank == print_rank)
                //if (j%100 == 0)
                //printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
            }
            
            delete[] T_prev;
            delete[] T_j;
            delete[] T_next;
            
        }
        
    } // end complex matrix block
}

void Locality::computeCondKPM(std::vector< std::vector<double> > &cheb_coeffs, SpMatrix &H, SpMatrix &dxH, Job_params jobIn, std::vector<int> current_index_reduction, int local_max_index, double* alpha_0){
    
    // The current-current correlation can be calculated with Chebyshev polynomials
    // via a 2-Dimensional sampling:
    // for the C_nm corresponding to sigma_xx:
    // C_nm = <0| dxH * T_n * dxH * T_m |0>
    // We define:
    //                        <0| * dxH = <alpha_0| (computed during matrix construction)
    //                         <alpha_0| * T_n = <alpha|
    //                       T_m * |0> = |beta>
    // With these definitions, we get:
    //                         C_nm = <alpha|dxH|beta>
    //
    // To speed up calculations, we pre-compute beta for all relavant m
    // Then we loop over n to compute the relavant <alpha|
    // and save the matrix elements C_nm into cheb_coeffs.
    //
    // How much calculations are we actually doing? Assume the polynomial order = p
    //
    // Well, we need to compute each |beta> :     thats p sparse matrix-vector products (SpMV)
    // We also need to compute each <alpha|dxH: thats 2*p SpMV
    // Finally, we need to take their inner products: thats p^2 vector-vector products
    // So the total time seems to scale like:
    //                                                    3*p SpMV + p^2 vector-vector
    // since p SpMV is just the DOS time, the total time can be thought of as:
    //                                                     3*Time_DoS + p^2 vector-vector
    //
    // This is not too bad if the p^2 vector-vector products are on the same order of time as
    // the DOS time, and we know that the DOS calculation for that size has succedded in the past.
    // But in the case that the vec-vec calculations dominate we will see that this method becomes
    // much much slower than DOS calculations!
    //
    // Also, since we are probably interested in really large polynomial order (C_nm with n,m ~ 1000s)
    // we probably cannot do the most efficient calculation (compute all beta and save) becaues of memory limitations.
    //
    // Therefore we try to parallize over the vector-vector products:
    //                     the tag "cond_poly_par" allows for parallelization over beta.
    //
    // This allows for us to only compute a few columns of the large C_nm at a time, and use MPI to parallelize
    // over the sets of columns
    //
    // This is much faster at large p due to the number of vector-vector products but alows us to parallelize
    // in memory at least. Assuming we have we choose to paralellize the polynomial order by N:
    //                                          1) Calculate (at worst) 3*p SpMV products (alpha and beta)
    //                                         2) Calculate p*(p/N) vector-vector products
    // Now assume we have M workers, and that the SpMV takes time T_sp and V-V takes time T_vv, with S samples
    //                Old Serial Time:
    //                                S serial calculations of time 3*p SpMV + p^2 V-V         = S*(3*p*T_sp  + p^2 T_vv)
    //
    //                New Serial Time (where we round up S*N/M to an integer > 0):
    //                                (S*N)/M serial calculations of time 3*p SpMV + (p^2)/N V-V    = (S/M)*(N*3*p*T_sp + p^2 T_vv)
    //
    // So we see that the best choice is always a larger M, and the biggest allowed M is S*N.
    // If one does a test run to try and calculate T_sp and T_vv for the problem at hand, we can find an optimal
    // choice of N and M to minimize the total serial run time for a given p.
    //
    // A guess of how to pick N would be to note that if 3*p*T_sp < p^2 T_vv => T_sp < p/3 * T_vv
    // If T_sp is much smaller than p/3 * T_vv, we can get free embarrasing parallelization  as long as
    // we keep N*T_sp much smaller than p/3 * T_vv. Since this relation is linear in p, we can always keep the
    // total compute time roughly equal across multiple runs as long as we scale the number of workers with p!
    //
    //
    // !!! THE FOLLOWING IS NOT IMPLEMENTED, BUT NOTES ON A DIFFERENT APPROACH FOR PARALELLIZATION !!!
    //
    // There is also a more involed way to parrellelize the calculation time and memory however:
    //
    // 1) Have one MPI rank compute all <alpha| in "chunks for a single target
    // 2) Have (maybe a different) MPI rank compute all dxH|beta> in "chunks" for a single target
    // 3) Then send out to pure compute-node ranks (i.e. ranks which do not need to know
    // anything about the physical system), or "calculators", a specific sub-block of C_nm = <alpha|dxH|beta>
    //
    // This would mean that compute time would scale roughly like p (as in the DoS case)
    // given that you are also scaling the number of cores like p
    //
    // For example:
    // A single worker-core does a p=p0 DoS-like problem in time T
    // now we want to do p=p0*N with M "calculator" cores:
    //
    // 1) 1 worker-core does a "chunk" of dxH|beta> up to p/N = p0 ~ time would take T (DoS scaling)
    // 2) the same worker-core does <alpha| up to p in N chunks    ~ time would take T*N (DoS scaling)
    // 3) as each <alpha| "chunk" finishes, the worker core sends out the <beta| "chunk" and
    //    the "alpha" chunk to a "calculator" core
    // 4) the worker-core continues this until all "calculator" cores are actively working,
    //    then it gets the next pair of chunks ready in memory and spools for the next finished calculator
    // 5) as a calculator-core finishes, it sends the C_nm data to the worker-core and waits for a new job
    //
    // So then the problem is to optimize the number of calculator cores to the size of the chunks
    // i.e. if setting up a chunk-job takes time T_chunk, we want to pick the number of calculators M
    // such that the time to complete a chunk-job of matrix-vector products on the calculator
    // takes T_chunk*M long.
    //
    // In principle this could be determined by timing by the user or by a small utility job in the code
    //
    // This is a HUGE amount of extra programming, and would require siginificant additions to the MPI part of the code
    // but in principle this option is available given someone has the time (and energy) to program and test it properly
    
    int magOn = jobIn.getInt("magOn");
    int poly_order = jobIn.getInt("poly_order");
    
    int num_targets = jobIn.getInt("num_targets");
    std::vector<int> target_list = jobIn.getIntVec("target_list");
    
    // 0 for Real, 1 for Complex
    int complex_matrix = H.getType();
    
    // p1 represents the poly indices that beta will run over
    // that is to say, we always do  the "full" C_nm matrix during the <alpha| loop
    
    int start_p1 = 0;
    int cond_poly_par_scale = 1;
    
    int cond_poly_par = jobIn.getInt("cond_poly_par");
    
    if (cond_poly_par == 1){
        cond_poly_par_scale = jobIn.getInt("cond_poly_par_scale");
        start_p1 = jobIn.getInt("cond_poly_par_start_p1");
    }
    
    if (complex_matrix == 0){
        
        // Starting vector for Chebyshev method is a unit-vector at the target orbital
        for (int t_count = 0; t_count < num_targets; ++t_count){
            
            double* beta_array;
            beta_array = new double[local_max_index*poly_order];
            
            int target_index = target_list[t_count] - current_index_reduction[target_list[t_count]];
            
            double T_prev[local_max_index];
            
            // Temporary vector for algoirthm ("previous" vector T_{j-1})
            for (int i = 0; i < local_max_index; ++i){
                T_prev[i] = 0.0;
            }
            
            // Here T_0 = T_prev
            T_prev[target_index] = 1.0;
            
            // Temporary vector for algorithm ("current" vector T_j)
            double T_j[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                T_j[i] = 0.0;
            }
            
            // Temporary vector for algorithm ("next" vector T_{j+1})
            double T_next[local_max_index];
            for (int i = 0; i < local_max_index; ++i){
                T_next[i] = 0.0;
            }
            
            // we need to save T_0 = T_prev into beta_array if we are starting from zero
            if (start_p1 == 0){
                for (int i = 0; i < local_max_index; ++i){
                    beta_array[i] = T_prev[i];
                }
            }
            
            // Now we compute T_1 = H*T_0 = H*T_prev
            
            H.vectorMultiply(T_prev, T_j, 1, 0);
            
            // we need to save T_1 into beta_array if we are starting from zero
            if (start_p1 == 0){
                
                for (int i = 0; i < local_max_index; ++i){
                    beta_array[local_max_index + i] = T_j[i];
                }
            }
            
            // we need to save T_1 into beta_array if we are starting from one
            // starting from one is stupid... but good to have this work correctly just in case!
            if (start_p1 == 1){
                
                for (int i = 0; i < local_max_index; ++i){
                    beta_array[i] = T_j[i];
                }
            }
            
            // right now:
            //                    T_prev[]     = T_0
            //                     T_j[]         = T_1
            //                    T_next[]  = 0
            
            // if we are not starting at j = 0 (or 1)
            // we need to do a loop of the Cheb. algorithm here    to get up to j = start_p1
            for (int skip_p1 = 2; skip_p1 < start_p1; ++skip_p1){
                
                // right now:
                //                T_prev[] = T_{ skip_p1 - 2 }
                //                T_j[]    = T_{ skip_p1 - 1 }
                
                // compute T_{skip_p1}:
                //                T_j = 2*H*T_{j-1} - T_{j-2}
                for (int c = 0; c < local_max_index; ++c){
                    T_next[c] = T_next[c] - T_prev[c];
                }
                
                // right now:
                //                T_prev[] = T_{ skip_p1 - 2 }
                //                T_j[]    = T_{ skip_p1 - 1 }
                //                T_next[] = T_{ skip_p1 }
                
                
                // now we set-up for the next loop
                for (int c = 0; c < local_max_index; ++c){
                    T_prev[c] = T_j[c];    //T_prev = T_j;
                    T_j[c] = T_next[c]; //T_j = T_next;
                }
                
                // right now:
                //                T_prev[] = T_{ skip_p1 - 1 }
                //                T_j[]    = T_{ skip_p1 }
                //                T_next[] = T_{ skip_p1 }
                
            }
            
            for (int j = 0; j < poly_order-1; ++j){
                
                // A quick reindexing of j in case we had start_p1 = 0 or 1
                if (j == 0){
                    // already on the T_2 value if start_p1 = 0
                    if (start_p1 == 0){
                        j = 2;
                    }
                    
                    // already on the T_1 value if start_p1 = 1
                    if (start_p1 == 1){
                        j = 1;
                    }
                }
                
                // right now:
                //                T_prev[] = T_{ j + start_p1 - 2 }
                //                T_j[]    = T_{ j + start_p1 - 1 }
                
                // compute T_{j + start_p1}:
                //                T_j = 2*H*T_{j-1} - T_{j-2}
                H.vectorMultiply(T_j, T_next, 2, 0);
                for (int c = 0; c < local_max_index; ++c){
                    T_next[c] = T_next[c] - T_prev[c];
                }
                
                // right now:
                //                T_prev[] = T_{ j + start_p1 - 2 }
                //                T_j[]    = T_{ j + start_p1 - 1 }
                //                T_next[] = T_{ j + start_p1 }
                
                // saves T_{j} into beta_array
                for (int i = 0; i < local_max_index; ++i){
                    beta_array[j*local_max_index + i] = T_next[i];
                }
                
                // now we set-up for the next loop and
                // reassign values from previous iteration:
                // before:
                //                T_prev[] = T_{ j + start_p1 -2 }
                //                T_j[]    = T_{ j + start_p1 -1 }
                //                T_next[] = T_{ j + start_p1 }
                // after:
                //                T_prev[] = T_{ j + start_p1 - 1 }
                //                T_j[]    = T_{ j + start_p1 }
                //                T_next[] = T_{ j + start_p1 }
                
                for (int c = 0; c < local_max_index; ++c){
                    T_prev[c] = T_j[c];    //T_prev = T_j;
                    T_j[c] = T_next[c];    //T_j = T_next;
                }
                
                // print every 100 steps on print rank
                //if (rank == print_rank)
                //if (j%100 == 0)
                //printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
                
            }    // end beta_array construction
            
            // start loop over <alpha|*T_n for C_nm construction
            for (int i = 0; i < local_max_index; ++i){
                T_prev[i] = alpha_0[t_count*local_max_index + i];
            }
            
            for (int p1 = 0; p1 < poly_order; ++p1){
                
                // temp_vec is = - <alpha|dxH ... because dxH.vectorMultiply defines dxH|alpha> and Transpose(dxH) = -dxH
                double temp_vec[local_max_index];
                dxH.vectorMultiply(T_prev, temp_vec,1,0);
                
                double temp_sum = 0;
                for (int i = 0; i < local_max_index; ++i){
                    temp_sum = temp_sum + temp_vec[i]*beta_array[p1*local_max_index + i];
                }
                
                cheb_coeffs[t_count][p1 + 0*poly_order*cond_poly_par_scale] = temp_sum;
                
            }
            
            // Temporary vector for algorithm ("current" vector T_j)
            for (int i = 0; i < local_max_index; ++i){
                T_j[i] = 0.0;
            }
            
            // This is OK because H is Hermitian, Transpose(H) = H, so <alpha|H can be computed with H|alpha>
            H.vectorMultiply(T_prev, T_j, 1, 0);
            for (int p1 = 0; p1 < poly_order; ++p1){
                
                // temp_vec is = - <alpha|dxH ... because dxH.vectorMultiply defines dxH|alpha> and Transpose(dxH) = -dxH
                double temp_vec[local_max_index];
                dxH.vectorMultiply(T_j, temp_vec,1,0);
                
                double temp_sum = 0;
                for (int i = 0; i < local_max_index; ++i){
                    temp_sum = temp_sum + temp_vec[i]*beta_array[p1*local_max_index + i];
                }
                
                cheb_coeffs[t_count][p1 + 1*poly_order*cond_poly_par_scale] = temp_sum;
                
            }
            
            // Temporary vector for algorithm ("next" vector T_j+1)
            for (int i = 0; i < local_max_index; ++i){
                T_next[i] = 0.0;
            }
            
            // Now loop algorithm up to poly_order to find all T values
            // double alpha2 = 2;
            
            // want to do: T_next = 2*H*T_j - T_prev;
            H.vectorMultiply(T_j, T_next, 2, 0);
            for (int c = 0; c < local_max_index; ++c){
                T_next[c] = T_next[c] - T_prev[c];
            }
            
            for (int j = 2; j < poly_order*cond_poly_par_scale; ++j){
                
                // reassign values from previous iteration
                for (int c = 0; c < local_max_index; ++c){
                    T_prev[c] = T_j[c];    //T_prev = T_j;
                    T_j[c] = T_next[c];    //T_j = T_next;
                }
                
                
                for (int p1 = 0; p1 < poly_order; ++p1){
                    
                    // temp_vec is = - <alpha|dxH ... because dxH.vectorMultiply defines dxH|alpha> and Transpose(dxH) = -dxH
                    double temp_vec[local_max_index];
                    dxH.vectorMultiply(T_j, temp_vec,1,0);
                    
                    // Now we do <alpha|dxH|beta> = <0|dxH T_n(H) dxH T_m(H)|0>
                    // where <alpha| = <0|dxH T_n(H)
                    // and   |beta > = T_m(H)|0>
                    double temp_sum = 0;
                    for (int i = 0; i < local_max_index; ++i){
                        temp_sum = temp_sum + temp_vec[i]*beta_array[p1*local_max_index + i];
                    }
                    cheb_coeffs[t_count][p1 + j*poly_order*cond_poly_par_scale] = temp_sum;
                    
                }
                
                // compute the next vector
                H.vectorMultiply(T_j, T_next, 2, 0);
                for (int c = 0; c < local_max_index; ++c){
                    T_next[c] = T_next[c] - T_prev[c];
                }
                
                // print every 100 steps on print rank
                //if (rank == print_rank)
                //if (j%100 == 0)
                //printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
            }
            
            delete beta_array;
            
        } // end of t_count for loop
        
        
        
    }    // end complex_matrix == 0 block
    
    
    
}

void Locality::computeEigen(std::vector<double> &eigvals, DMatrix &eigvecs, DMatrix &kpm_dos, DMatrix &M_xx, DMatrix &M_yy, DMatrix &M_xy, SpMatrix &H, SpMatrix &dxH, SpMatrix &dyH, Job_params jobIn, std::vector<int> current_index_reduction, int local_max_index){
    
    int d_weights = jobIn.getInt("d_weights");
    int d_vecs = jobIn.getInt("d_vecs");
    int d_cond = jobIn.getInt("d_cond");
    int d_kpm_dos = jobIn.getInt("d_kpm_dos");
    int p = jobIn.getInt("poly_order");
    
    int d_type = jobIn.getInt("d_type");
    int debugPrint = 0;
    
    // H.eigenSolve should always be used now
    if (false){
        /*
         Eigen::MatrixXd H_dense = Eigen::MatrixXd::Zero(local_max_index,local_max_index);
         H.denseConvert(H_dense);
         //printf("Running EigenSolver... \n");
         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H_dense,Eigen::EigenvaluesOnly);
         //printf("EigenSolver complete! \n");
         
         eigvecs.setup(local_max_index, local_max_index);
         
         
         Eigen::VectorXd::Map(&eigvals[0], local_max_index) = es.eigenvalues();
         */
    } else {
        if (d_type == 0){ // full solve
            H.eigenSolve(eigvals, eigvecs,'V','A', 0, local_max_index);
        } else if (d_type == 1){ // only eigenvalues
            H.eigenSolve(eigvals, eigvecs,'N','A', 0, local_max_index);
        } else if (d_type == 2){ // only solve for eigenvalues within 4 indices of center of matrix
            int center_index = (int) (local_max_index/2);
            H.eigenSolve(eigvals, eigvecs,'N','I', center_index-3, center_index+4);
        }
        
        if (debugPrint == 1)
            eigvecs.debugPrint();
        
        /*
         Eigen::MatrixXd H_dense = Eigen::MatrixXd::Zero(local_max_index,local_max_index);
         H.denseConvert(H_dense);
         //printf("Running EigenSolver... \n");
         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H_dense);
         //printf("EigenSolver complete! \n");
         
         eigvecs.setup(local_max_index, local_max_index);
         
         double* eigvecs_ptr;
         eigvecs_ptr = eigvecs.allocRealVal();
         
         
         Eigen::VectorXd::Map(&eigvals[0], local_max_index) = es.eigenvalues();
         Eigen::MatrixXd::Map(&eigvecs_ptr[0], local_max_index, local_max_index) = es.eigenvectors();
         */
        
        // now we compute <n|j|m>, the current correlation matrix
        // for all eigenvectors |n>,|m>, and for j_x and j_y (via dxH, dyH)
        
        // Probably need to rewrite as sparse-matrix x dense-matrix multiplication for speedup
        
        if (d_kpm_dos > 0){
            
            int N = local_max_index;
            DMatrix cheb_eigvals;
            cheb_eigvals.setup(N,p);
            double* v_cheb;
            v_cheb = cheb_eigvals.allocRealVal();
            
            // COLUMN MAJOR ORDER (for Fortran MKL subroutines)
            //
            // Compute Cheby. coeffs. T_m(lambda_i) as a matrix
            
            for (int i = 0; i < N; ++i){
                // m = 0
                v_cheb[0*N + i] = 1;
                v_cheb[1*N + i] = eigvals[i];
            }
            
            // Chebyshev iteration on eigenvalues
            for (int m = 2; m < p; ++m){
                // start cheb iteration p
                for (int i = 0; i < N; ++i){
                    v_cheb[m*N + i] = 2*eigvals[i]*v_cheb[(m-1)*N + i] - v_cheb[(m-2)*N + i];
                }
            }
            
            cheb_eigvals.matrixMultiply(kpm_dos,cheb_eigvals, 1.0/local_max_index, 0, 'T', 'N');
            
        }
        
        if (d_cond > 0){
            
            //int p = jobIn.getInt("poly_order");
            // Compute J_x, J_y
            DMatrix J_x;
            DMatrix temp_matrix_x;
            
            if (debugPrint == 1)
                dxH.debugPrint();
            
            dxH.denseMatrixMultiply(temp_matrix_x, eigvecs, 1.0, 0.0);
            
            if (debugPrint == 1)
                temp_matrix_x.debugPrint();
            
            eigvecs.matrixMultiply(J_x, temp_matrix_x, 1.0, 0.0,'T','N');
            
            if (debugPrint == 1)
                J_x.debugPrint();
            //printf("J_x construction complete \n")
            
            
            int N = local_max_index;
            DMatrix cheb_eigvals;
            cheb_eigvals.setup(N,p);
            double* v_cheb;
            v_cheb = cheb_eigvals.allocRealVal();
            
            // COLUMN MAJOR ORDER (for Fortran MKL subroutines)
            //
            // Compute Cheby. coeffs. T_m(lambda_i) as a matrix
            
            for (int i = 0; i < N; ++i){
                // m = 0
                v_cheb[0*N + i] = 1;
                v_cheb[1*N + i] = eigvals[i];
            }
            
            // Chebyshev iteration on eigenvalues
            for (int m = 2; m < p; ++m){
                // start cheb iteration p
                for (int i = 0; i < N; ++i){
                    v_cheb[m*N + i] = 2.0*eigvals[i]*v_cheb[(m-1)*N + i] - v_cheb[(m-2)*N + i];
                }
            }
            
            // Finish calculating the M_ij tensors
            
            // J_xx element-wise via vdMul in MKL
            DMatrix J_xx;
            J_x.eleMatrixMultiply(J_xx, J_x, 1.0, 0.0, 'N', 'N');
            
            if (debugPrint == 1)
                J_xx.debugPrint();
            
            DMatrix temp_mat_xx;
            J_xx.matrixMultiply(temp_mat_xx, cheb_eigvals, 1.0, 0.0);
            
            if (debugPrint == 1)
                temp_mat_xx.debugPrint();
            
            //cheb_eigvals.matrixMultiply(M_xx,temp_mat_xx, 1.0/local_max_index, 0, 'T', 'N');
            cheb_eigvals.matrixMultiply(M_xx,temp_mat_xx, 1.0/((double) local_max_index), 0.0, 'T', 'N');
            
            
            if (debugPrint == 1)
                M_xx.debugPrint();
            
            if (d_cond > 1){
                
                DMatrix J_y;
                DMatrix temp_matrix_y;
                dyH.denseMatrixMultiply(temp_matrix_y, eigvecs, 1.0, 0.0);
                eigvecs.matrixMultiply(J_y, temp_matrix_y, 1, 0,'T','N');
                //printf("J_y construction complete \n");
                
                // J_yy, J_xy, element-wise via vdMul in MKL
                DMatrix J_yy;
                J_y.eleMatrixMultiply(J_yy, J_y, 1.0, 0.0, 'N', 'N');
                DMatrix J_xy;
                J_x.eleMatrixMultiply(J_xy, J_y, 1.0, 0.0, 'N', 'N');
                
                //printf("J_ij construction complete \n");
                
                DMatrix temp_mat_yy;
                J_yy.matrixMultiply(temp_mat_yy, cheb_eigvals, 1, 0);
                cheb_eigvals.matrixMultiply(M_yy,temp_mat_yy, 1.0/((double) local_max_index), 0.0, 'T', 'N');
                
                DMatrix temp_mat_xy;
                J_xy.matrixMultiply(temp_mat_xy, cheb_eigvals, 1, 0);
                cheb_eigvals.matrixMultiply(M_xy,temp_mat_xy, 1.0/((double) local_max_index), 0.0, 'T', 'N');
                
            }
            
        }
    }
    
}

void Locality::computeEigenComplex(std::vector<double> &eigvals, DMatrix &eigvecs,  DMatrix &kpm_dos, DMatrix &M_xx, DMatrix &M_yy, DMatrix &M_xy, SpMatrix &H, SpMatrix &dxH, SpMatrix &dyH, Job_params jobIn, std::vector<int> current_index_reduction, int local_max_index){
    
    int d_weights = jobIn.getInt("d_weights");
    int d_vecs = jobIn.getInt("d_vecs");
    int d_cond = jobIn.getInt("d_cond");
    int d_kpm_dos = jobIn.getInt("d_kpm_dos");
    int chiral_on = jobIn.getInt("chiral_on");
    int p = jobIn.getInt("poly_order");
    int d_type = jobIn.getInt("d_type");
    
    int debugPrint = 0;
    
    // H.eigenSolve should always be used now...
    if (false){
        /*
         Eigen::MatrixXcd H_dense = Eigen::MatrixXcd::Zero(local_max_index,local_max_index);
         H.denseConvert(H_dense);
         //printf("Running EigenSolver... \n");
         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(H_dense,Eigen::EigenvaluesOnly);
         
         //printf("EigenSolver complete! \n");
         
         Eigen::VectorXd::Map(&eigvals[0], local_max_index) = es.eigenvalues();
         */
    } else {
        
        if (d_type == 0){ // full solve
            H.eigenSolve(eigvals, eigvecs,'V','A', 0, local_max_index);
        } else if (d_type == 1){ // only eigenvalues
            H.eigenSolve(eigvals, eigvecs,'N','A', 0, local_max_index);
        } else if (d_type == 2){ // only solve for eigenvalues within 4 indices of center of matrix
            int center_index = (int) (local_max_index/2);
            H.eigenSolve(eigvals, eigvecs,'N','I', center_index-3, center_index+4);
        }
        
        if (debugPrint == 1)
            eigvecs.debugPrint();
        
        /*
         Eigen::MatrixXcd H_dense = Eigen::MatrixXcd::Zero(local_max_index,local_max_index);
         H.denseConvert(H_dense);
         //printf("Running EigenSolver... \n");
         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(H_dense);
         
         //printf("EigenSolver complete! \n");
         
         eigvecs.setup(local_max_index, local_max_index, 1);
         
         std::complex<double>* eigvecs_ptr;
         eigvecs_ptr = eigvecs.allocCpxVal();
         
         
         Eigen::VectorXd::Map(&eigvals[0], local_max_index) = es.eigenvalues();
         Eigen::MatrixXcd::Map(&eigvecs_ptr[0], local_max_index, local_max_index) = es.eigenvectors();
         */
        
        
        // now we compute <n|j|m>, the current correlation matrix
        // for all eigenvectors |n>,|m>, and for j_x and j_y (via dxH, dyH)
        
        // Probably need to rewrite as sparse-matrix x dense-matrix multiplication for speedup
        
        if (d_kpm_dos > 0){
            int N = local_max_index;
            DMatrix cheb_eigvals;
            cheb_eigvals.setup(N,p);
            double* v_cheb;
            v_cheb = cheb_eigvals.allocRealVal();
            
            // COLUMN MAJOR ORDER (for Fortran MKL subroutines)
            //
            // Compute Cheby. coeffs. T_m(lambda_i) as a matrix
            
            for (int i = 0; i < N; ++i){
                // m = 0
                v_cheb[0*N + i] = 1;
                v_cheb[1*N + i] = eigvals[i];
            }
            
            // Chebyshev iteration on eigenvalues
            for (int m = 2; m < p; ++m){
                // start cheb iteration p
                for (int i = 0; i < N; ++i){
                    v_cheb[m*N + i] = 2*eigvals[i]*v_cheb[(m-1)*N + i] - v_cheb[(m-2)*N + i];
                }
            }
            cheb_eigvals.matrixMultiply(kpm_dos,cheb_eigvals, 1.0/local_max_index, 0, 'T', 'N');
        }
        
        if (d_cond > 0){
            
            //int p = jobIn.getInt("poly_order");
            // Compute J_x, J_y
            
            // J_x = eigvecs*dxH*eigvecs
            
            DMatrix J_x;
            DMatrix temp_matrix_x;
            std::complex<double> a_cpx(1.0,0.0);
            std::complex<double> b_cpx(0.0,0.0);
            
            if (debugPrint == 1)
                dxH.debugPrint();
            
            dxH.denseMatrixMultiply(temp_matrix_x, eigvecs, a_cpx, b_cpx);
            
            if (debugPrint == 1)
                temp_matrix_x.debugPrint();
            
            eigvecs.matrixMultiply(J_x, temp_matrix_x, a_cpx, b_cpx,'C','N');
            
            if (debugPrint == 1)
                J_x.debugPrint();
            
            //printf("J_x construction complete \n")
            
            
            int N = local_max_index;
            DMatrix cheb_eigvals;
            cheb_eigvals.setup(N,p,1);
            std::complex<double>* v_cheb;
            v_cheb = cheb_eigvals.allocCpxVal();
            
            // COLUMN MAJOR ORDER (for Fortran MKL subroutines)
            //
            // Compute Cheby. coeffs. T_m(lambda_i) as a matrix
            
            for (int i = 0; i < N; ++i){
                // m = 0
                v_cheb[0*N + i] = 1;
                v_cheb[1*N + i] = eigvals[i];
            }
            
            // Chebyshev iteration on eigenvalues
            for (int m = 2; m < p; ++m){
                // start cheb iteration p
                for (int i = 0; i < N; ++i){
                    v_cheb[m*N + i] = 2.0*eigvals[i]*v_cheb[(m-1)*N + i] - v_cheb[(m-2)*N + i];
                }
            }
            // Finish calculating the the M_ij tensors
            
            // J_xx element-wise via vdMul in MKL
            // J_xx = J_x.*J_x, that is to say: J_xx^{ij} = (J_x^{ij})^2
            DMatrix J_xx;
            J_xx.setup(local_max_index,local_max_index,1);
            J_x.eleMatrixMultiply(J_xx, J_x, 1, 0, 'C', 'N');
            
            if (debugPrint == 1)
                J_xx.debugPrint();
            
            // temp_mat_xx = cheb_eigvals*J_xx*cheb_eigvals
            DMatrix temp_mat_xx;
            J_xx.matrixMultiply(temp_mat_xx, cheb_eigvals, a_cpx, b_cpx);
            
            if (debugPrint == 1)
                temp_mat_xx.debugPrint();
            
            cheb_eigvals.matrixMultiply(M_xx,temp_mat_xx, a_cpx/((std::complex<double>)local_max_index), b_cpx, 'T', 'N');
            
            if (debugPrint == 1)
                M_xx.debugPrint();
            
            if (d_cond > 1){
                
                DMatrix J_y;
                DMatrix temp_matrix_y;
                dyH.denseMatrixMultiply(temp_matrix_y, eigvecs, a_cpx, b_cpx);
                eigvecs.matrixMultiply(J_y, temp_matrix_y, a_cpx, b_cpx,'T','N');
                //printf("J_y construction complete \n");
                
                // J_yy, J_xy, element-wise via vdMul in MKL
                DMatrix J_yy;
                J_y.eleMatrixMultiply(J_yy, J_y, 1, 0, 'C', 'N');
                DMatrix J_xy;
                J_x.eleMatrixMultiply(J_xy, J_y, 1, 0, 'C', 'N');
                
                //printf("J_ij construction complete \n");
                
                DMatrix temp_mat_yy;
                J_yy.matrixMultiply(temp_mat_yy, cheb_eigvals, a_cpx, b_cpx);
                cheb_eigvals.matrixMultiply(M_yy,temp_mat_yy, a_cpx/((std::complex<double>)local_max_index), b_cpx, 'T', 'N');
                
                DMatrix temp_mat_xy;
                J_xy.matrixMultiply(temp_mat_xy, cheb_eigvals, a_cpx, b_cpx);
                cheb_eigvals.matrixMultiply(M_xy,temp_mat_xy, a_cpx/(std::complex<double>)local_max_index, b_cpx, 'T', 'N');
                
            }
        }
    }
    
}

void Locality::computeMatrixResponse(Job_params jobIn, std::vector<double> eigvals, DMatrix& eigvecs, SpMatrix& matIn, DMatrix& matOut){
    
    // mat_1 = eigvecs*matIn*eigvecs
    
    DMatrix mat_1;
    DMatrix temp_mat_1;
    std::complex<double> a_cpx(1,0);
    std::complex<double> b_cpx(0,0);
    matIn.denseMatrixMultiply(temp_mat_1, eigvecs, a_cpx, b_cpx);
    eigvecs.matrixMultiply(mat_1, temp_mat_1, a_cpx, b_cpx,'T','N');
    
    
    int N = matIn.getNumRows();
    int p = jobIn.getInt("poly_order");
    
    DMatrix cheb_eigvals;
    cheb_eigvals.setup(N,p,1);
    std::complex<double>* v_cheb;
    v_cheb = cheb_eigvals.allocCpxVal();
    
    // COLUMN MAJOR ORDER (for Fortran MKL subroutines)
    //
    // Compute Cheby. coeffs. T_m(lambda_i) as a matrix
    
    for (int i = 0; i < N; ++i){
        // m = 0
        v_cheb[0*N + i] = 1;
        v_cheb[1*N + i] = eigvals[i];
    }
    
    // Chebyshev iteration on eigenvalues
    for (int m = 2; m < p; ++m){
        // start cheb iteration p
        for (int i = 0; i < N; ++i){
            v_cheb[m*N + i] = ((std::complex<double>)2.0)*eigvals[i]*v_cheb[(m-1)*N + i] - v_cheb[(m-2)*N + i];
        }
    }
    
    // matOut = cheb_eigvals*( eigvecs*matIn*eigvecs )*cheb_eigvals
    //          = cheb_eigvals*(         mat_1         )*cheb_eigvals
    DMatrix temp_mat_2;
    mat_1.matrixMultiply(temp_mat_2, cheb_eigvals, a_cpx, b_cpx);
    cheb_eigvals.matrixMultiply(matOut,temp_mat_2, a_cpx/(std::complex<double>)N, b_cpx, 'T', 'N');
    
}

double Locality::peierlsPhase(double x1, double x2, double y1, double y2, double B_in){
    double preFactor = 3.14e-5; // This should be changed to be equal to (2 Pi)/(Flux Quanta), with Flux Quanta = h/e, in units of (T*Ang^2)^-1
    double phase = preFactor*B_in*(x2 - x1)*(0.5)*(y1 + y2);
    return phase;
}

double Locality::onSiteE(double x, double y, double z, double E_in){
    // E a constant in the Z direction, return z*E (physical prefactors missing!!!)
    return z*E_in;
}

std::vector< std::vector<double> > Locality::getReciprocal(int s){
    // calculate the monolayer reciprocal lattice vector of a given layer
    
    std::vector< std::vector<double> > orig_a = sdata[s].a;
    double theta = angles[s];
    
    std::vector< std::vector<double> > a;
    a.resize(2);
    a[0].resize(2);
    a[1].resize(2);
    
    
    for (int i = 0; i < 2; ++i){
        a[i][0] = orig_a[i][0]*cos(theta) - orig_a[i][1]*sin(theta);
        a[i][1] = orig_a[i][0]*sin(theta) + orig_a[i][1]*cos(theta);
    }
    
    std::vector< std::vector<double> > b;
    b = getReciprocal(a);
    return b;
}


std::vector< std::vector<double> > Locality::getReciprocal(std::vector< std::vector<double> > a_in){
    
    // Here we assume a is a 2x2 matrix
    std::vector< std::vector<double> > a;
    a.resize(3);
    a[0].resize(3);
    a[0][0] = a_in[0][0];
    a[0][1] = a_in[0][1];
    a[0][2] = 0.0;
    a[1].resize(3);
    a[1][0] = a_in[1][0];
    a[1][1] = a_in[1][1];
    a[1][2] = 0.0;
    a[2].resize(3);
    a[2][0] = 0.0;
    a[2][1] = 0.0;
    a[2][2] = 1.0;
    
    std::vector< std::vector<double> > b;
    b.resize(2);
    for (int i = 0; i < 2; ++i){
        b[i].resize(2);
    }
    
    double denom1 = 0.0;
    double denom2 = 0.0;
    //double denom3 = 0.0;
    
    for (int j = 0; j < 3; ++j) {
        denom1 += a[0][j]*crossProd(a[1],a[2],j);
        denom2 += a[1][j]*crossProd(a[2],a[0],j);
        //denom3 += a[2][j]*crossProd(a[0],a[1],j);
    }
    
    for (int k = 0; k < 2; ++k){
        b[0][k] = 2*M_PI*crossProd(a[1],a[2],k)/denom1;
        b[1][k] = 2*M_PI*crossProd(a[2],a[0],k)/denom2;
        //b[2][k] = 2*M_PI*crossProd(a[0],a[1],k)/denom3;
    }
    
    /*
     for (int l = 0; l < 3; ++l){
     printf("b[%d] = [%lf, %lf, %lf] \n", l, b[l][0], b[l][1], b[l][2]);
     }
     */
    
    return b;
}

double Locality::crossProd(std::vector<double> x, std::vector<double> y, int dim){
    
    if (dim == 0){
        return ( x[1]*y[2] - x[2]*y[1] );
    } else if (dim == 1) {
        return ( x[2]*y[0] - x[0]*y[2] );
    } else if (dim == 2) {
        return ( x[0]*y[1] - x[1]*y[0] );
    } else {
        return 0;
    }
    
}

/*
 void Locality::writeBufferToFile(double* buff, int length, std::string file){
 // Writes a buffer of doubles to a file in binary
 
 // Need ios library, and some other bugs left to be sorted out
 
 std::ofstream fout(file.c_str(), std::ios::binary);
 fout.write(reinterpret_cast<char*>(buff), std::streamsize(sizeof(double)*length));
 fout.close();
 
 }
 */

void Locality::save(){
    
    if (rank != root){
        int jobCount = solverTimes.size() / 2;
        double avgTime = 0;
        double maxTime = 0;
        double minTime = difftime(solverTimes[1],solverTimes[0]);
        
        for (int i = 0; i < jobCount; ++i){
            double tempTime = difftime(solverTimes[2*i + 1],solverTimes[2*i]);
            if (tempTime > maxTime)
                maxTime = tempTime;
            if (tempTime < minTime)
                minTime = tempTime;
            avgTime += tempTime/( (double) jobCount);
        }
        
        printf(    "=======- RANK %d TIMING -======= \n"
               "Number of Jobs: %d \n"
               "Avg Solve Time: %lf sec \n"
               "Max Solve Time: %lf sec \n"
               "Min Solve Time: %lf sec \n"
               "================================ \n",
               rank, jobCount, avgTime, maxTime, minTime);
    }
    
    if (rank == print_rank){
        printf(    "=======- TIMING RESULTS -======= \n"
               "Construct Geom: %lf sec \n"
               "Solver        : %lf sec \n"
               "================================ \n",
               difftime(constructEnd, constructStart),
               difftime(solveEnd, solveStart));
    }
}

void Locality::printTiming(std::vector< Job_params > results){
    
    std::vector<std::vector<double> > times;
    std::vector<std::string> tags = results[0].getStringVec("cpu_time_type");
    
    int num_r = results.size();
    int num_t = tags.size();
    
    // should eventually add some catch/throws to make sure every result has same # of times and types
    
    for (int r = 0; r < num_r; ++r){
        std::vector<double> temp_times;
        times.push_back(results[r].getDoubleVec("cpu_time"));
    }
    
    std::vector<double> avg_times;
    avg_times.resize(num_t);
    
    for (int r = 0; r < num_r; ++r){
        for (int t = 0; t < num_t; ++t){
            double temp_diff = times[r][t];
            avg_times[t] = avg_times[t] + temp_diff/((double) num_r);
        }
    }
    printf("=========- AVG TIMINGS -========\n");
    for (int t = 0; t < num_t; ++t){
        int text_length = tags[t].length();
        printf("%s",tags[t].c_str());
        if (text_length < 14) {
            for (int s = 0; s < 14 - text_length; ++s){
                printf(" ");
            }
        }
        printf(": %lf sec \n",avg_times[t]);
    }
    printf("================================\n");
    
}

void Locality::finMPI(){
    //printf("rank %d finalizing MPI. \n", rank);
    
    MPI_Finalize();
    
}
