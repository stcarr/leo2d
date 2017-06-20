/* 
 * File:   job_params.h
 * Author: Stephen
 * 
 * Created on August 24, 2016, 2:46 PM
 */

#ifndef JOB_PARAMS_H
#define JOB_PARAMS_H

#include <string>
#include <vector>

using namespace std;

class Job_params {
    private:
		
		vector<int> int_params;
		vector<string> int_param_tags;
		
		vector<double> double_params;
		vector<string>  double_param_tags;
		
		vector<string> string_params;
		vector<string> string_param_tags;
		
		vector< vector<int> > int_vec_params;
		vector<string> int_vec_param_tags;

		vector< vector<double> > double_vec_params;
		vector<string> double_vec_param_tags;

		vector< vector<string> > string_vec_params;
		vector<string> string_vec_param_tags;

		vector< vector< vector<int> > > int_mat_params;
		vector<string> int_mat_param_tags;

		vector< vector< vector<double> > > double_mat_params;
		vector<string> double_mat_param_tags;			
	
	public:
        Job_params();
        Job_params(Job_params& orig);
        ~Job_params();
		
		void setParam(string, int);
		void setParam(string, double);
		void setParam(string, string);

		void setParam(string, vector<int>);
		void setParam(string, vector<double>);
		void setParam(string, vector<string>);
		
		void setParam(string, vector< vector<int> >);
		void setParam(string, vector< vector<double> >);
		
		int getInt(string);
		double getDouble(string);
		string getString(string);
		
		vector<int> getIntVec(string);
		vector<double> getDoubleVec(string);
		vector<string> getStringVec(string);
		
		vector< vector<int> > getIntMat(string);
		vector< vector<double> > getDoubleMat(string);
		
		vector<string> getParamTags(string);
		
		void printParams();
		
		int recvSpool();
		
		void sendParams(int);
		void recvParams(int);
		
		void sendTag(string, int);
		string recvTag(int);
		
		void sendInt(string,int,int);
		void recvInt(int);
		void sendDouble(string,double,int);
		void recvDouble(int);
		void sendString(string,string,int);
		void recvString(int);
		
		void sendIntVec(string,vector<int>,int);
		void recvIntVec(int);
		void sendDoubleVec(string,vector<double>,int);
		void recvDoubleVec(int);
		void sendStringVec(string,vector<string>,int);
		void recvStringVec(int);	

		void sendIntMat(string,vector< vector<int> >,int);
		void recvIntMat(int);
		void sendDoubleMat(string,vector< vector<double> >,int);
		void recvDoubleMat(int);	
		
};

#endif /* JOB_PARAMS_H */