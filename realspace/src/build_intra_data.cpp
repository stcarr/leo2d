#include "intralayer_coupling.h"
#include <iostream>
#include <fstream>
#include <string>

int main()
{
std::ofstream fout("tmdc_mose2");
int L = 3;
int orbitals = 11;

double A[2][2] = { { 1, -.5}, {0, .86602540378}};
double alpha = 3.32;
double basis[11][3] = 
{
{ 0, 0, 0},
{ 0, 0, 0},
{ 0, 0, 0},
{ 0, 0, 0},
{ 0, 0, 0},
{ 0.5, 0.28867513459, -0.50451807228},
{ 0.5, 0.28867513459, -0.50451807228},
{ 0.5, 0.28867513459, -0.50451807228},
{ 0.5, 0.28867513459, 0.50451807228},
{ 0.5, 0.28867513459, 0.50451807228},
{ 0.5, 0.28867513459, 0.50451807228} };

// z is fixed to be 1

fout << L << ", " << orbitals << ", ";
for (int i = 0; i < orbitals*orbitals-2; i++)
fout << "0, ";
fout << "\n";
for (int i = 0; i < 2*L+1; i++)
	for (int j =  0; j < 2*L+1; j++)
{
		for (int k = 0; k < orbitals; k++)
			for (int l = 0; l < orbitals; l++)
			{
				double r[3];
				double n1 = i-L;
				double n2 = j-L;
				r[0] = A[0][0]*n1 + A[0][1]*n2 + basis[k][0];
				r[1] = A[1][0]*n1 + A[1][1]*n2 + basis[k][1];
				r[2] = basis[k][2];
				for (int t = 0; t < 3; t++)
					r[t] = alpha*r[t];
				fout << intralayer_term(r[0],r[1],r[2],alpha*basis[l][0],alpha*basis[l][1],alpha*basis[l][2],k,l,3);
				fout << ", ";

			}
fout << std::endl;
}
fout.close();




/*
For MoS2:

ALPHA = 3.18
UNITCELL1 = 1 0 0
UNITCELL2 = -0.5 0.86602540378 0
UNITCELL3 = 0 0 1
ORBITAL POSITIONS:
1: 0 0 0
2: 0 0 0
3: 0 0 0
4: 0 0 0
5: 0 0 0
6: 0.5 0.28867513459 -0.50451807228
7: 0.5 0.28867513459 -0.50451807228
8: 0.5 0.28867513459 -0.50451807228
9: 0.5 0.28867513459 0.50451807228
10: 0.5 0.28867513459 0.50451807228
11: 0.5 0.28867513459 0.50451807228

for WSe2 it is the same but ALPHA = 3.32
*/







return 0;
}

