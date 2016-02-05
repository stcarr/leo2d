#include <slepceps.h>



Mat A;
EPS eps;

char help[] = "some small description can go here";
SlepcInitialize(&argc,&argv,(char*)0,help);

EPSCreate(PETSC_COMM_WORLD,&eps);

EPSSetOperators(eps,A,NULL);

EPSSetProblemType(eps,EPS_HEP); // sets as hermitian



PetscReal E1,E2; // energy range of interest

EPSSetWhichEigenpairs(eps, EPS_ALL);
EPSSetInterval(eps,E1,E2);


EPSSetFromOptions(eps);

PetscInt nconv;

EPSGetConverged(eps,&nconv); // how many converged eigenpairs

Vec xr,xi; // real part and imaginary part of eigenvector
PetscScalar kr,ki; // real part, imaginary part of eigenvalue
PetscInt i; // which eigenpair

EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);




