
PetscErrorCode ierr;

PetscInt N; // matrix size is N x N
Mat A;
PetscInt nnz[N]; // nnz[i] is the number of column entries in row i

char help[] = "what this program does in brief can go here."

PetscInitialize(&argc,&args,(char*)0,help);

ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
ierr = MatSetType(A,MATSEQAIJ);CHKERRQ(ierr);
ierr = MatSetSizes(A,N,N,N,N);CHKERRQ(ierr);

ierr = MatSeqAIJSetPreallocation(A,NULL,nnz);


// ******

// some loop here . . .

PetscInt m; // number of rows being added
PetscInt idxm[m]; // row index values
PetscInt n; // number of cols being added
PetscInt idxn[n]; // col index values
PetscScalar v[n]; // entry values

// typically you want m = 1, because v corresponds to the columns
// from what I understand. So if you want the same row on several rows
// you can use m > 1, but otherwise just m = 1

ierr = MatSetValues(A, m, idxm, n, idxn, v, INSERT_VALUES); 
// loop ends

// *******


ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
ierr = PetscLogStagePop();CHKERRQ(ierr);


ierr = MatDestroy(&A);CHKERRQ(ierr);


ierr = PetscFinalize();

