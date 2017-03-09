#ifdef __NAG__
#include <nag.h>
#include <nag_stdlib.h>
#include <nag_string.h>
#include <nagf04.h>
#include <nagx04.h>
#include <nagf11.h>
#include <stdio.h>
#include "sim2d.h"

#define TDA 400
#define ARRGROW 2
#define NMAX 5
#define SCALE 1000.0
/**
 * Direct NAG solver only used for debuging
 *
 */
int
DNAG (double *ae, double *aw, double *an, double *as, double *ap,
        double *b, double *u, int ni, int nj, int maxit, int method)
{
    double errbnd, rcond;
    Integer exit_status, i, j, n, nrhs, pda, pdb;
    n = ni;
    nrhs = 1;
    // arrays
    double a[NMAX*NMAX][TDA]; 

    // Nag Types
    NagError fail;
    Nag_OrderType order;
    
#ifdef NAG_COLUMN_MAJOR
#define A(I,J) a[(J-1)*pda + I - 1]
#define B(I,J) b[(J-1)*pdb + I - 1]
    order = Nag_ColMajor;
#else
#define A(I,J) a[(I-1)*pda + J - 1]
#define B(I,J) b[(I-1)*pdb + J - 1]
    order = Nag_RowMajor;
#endif
    
    exit_status = 0;
    INIT_FAIL(fail);
    // allocate memory
#ifdef NAG_COLUMN_MAJOR
    pda = n;
    pdb = n;
#else
    pda = n;
    pdb = nrhs;
#endif
    for (i = 0; i < n*n; i++) {
        for (j = 0; j < n*n; j++) {
            a[i][j] = 0.;
        }
    }
    for (i = 0; i < n*n; i++) {
        for (j = 0; j < n*n; j++) {
        if (i == j)
            a[i][j] = ap[i]*SCALE;
        if (i == j-1)
            a[i][j] = aw[i]*SCALE;
        if (i == j+1)
            a[i][j] = ae[i]*SCALE;
        if (i-1 == j)
            a[i][j] = as[i]*SCALE;
        if (i+1 == j)
            a[i][j] = an[i]*SCALE;
        }
    }
    // Before you solve equations write out the 
    // matrix and vector AX=B
    Vprintf("f04arc Test Output\n");
    Vprintf("%ld", n);
    if (n < 1 || n > NMAX) {
        Vfprintf(stderr,"n is out of range: n = %5ld\n",n);
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < n; i++) { 
        for (j = 0; j < n; j++)
            Vprintf("%lf ", a[i][j]);
        Vprintf("\n");
    }
    for (i = 0; i < n; i++)
        Vprintf("%lf\n",b[i]);

    // Solve the equations AX = B for X 
    // Computes the solution and error-bound to a real system of
    // linear equations
    f04arc(n, &a[0][0], TDA, b, u, 
            NAGERR_DEFAULT);
END:
    if (a) NAG_FREE(a);
    if (b) NAG_FREE(b);
    return exit_status;
}

extern Csim2d *csim2d;
static int max_size;



/**
 * f11dcc solver from NAG library.
 * for more look in NAG documentation Mark 8.
 * This is not optimized due transforming of the formats.
 *
 * \param ae coefficient of A matrix east node
 * \param aw coefficient of A matrix west node 
 * \param an coefficient of A matrix north node
 * \param as coefficient of A matrix south node
 * \param ap coefficient of A matrix P node
 * \param bm coefficient of b vector right side
 * \param ni size of the A matrix in x direction
 * \param nj size of the A matrix in y direction
 * \param enum_meth 
 * \returns EXIT_SUCCESS on solution
 * \returns u vector solution of the linear equation 
 */ 
int
NAGSOL(double *ae, double *aw, double *an, double *as, double *ap, double *bm,double *u, int ni, int nj, char* enum_meth, char* enum_piv, char* enum_fact, int maxiter)
{
    FILE *fp;
    double *a=0; 
    double *b=0;
    double *x=0;
    double rnorm, dtol, tol;
    Integer *irow, *icol;
    Integer *istr=0, *idiag, *ipivp=0, *ipivq=0;
    Integer i, j, m, n, nnzc;
    Integer lfill, npivm;
    Integer maxitn;
    Integer itn;
    Integer nnz = 0;
    Integer num;
    int Nx;
    int Index;

    Nag_SparseNsym_Method method;
    Nag_SparseNsym_Piv pstrat;
    Nag_SparseNsym_Fact milu;
    Nag_Sparse_Comm comm;

    Nag_SparseNsym_Zeros zero;
    Nag_SparseNsym_Dups dup;

    lfill = 0;
    tol = 1.0e-5;
    m = 4;
    dtol = 0.0;

    maxitn = maxiter;
    //char char_enum[20];
    n =(ni-2)*(nj-2);

    // allocate matrix a
    //printf("size of matrix: %i \n", sizeof(double)*n*n);
    //printf("size of n: %i\n", n);
    
    
    
    Nx = (csim2d->getDimension())->Nx;
    Index = (Nx+1); 
    
    // find out how many non-zero entries
    for (i = 1; i <= ni-2; i++) {
        for (j = 1; j <= nj-2; j++) {
            int k = Index-Nx-2*(i-1);

            // this dangerous due the problem of comparing
            // null with real number.
	    if (as[Index] != 0.0) {
				nnz++;
            }
            if (aw[Index] != 0.0) {
				nnz++;
            }
	    if (ap[Index] != 0.0) {
				nnz++;
            }
	    if (ae[Index] != 0.0) {
				nnz++;
            }
            if (an[Index] != 0.0) {
				nnz++;
            }
            Index++;
        }
        Index +=2;
    }
	
    // choose appropriate solver
    if (!strcmp(enum_meth, "RGMRES"))
        method = Nag_SparseNsym_RGMRES;
    else if(!strcmp(enum_meth, "CGS"))
        method = Nag_SparseNsym_CGS;
    else if(!strcmp(enum_meth, "BiCGSTAB"))
        method = Nag_SparseNsym_BiCGSTAB;
    else {
        Vprintf("Unrecognized string for method solver.\n");
        return EXIT_FAILURE;
    }
    
    // choose appropriate pivoting strat.
    if (!strcmp(enum_piv, "NoPiv"))
        pstrat = Nag_SparseNsym_NoPiv;
    else if (!strcmp(enum_piv, "UserPiv"))
        pstrat = Nag_SparseNsym_UserPiv;
    else if (!strcmp(enum_piv, "PartialPiv"))
        pstrat = Nag_SparseNsym_PartialPiv;
    else if (!strcmp(enum_piv, "CompletePiv"))
        pstrat = Nag_SparseNsym_CompletePiv;
    else {
        Vprintf("Unrecognised string for pstrat");
        return EXIT_FAILURE;
    }
    
    // choose appropriate factorization
    if (!strcmp(enum_fact, "ModFact"))
        milu = Nag_SparseNsym_ModFact;
    else if(!strcmp(enum_fact, "UnModFact"))
        milu = Nag_SparseNsym_UnModFact;
    else {
        Vprintf("Unrecognised string for milu enum representation.\n");
        return EXIT_FAILURE;
    }
    
    // allocate matrix
    num = 2*nnz;
    istr = NAG_ALLOC(n+1, Integer);
    idiag = NAG_ALLOC(n, Integer);
    ipivp = NAG_ALLOC(n, Integer);
    ipivq = NAG_ALLOC(n, Integer);
    x = NAG_ALLOC(n, double);
    b = NAG_ALLOC(n, double);
    irow  = NAG_ALLOC(num, Integer);
    icol  = NAG_ALLOC(num, Integer);
	a     = NAG_ALLOC(num, double);

	if (!a || !icol || !irow || !x || !b || !istr || !idiag || !ipivp || !ipivq) { 
        Vprintf("Allocation failure for nag matrix.\n");
        return EXIT_FAILURE;
    }
	//Vprintf("Number of non-zero: %ld\n",nnz);
	nnz = 0;
	Index = (Nx+1); 
	// convert the to NAG COO format
    for (i = 1; i <= ni-2; i++) {
        for (j = 1; j <= nj-2; j++) {
            Integer k = Index-Nx-2*(i-1);

            // this dangerous due the problem of comparing
            // null with real number.
	    if (as[Index] != 0.0 && k-Nx-2 > 0) {
				a[nnz] = as[Index];
				irow[nnz] = k; icol[nnz] = k-Nx+2;
				nnz++;
			}
            if (aw[Index] != 0.0 && k-1 > 0) {
				a[nnz] = aw[Index];
				irow[nnz] = k; icol[nnz] = k-1;
				nnz++;	
			}
	    if (ap[Index] != 0.0) {
				a[nnz] = ap[Index];
				irow[nnz] = k; icol[nnz] = k;
				nnz++;
			}
	    if (ae[Index] != 0.0 && k+1 < n) {
				a[nnz] = ae[Index];
				irow[nnz] = k; icol[nnz] = k+1;
				nnz++; 
			}
            if (an[Index] != 0.0 && k+Nx+2 < n) {
				a[nnz] = an[Index];
				irow[nnz] = k; icol[nnz] = k+Nx-2;
				nnz++;
			}
            Index++;
        }
        Index +=2;
    }
	/*for (int i = 0; i < nnz; i++)
		Vprintf("irow[%ld]=%ld, icol[%ld]=%ld, nnz=%ld\n", i, irow[i], i, icol[i],nnz);*/

    // Read in right-hand side vector b and initialize approximate solution
    int k = 0;
    Index = Nx + 1;
    for (i = 1; i <= ni-2; i++) {
        for (j = 1; j <= nj-2; j++) {
            b[k] = bm[Index];
            x[k] = u[Index];
            Index++;
            k++;
        }
        Index += 2;
    }
/*
	fp = fopen("amat", "w");
	// Output matrix and b vector for debugging.
	for (i = 1; i <= nnz; ++i)
		fprintf(fp, " %8ld%16.4e%8ld%8ld\n",i,a[i-1],irow[i-1],icol[i-1]);
	fprintf(fp,"\n");
	fclose(fp);
*/
	dup = Nag_SparseNsym_SumDups;
	zero = Nag_SparseNsym_RemoveZeros;
	

	// order the elements in the matrix 
	f11zac(n, &nnz, a, irow, icol, dup, zero, istr, NAGERR_DEFAULT);

    // Calculate incomplete LU factorization
	//printf("n:%i,nnz:%i,num:%i\n",n,nnz,num);
	f11dac(n, nnz, &a, &num, &irow, &icol, lfill, dtol,
          pstrat, milu, ipivp, ipivq, istr, idiag, &nnzc, &npivm, NAGERR_DEFAULT);
	
    // Solve Ax=b using F11DCC
    f11dcc(method, n, nnz, a, num, irow, icol, ipivp, ipivq, istr,
            idiag, b, m, tol, maxitn, x, &rnorm, &itn, &comm, NAGERR_DEFAULT);
    Vprintf("%s%10ld%s\n","Converged in",itn, " iterations");
    Vprintf("%s%16.3e\n","Final residual norm =",rnorm);

    // Output x
    k = 0;
    Index = Nx + 1;
    for (int i = 1; i <= (ni-2); i++) {
        for (j = 1; j <= (nj-2); j++) {
            u[Index] = x[k];
            Index++;
            k++;
        }
        Index += 2;
    }

    // free memory
    if (istr)  NAG_FREE(istr);
    if (idiag) NAG_FREE(idiag);
    if (ipivp) NAG_FREE(ipivp);
    if (ipivq) NAG_FREE(ipivq);
    if (irow)  NAG_FREE(irow);
    if (icol)  NAG_FREE(icol);
    if (a)     NAG_FREE(a);
    if (x)     NAG_FREE(x);
    if (    NAG_FREE(b);

    return EXIT_SUCCESS;
}

#endif 
