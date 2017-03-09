//This routine incorporates the Incomplete Cholesky preconditioned
//Conjugate Gradient Solver for symmetric matrices in 2d problems
//with 5 diagonal matrix structure.

// translation of fortran code by Ismet Demirdzic, Sarajevo, 1991.	
#include <iostream>
#include <math.h>
#include "iccg.h"

/**
 * \param ae coefficient for east node
 * \param aw coefficient for west node
 * \param an coefficient for north node
 * \param as coefficient for south node
 * \param ap coefficient for P-node present node
 * \param q  coefficient 
 * \param fi
 * \param Nx size of the matrix in x direction
 * \param Ny size of the matrix in y direction
 * \param NS
 * \param ltest if value is 1 print sweep rsm value
 * \param resmax 
 * \param numiter returns number of iterations done by solver
 * \param resf residue of the solutions
 */
void iccg(double *ae, double *aw, double *an, double *as, double *ap, double *q, double *fi, int Nx, int Ny, int NS, int ltest, double resmax, int *numiter, float *resf)
{
	int i,j,Index, l;  
	
	int numX = Nx - 2;
	int numY = Ny - 2;
	int Nxy  = Nx * Ny;
	// allocate the working arrays
	double *pk  = (double *) malloc( Nxy * sizeof(double) ); 
	double *zk  = (double *) malloc( Nxy * sizeof(double) );
	double *d   = (double *) malloc( Nxy * sizeof(double) );
	double *res = (double *) malloc( Nxy * sizeof(double) );
	
	// initialize working arrays
	for (i = 0; i < Nxy; i++) {
		pk[i] = 0.0;
		zk[i] = 0.0;
		d[i]  = 0.0;
		res[i] = 0.0;
	}
	
	// calculate initial residual vector and the norm
	Index = Nx + 1;
	double res0 = 0;
	for (j = 0; j < numY; j++) {
		for (i = 0; i < numX; i++) {
			res[Index] = q[Index] - ap[Index]*fi[Index] - ae[Index]*fi[Index+1] - aw[Index]*fi[Index-1] -
				an[Index]*fi[Index+Nx] - as[Index]*fi[Index-Nx];
			res0 += fabs(res[Index]); 
			Index++;
		}
		Index += 2;
	}
	
	// if ltest = 1, print norm
	if (ltest) printf("SWEEP, RES0 = %lf\n", res0);
	
	Index = Nx + 1;
	//calculate elements of diagonal preconditioning matrix
	for (j = 0; j < numY; j++) {
		for (i = 0; i < numX; i++) {
			d[Index] = 1./(ap[Index] - aw[Index]*aw[Index]*d[Index-1] - as[Index]*as[Index]*d[Index-Nx]);
			Index++;
		}
		Index += 2;
	}
	
	double s0 = 1E+20;
	
	//start inner iterations
	for (l = 0; l < NS; l++) {
	
	
		//solve for zk[Index] -- forward substitution
		Index = Nx + 1;
		for (j = 0; j < numY; j++) {
			for (i = 0; i < numX; i++) {
				zk[Index] = (res[Index] - aw[Index]*zk[Index-1] - as[Index]*zk[Index-Nx])*d[Index];
				Index++;
			}
			Index += 2;
		}
		
		Index = Nx + 1;
		for (j = 0; j < numY; j++) {
			for (i = 0; i < numX; i++) {
				zk[Index] /= (d[Index]+1E-30);
				Index++;
			}
			Index += 2;
		}
	
		//backward substitution calculate scalar product sk
		double sk = 0.0;
		Index = Nxy - 2 - Nx;
		for (j = 0; j < numY; j++) {
			for (i = 0; i < numX; i++) {
				zk[Index] = (zk[Index] - ae[Index]*zk[Index+1] - an[Index]*zk[Index+Nx])*d[Index];
				sk += res[Index]*zk[Index];
				Index--;
			}
			Index -= 2;
		}
		
		//calculate beta
		double bet = sk / s0;
		
		//calculate new search vector pk
		Index = Nx + 1;
		for (j = 0; j < numY; j++) {
			for (i = 0; i < numX; i++) {
				pk[Index] = zk[Index] + bet * pk[Index];
				Index++;
			}
			Index += 2;
		}
		
		//calculate scalar product (pk.A pk ) and alpha (overwrite zk)
		Index = Nx + 1;
		double pkApk = 0.0;
		for (j = 0; j < numY; j++) {
			for (i = 0; i < numX; i++) {
				zk[Index] = ap[Index]*pk[Index] + ae[Index]*pk[Index+1] + aw[Index]*pk[Index-1] +
					an[Index]*pk[Index+Nx] + as[Index]*pk[Index-Nx];
				pkApk += pk[Index] * zk[Index];
				Index++;
			}
			Index += 2;
		}
		
		double alf = sk / pkApk;
		
		//calculate variable correction, new residual vector, and norm
		double resl = 0.;
		Index = Nx + 1;
		for(j = 0; j < numY; j++) {
			for(i = 0; i < numX; i++) {
				fi[Index]  += alf*pk[Index];
				res[Index] -= alf*zk[Index];
				resl += fabs(res[Index]);
				Index++;  
			}
			Index += 2;
		}
		
		s0 = sk;
		
		//check convergence
		double rsm = resl / (res0+1E-30);
		if (ltest) printf(" SWEEP, RESL = %lf, RSM = %lf\n", resl, rsm);
		if (resl < resmax) {
			*numiter = l;
			*resf = resl;
			free(pk);
			free(zk);
			free(d);
			free(res);
			return;
		} 
		
	}
	
}


