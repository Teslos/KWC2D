/*
*	After trying to solve concentration equation using different solvers, 
*   decision is made to return to iterative Gauss-Seidel method which 
*	showed as stable for concentration solutions, until better method is found.
*
*   Simpler test with the Gauss-Seidel method using normal diffusion equation
*	solves theta equation. So Gauss-Seidel method is working correctly. 
*/
#include <iostream>
#include <math.h>

#include "sim2d.h"

// solvers
#include "cghs.h"
#include "bicgsq.h"
#include "bicgstab.h"
#include "gmres.h"
#include "cblas.h"

#define CGHS     1
#define BICGSQ   2
#define BICGSTAB 3
#define GMRES    4
#define PI       3.14159
#define PI2      6.28319
#define FSMALL   1e-6
//#define DEBUG 1

double EPS  =  1e-6;

using namespace std;

extern Csim2d *csim2d;
struct FiveDiagMatrix {
    int N;
    double *ap, *ae, *aw, *an, *as;
    FiveDiagMatrix( int mN, double *Ap, double *Ae, double *Aw,
            double *An, double *As) : 
        N(mN), ap(Ap), ae(Ae), aw(Aw), an(An), as(As){}
};

void 
mult(const FiveDiagMatrix &A, const double *v, double *w)
{
    int i,j, Index;
    int Nx = (csim2d->getDimension())->Nx;
    int Ny = (csim2d->getDimension())->Ny;
    Index = Nx + 1;
    for (i = 0; i < Nx * Ny; i++) w[i] = 0.; // restart to zero
    for (i = 0; i < Ny-2; i++) {
        for (j = 0; j < Nx-2; j++) {
           w[Index] = A.ap[Index]*v[Index]+A.ae[Index-1]*v[Index-1]+
               A.aw[Index+1]*v[Index+1]+A.an[Index+Nx]*v[Index+Nx]
               +A.as[Index-Nx]*v[Index-Nx];
           // go to next cell in x-direction
           Index++;
        }
        // finished in x-direction go 
        // further in y-direction one above
        Index += 2;
    }
}

struct DiagMatrix {
    int N;
    double *ap, *ae, *aw, *an, *as, *ane, *anw, *ase, *asw;
    DiagMatrix( int mN, double *Ap, double *Ae, double *Aw,
            double *An, double *As, double *Ane, double *Anw, double *Ase,
            double *Asw) : 
        N(mN), ap(Ap), ae(Ae), aw(Aw), an(An), as(As), ane(Ane),
    anw(Anw), ase(Ase), asw(Asw){}
};

void 
mult(const DiagMatrix &A, const double *v, double *w)
{
    int i,j, Index;
    // dimension of the matrix are smaller because 
    // of boundary cells.
    int Nx = (csim2d->getDimension())->Nx;
    int Ny = (csim2d->getDimension())->Ny;
    Index = (Nx + 1);
    for (i = 0; i < Nx * Ny; i++) w[i] = 0.; // restart to zero
    for (i = 0; i < Ny-2; i++) {
        for (j = 0; j < Nx-2; j++) {
           w[Index] = A.ap[Index]*v[Index]+A.ae[Index-1]*v[Index-1]+
               A.aw[Index+1]*v[Index+1]+A.an[Index+Nx]*v[Index+Nx]
               +A.as[Index-Nx]*v[Index-Nx]+A.ane[Index+Nx-1]*v[Index+Nx-1]
               +A.anw[Index+Nx+1]*v[Index+Nx+1]+A.ase[Index-Nx-1]*v[Index-Nx-1]
               +A.asw[Index-Nx+1]*v[Index-Nx+1];
           // go to next cell in x-direction
           Index++;
        }
        // finished in x-direction go 
        // further in y-direction one above
        Index += 2;
    }
}
/* new conjugate based solver */
int CGM(double *ae, double *aw, double *an, double *as, double *ap, double *b, double *u, int ni, int nj, int maxit, int method)
{
    int n, i, j, ij; // counters
    int m = 5;       // only for GMRES solver
    int M = ni*nj;
    int its;
    
    // solve equation Au=b
    FiveDiagMatrix A(ni,ap,ae,aw,an,as);
    //mult(A,x,b);
    its = -1;
    switch (method)
    {
        case CGHS:
        its = cghs(M,A,b,u,EPS);
        break;
        case BICGSQ:
        //its = bicgsq(M,A,b,u,EPS);
        break;
        case BICGSTAB:
        //its = bicgstab(M,A,b,u,EPS);
        break;
        case GMRES:
        its = gmres(m,M,A,b,u,EPS);
        break;
        default:
        cerr << "Error running conjugate solver" << endl;
        break;
    }
    return its;
    cout << "finished in " << its << " iterations" << 
        " difference is " << EPS << endl;
    //for (int i=0; i < M; i++) x[i] = u[i];
    
}

void CGM(double *ae, double *aw, double *an, double *as, double *ane, double *anw, double *ase, double *asw,
        double *ap, double *b, double *u, int ni, int nj, int maxit, int method)
{
    int n, i, j, ij; // counters
    int m = 5;       // only for GMRES solver
    int M = ni*nj;
    int its;
    
    // solve equation Au=b
    DiagMatrix A(ni,ap,ae,aw,an,as,ane,anw,ase,asw);
    //mult(A,x,b);
    its = -1;
    switch (method)
    {
        case CGHS:
        its = cghs(M,A,b,u,EPS);
        break;
        case BICGSQ:
        //its = bicgsq(M,A,b,u,EPS);
        break;
        case BICGSTAB:
        //its = bicgstab(M,A,b,u,EPS);
        break;
        case GMRES:
        its = gmres(m,M,A,b,u,EPS);
        break;
        default:
        cerr << "Error running conjugate solver" << endl;
        break;
    }
    cout << "finished in " << its << " iterations" << 
        " difference is " << EPS << endl;
    //for (int i=0; i < M; i++) x[i] = u[i];
    
}

//This routine incorporates the Incomplete Cholesky preconditioned
//Conjugate Gradient Solver for symmetric matrices in 2d problems
//with 5 diagonal matrix structure.

// translation of fortran code by Ismet Demirdzic, Sarajevo, 1991.	

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

void slapsol(double *ae, double *aw, double *an, double *as, double *ap, double *b, double *u, int ni, int nj, int maxit)
{
    char *method;
    slapsol_(ae,aw,an,as,ap,b,u,&ni,&nj,&maxit,method);
}

int 
printout(double *ae, double *aw, double *an, double *as, double *ap, double *b, int ni, int nj)
{
    int i,j;
    int nnz=0;
    int n = ni;
    int kl = 2;
    int ku = 2;
    int nrhs = 1;
    FILE *matout;
    matout = fopen("matout", "w+");
    if (matout == NULL) {
        fprintf(stderr,"Can not open file matout\n");
    }
    // prints out the matrix in the nag format
    fprintf(matout, "%d %d %d %d\n", n, kl, ku, nrhs);
    int Nx = (csim2d->getDimension())->Nx;
    int Index = Nx + 1;
    for (i=1; i <=n-2; i++) {
        for (j=1; j<=n-2; j++) {
            int k = Index-Nx-2*(i-1);
            if (as[Index-Nx] != 0.0) {
                fprintf(matout, "%5lf %d %d \n",as[Index-Nx], k, k-(Nx-2));
                nnz++;
            }
            if (ae[Index-1] != 0.0) {
                fprintf(matout, "%5lf %d %d \n",ae[Index-1], k, k-1);
                nnz++;
            }
            if (ap[Index] != 0.0) {
                fprintf(matout, "%5lf %d %d \n",ap[Index], k, k);
                nnz++;
            }
            if (aw[Index+1] != 0.0) {
                fprintf(matout, "%5lf %d %d \n",aw[Index+1], k, k+1);
                nnz++;
            }
            if (an[Index+Nx] != 0.0) {
                fprintf(matout, "%5lf %d %d \n",an[Index+Nx], k, k+Nx-2);
                nnz++;
            }
            Index++;
        }
        Index +=2;
    }
    fprintf(matout,"Non zero elements: %d\n",nnz);
    // write out b vector
    Index = Nx+1;
    for (i=1; i<=n-2; i++) {
        for (j=1; j<=n-2; j++) {
            fprintf(matout, "%5lf \n",b[Index]);
            Index++;
        }
        Index +=2;
    }
    
    fclose(matout);
	return 0;
}
    
int 
Csim2d::alloc_c()
{
    ap = new double[DimInfo->Nxy];
    if (ap == NULL) {
        cerr << "Error allocating memory for ap coeff." << endl;
        return -1;
    }
    ae = new double[DimInfo->Nxy];
    if (ae == NULL) {
        cerr << "Error allocating memory for ae coeff." << endl;
        return -1;
    }
    aw = new double[DimInfo->Nxy];
    if (aw == NULL) {
        cerr << "Error allocating memory for aw coeff." << endl;
        return -1;
    }
    an = new double[DimInfo->Nxy];
    if (an == NULL) {
        cerr << "Error allocating memory for an coeff." << endl;
        return -1;
    }
    as = new double[DimInfo->Nxy];
    if (as == NULL) {
        cerr << "Error allocating memory for as coeff." << endl;
        return -1;
    }
    ane = new double[DimInfo->Nxy];
    if (ane == NULL) {
        cerr << "Error allocating memory for ane coeff." << endl;
        return -1;
    }
    anw = new double[DimInfo->Nxy];
    if (anw == NULL) {
        cerr << "Error allocating memory for anw coeff." << endl;
        return -1;
    }
    ase = new double[DimInfo->Nxy];
    if (ase == NULL) {
        cerr << "Error allocating memory for ase coeff." << endl;
        return -1;
    }
    asw = new double[DimInfo->Nxy];
    if (asw == NULL) {
        cerr << "Error allocating memory for asw coeff." << endl;
        return -1;
    }
    b = new double[DimInfo->Nxy];
    if (b == NULL) {
        cerr << "Error allocating memory for b coeff." << endl;
        return -1;
    }
    u = new double[DimInfo->Nxy];
    if (u == NULL) {
        cerr << "Error allocating memory for u coeff." << endl;
        return -1;
    }
    return 0;
}

void
Csim2d::resetv(void)
{
    // set all vectors to zero
    for (int i = 0; i < DimInfo->Nxy; i++) {
        ae[i] = aw[i] = an[i] = as[i] = ap[i] = 0.;
        ane[i] = anw[i] = ase[i] = asw[i] = b[i] = u[i] = Theta[i] = 0.;
    }
}

void
Csim2d::dealloc_c()
{
    delete [] ap;
    delete [] ae;
    delete [] aw;
    delete [] an;
    delete [] as;
    delete [] ane;
    delete [] anw;
    delete [] ase;
    delete [] asw;
    delete [] b;
    delete [] u;
}

inline double 
Csim2d::Igamma(double zeta)
{
    if ((0 <= zeta) && (zeta < (1.0/gamma)))
        return gamma;
    else
        return (1./zeta);
}

inline double 
Csim2d::Agamma(double zeta)
{
    if (0. <= zeta && zeta < 1./gamma)
        return gamma/2 * zeta*zeta;
    else
        return zeta - 1./(2.*gamma);
}

// subtitute now returns correction term
// this is described in eq (A.29) Acta Materialia 51 (2003) 6035–6058.
inline double 
Csim2d::Subt(double Theta)
{
    if (-PI <= Theta && Theta < PI)
        return 0.0;
    else if (Theta < -PI)
        return PI2;
    else if (Theta > PI)
        return -PI2;
    return Theta; // this is liquid 
}

inline double 
Csim2d::P(double w)
{
    double t = exp(-beta * w);
	//cout << " P value is : " << 1.<< endl;
    return 1. - t + (mu / epsilonL) * t;
	//return 1.;
}

int
Csim2d::solve2d(void)
{
    // constants
    double const normConst = 1.e-8;
    int const maxit = 1000;
    int i,j;
	int type = 2;      // calculate full KWC model
    float x_maxPos = 0, y_maxPos = 0; 
    // coefficients of the "stiffness" matrix
    double atc = 0., center = 0.,
		   dpxl = 0., dpxu = 0., dpxd = 0., dpxr = 0., dpyr = 0., 
		   dpyl = 0., dpyu = 0., dpyd = 0., normr = 0., 
		   norml = 0., normu = 0., normd = 0.;
    double poly = 0.0,
           phi_x = 0.0,
           phi_y = 0.0,
           phi_xx = 0.0,
           phi_yy = 0.0,
           phi_xy = 0.0, co = 0., si = 0., thx = 0., thy = 0., anis = 0.;
    double FSold = 0.;
    double FSNew = 0.;
	double CNew  = 0.;
	double Tnew  = 0.;
    double laplace = 0.;
    double norm2 = 0.;
    double norm2R = 0.;
	double error = 0.;
	
    // calculate values
    double deltaR = 1. / Delta;
    double deltaQR = 72. /(Delta * Delta);
    double dt_mobil = TimeInfo->tWidth * mobil;
    double a12 = 0.6267 * 5. / 24. ;
    a12 *= Delta * ML * (1. - k) / DC[1];
    double VMax = V_control;

    // calculate reciprocal quadrats
    double XR = 1. / LengthX;
    double twoXR = 0.5 * XR;
    double YR = 1. / LengthY;
    double TR = 1. / TimeInfo->tWidth;

    double XRQuadrat = XR * XR;
    double YRQuadrat = YR * YR;

    double fourXRQuadrat = XRQuadrat * 0.25;
    double atcm = 0.5 * ATC * Delta * XR * TR;

    double epsilonLQ = epsilonL * epsilonL;
    double alfaQ = alfa * alfa;

    int Nx = DimInfo->Nx;
    int Ny = DimInfo->Ny;
    int Nxy = DimInfo->Nxy;

    double Dsk = k * DC[0];
    double Dls = DC[1] - Dsk;
    double p = 1. - k;
    
    double KV1 = TimeInfo->tWidth;
    KV1 *= XRQuadrat;
    KV1 *= my;

    double kp1 = -LengthX * LengthX * Rho * TR;
    double kp1R = 1. / kp1;
    double kp2 = TimeInfo->tWidth / (Rho * LengthX);

    // summing up all liquid contributions for averaging liquid conc */
    // reset the old values
    double C_LA = 0.;
    double S_FL = 0.;
    double S_d_FS = 0.;
    int l      = 0;

    // reset all the coefficients
    resetv();
    
    int Index = Nx + 1;
    for (i = 0; i < numY; i++) {
        double C_LAR = 0.;
        double S_LAR = 0.;
        for (j = 0; j < numX; j++) {
            Index = Nx * (i + 1) + j + 1;
            if (FL[Index] > .5) {
                C_LAR += CellInfo[Index].Conc;
                l++;
            }
        }
        C_LA += C_LAR;
        S_FL += S_LAR;
    }
    C_LA /= (double) l;
    S_FL /= (double) l;
	
	// first test values
	/*fprintf(stdout,"FLmx, FLmy\n");
	for (i = 0; i < Nxy; i++) {
		fprintf(stdout, "%i: %lf %lf\n", i, FLmx[i], FLmy[i]);
	}*/
	
	
	// start point for flmx: right down and then one row up 
	Index = 2 * Nx - 1;
    for (i = 0; i < Ny; i++) {
        FLmx[Index] = .5 * (FL[Index] + FL[Index - 1]);
        Index += Nx;
    }

	// start point for flmy: left up and then right up
	// correction: 19.02.14 Ny is now numY
    Index = Nx * (numY + 1) + 1;
    for (i = 0; i < Nx; i++) {
        FLmy[Index] = .5 * (FL[Index] + FL[Index - Nx]);
        Index++;
    }
	
    // calculate fraction liquid at middle of the cell edge(x-dir).
    for (i = 1; i < Nxy; i++) {
        FLmx[i] = .5 * (FL[i] + FL[i-1]);
        FSmx[i] = 1.-FLmx[i];
    }

    // calculate fraction liquid at middle of the cell edge(y-dir). 
    for (i = Nx+1; i < Nxy; i++) {
        FLmy[i] = .5 * (FL[i] + FL[i - Nx]);
        FSmy[i] = 1.-FLmy[i];
    }

	// write test values
	/*for (i = 0; i < Nxy; i++) {
		fprintf(stdout, "%i %lf %lf\n", i, FLmx[i], FLmy[i]);
	}
		
	exit(1);*/
    // calculate gradient of theta at middle of cell edge Tmx:
    // absolute value of Tmx is calculated.
    for (i = 1; i < Nxy; i++) {
        Tmx[i] = CellInfo[i].Theta-CellInfo[i-1].Theta;
        Cmx[i] = Subt(Tmx[i]);
        Tmx[i] += Cmx[i];
		Tmx[i] *= XR;
        Tmx[i] = fabs(Tmx[i]);
		//cout << "Tmx[" << i << "]=" << Tmx[i] << " Theta[" << i << "]=" << CellInfo[i].Theta 
		//		<< " Theta[" << i-1 << "]=" << CellInfo[i-1].Theta << " XR="<< XR << endl;
    }
	
    // calculate gradient of theta at middle of cell edge Tmy: 
    // absolute value of Tmy is calculated.
    for (i = Nx+1; i < Nxy; i++) {
        Tmy[i] = CellInfo[i].Theta-CellInfo[i-Nx].Theta;
        Cmy[i] = Subt(Tmy[i]);
        Tmy[i] += Cmy[i];
		Tmy[i] *= XR;
        Tmy[i] = fabs(Tmy[i]);
		//cout << "Tmy[" << i <<"]=" << Tmy[i]<< " Theta[" << i << "]=" << CellInfo[i].Theta << endl;
    }
	//exit(-1);
    
    // friction term as function of frac-liquid in calculation of 
    // implicit solver for Velocity.
    for (i = 0; i < Nxy; i++) {
        // XXX leaving for new implementation
    }

    //cpFL(DimInfo->Nx-2, 0); // on left (negative x axis)
    
    // iterative solution for velocity
    C_liquid = C_LA;

    l = 0;
    int m = 0;
    int n = 0;
    float stheta = 0;
    float epstheta = 0;
	
	do {
		Index = (Nx + 1);
		F_error = 0.; 
		//cout << "Mobil=" << mobil << " Gamma="<<Gamma
		//    <<" dtime=" << TimeInfo->tWidth <<" XRQuadrat="<<XRQuadrat<<" Delta="<<Delta<<endl;
		// calculate FS, Concentration 
		// calculate starting points for iteration of phase field
		// XXX implementation of moving frame
		for (i = 0; i < numY; i++) {
			for (j = 0; j < numX; j++) {
				FSold = CellInfo[Index].FracSol;
				float TM = Tmelt - ML * CellInfo[Index].Conc;
				float U  = (TM - (T_D + i * TGrad));
			
				// correction 30.01.14 poly was not defined 
				poly = FSold * ( 1. - FSold );
				// redefine sigma using old parameters
				sigma = tauPhi * 30 * U * deltaR;
				float Mpoly = FSold - 0.5 + sigma * poly;
				FSNew = CellInfo[Index+1].FracSol - FSold;
				FSNew += CellInfo[Index-1].FracSol - FSold;
				FSNew += CellInfo[Index+Nx].FracSol - FSold;
				FSNew += CellInfo[Index-Nx].FracSol - FSold;
				laplace = CellInfo[Index+1+Nx].FracSol - FSold;
				laplace += CellInfo[Index-1-Nx].FracSol - FSold;
				laplace += CellInfo[Index+Nx-1].FracSol - FSold;
				laplace += CellInfo[Index-Nx+1].FracSol - FSold;
				FSNew = 2. * FSNew + 0.5*laplace;
				FSNew /=3.;
				FSNew *= XRQuadrat;
				laplace = FSNew;
				// for now use laplace 8-star 
				if ( poly > (.01)) {
					phi_x = CellInfo[Index+1].FracSol - CellInfo[Index-1].FracSol;
					phi_x *= twoXR;
					phi_y = CellInfo[Index+Nx].FracSol - CellInfo[Index-Nx].FracSol;
					phi_y *= twoXR;
					norm2 = phi_x * phi_x + phi_y * phi_y;
					norm2R = 1. / norm2;
					phi_xx = CellInfo[Index+1].FracSol + CellInfo[Index-1].FracSol
						-2. * FSold;
					phi_xx *= XRQuadrat;
					phi_yy = CellInfo[Index+Nx].FracSol + CellInfo[Index-Nx].FracSol
						-2. * FSold;
					phi_yy *= XRQuadrat;
					phi_xy = CellInfo[Index+Nx+1].FracSol + CellInfo[Index-Nx-1].FracSol
						- CellInfo[Index+Nx-1].FracSol - CellInfo[Index-Nx+1].FracSol;
					phi_xy *= fourXRQuadrat;
					// correction 31.01.14 added phi_y to second term
					si = 4. *(phi_x * phi_x * phi_x * phi_y - phi_y * phi_y * phi_y * phi_x)
						* norm2R * norm2R;
					co = 1. - 8. * phi_x * phi_x * phi_y * phi_y * norm2R * norm2R;
					thx = phi_x * phi_xy - phi_y * phi_xx;
					thx *= norm2R;
					// correction 31.01.14 changed phi_yy to phi_xy
					thy = phi_x * phi_yy - phi_y * phi_xy;
					thy *= norm2R;
					anis  = co * (2. + epssigma * co) * laplace;
					// anis  = co * (2. + epssigma * co) * (phi_xx + phi_yy);
					anis -= 8. * si * (1. + epssigma * co) * (thx * phi_x + thy * phi_y);
					anis -= 16. * (co + epssigma * (co * co - si * si)) * (thy * phi_x - thx * phi_y);
					anis *= epssigma;
                
				}
				else { // no anisotropy
					anis = 0;
					co   = 0;
				}
			
				if ( type == 1 ) {
					FSNew += anis;
					FSNew += poly * (FSold - 0.5) * deltaQR;
					FSNew *= Gamma;
					FSNew += U * 30. * poly * poly * deltaR;
					dt_mobil  = 1. / (mobil * (epskin * co + 1.));
					dt_mobil += a12 * CellInfo[Index].Conc * (epssigma * co + 1.)
						*(epssigma * co + 1.);
					dt_mobil = TimeInfo->tWidth / dt_mobil;
					FSNew *= dt_mobil;
					FSNew += FSold;  
					FracSolNew[Index] = FSNew;
					
				} 
				else {
					// theta contributions 
					anis = 0.;  // for now anisotropy is zero
					FSNew += anis;
					// multiply laplace operator: alpha^2 * nabla^2 Phi
					FSNew *= mobil * Gamma * tauPhi; // this is equal to alpha^2 in KWC model
					//printf("alphaSQ: %lf\n", mobil * Gamma * tauPhi);
					FSNew += FSold * tauPhi * TR;
					// gradient theta field (here evaluate Tmx at center of cell) check if this is ok!!
					double gTmx = .5 * (Tmx[Index] + Tmx[Index+1]);
					double gTmy = .5 * (Tmy[Index] + Tmy[Index+Nx]);
					stheta = sqrt(gTmx*gTmx+gTmy*gTmy);
					epstheta = stheta * stheta;
					epstheta *= epsilonLQ;
					stheta *= 2*s;
					if (Mpoly > 0) {
						theta_over[Index] = tauPhi*TR + FSold*Mpoly + stheta + epstheta;
						FSNew += FSold * Mpoly;
					}
					else {
						theta_over[Index] = tauPhi*TR - (1.-FSold)*Mpoly + stheta + epstheta;
					}
					FSNew /= theta_over[Index];
				
				
					// Estimation of maximum relativ error againg last iteration step
					error = fabs(( FSNew - FracSolNew[Index] ) / (FSNew+1E-20));
					//fprintf(stdout, " Error : %f F[%i]=%lf\n", error, Index, FSNew);
					F_error = (F_error > error) ? F_error : error;
				
					FracSolNew[Index] = FSNew;
				}
				Index++;
			}
			Index += 2;
		}
		// Set boundary 24.06.14, maybe not necessary it is called in main routine
		// Extruder();
		n++;
	} while (F_error > d_C && n < .5 * sqrt(Nxy) && type == 2);
	NF = n;
	n = 0;
	//exit(-1);
    // solve equations using conjugate gradient methods (GMRES)
    //CGM(ae, aw, an, as, ane, anw, ase, asw, ap, b, u, Nx, Ny, maxit, GMRES);
#ifndef DEBUG
    //CGM(ae, aw, an, as, ap, b, u, Nx, Ny, maxit, CGHS);
    //this->OutMathematica(1);
    //exit(-1);
#endif
//#else
//    NC = DNAG(ae, aw, an, as, ap, b, u, Nx, Ny, maxit, 0);
//#endif
	
    // copy for now results in FracSolNew
    Index = (Nx+1);
    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            //FracSolNew[Index] = u[Index];
            S_d_FS += FracSolNew[Index];
            //cout << "u["<<Index<<"]="<<u[Index]<< endl;
            Index++;
        }
        Index += 2;
    }
    
    // calculate average fraction solid
    S_d_FS /= (double)(numY * numX);
    float Xmax_Pos = x_maxPos;
    float Ymax_Pos = y_maxPos;
    Utip = Tmelt - ML*C_0;
    Utip -= T_D + y_maxPos * TGrad;
    Utip /= ML * C_0 * p;
    Otip = Utip / (1. + p * Utip);
    // temperature reduction latent heat * d_FS 
    T_D = T_0 - ( TimeInfo->time * TPoint -
            n_move * TGrad );
    FS = S_d_FS;

    // calculate fraction solid at middle of the cell edge(x-dir).
    for (i = 1; i < Nxy; i++) {
        FSmx[i] = .5 * (FracSolNew[i] + FracSolNew[i-1]);
    }

    // calculate fraction solid at middle of the cell edge(y-dir). 
    for (i = Nx+1; i < Nxy; i++) {
        FSmy[i] = .5 * (FracSolNew[i] + FracSolNew[i - Nx]);
    }

    // calculate phase field for orientation
    // 1. part calculate constants
    Index = (Nx+1);
    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            Dmx[Index] = FSmx[Index] * FSmx[Index] * (s*Igamma(Tmx[Index])+epsilonLQ);
            Dmy[Index] = FSmy[Index] * FSmy[Index] * (s*Igamma(Tmy[Index])+epsilonLQ);
			//cout << "Dmx["<<Index<<"]="<<Dmx[Index]<< endl;
            //cout << "Dmy["<<Index<<"]="<<Dmy[Index]<< endl;
			// testing 
			//Dmx[Index] = 0.5E-4;
			//Dmy[Index] = 0.5E-4;
			Index++;
        }
        Index += 2;
    }
	//exit(-1);
	double gTmx = 0.;
	double gTmy = 0.;
	Index = (Nx+1);
    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
			// add small constant to fract.solid 
			FracSolNew[Index] += FSMALL;
            Dmbar[Index] = (Dmx[Index]+Dmx[Index+1]+Dmy[Index]+Dmy[Index+Nx]);
			Dmbar[Index] *= XRQuadrat;
			// gradient at center 
			gTmx = .5 * (Tmx[Index] + Tmx[Index+1]);
			gTmy = .5 * (Tmy[Index] + Tmy[Index+Nx]);
			
			Dmbar[Index] += P(epsilonL*sqrt(gTmx*gTmx + gTmy*gTmy)) *
					FracSolNew[Index]*FracSolNew[Index]*TR*tauTheta;
			//cout << "I have tshirt"<<endl;
			Pbar[Index] = P(epsilonL*sqrt(gTmx*gTmx + gTmy*gTmy)) 
						* FracSolNew[Index]*FracSolNew[Index];
			//testing
			//Dmbar[Index] += TR;
			T_over[Index] = 1. / Dmbar[Index];
			//cout << "Dmbar["<<Index<<"]="<<Dmbar[Index]<< endl;
            //cout << "T_over["<<Index<<"]="<<T_over[Index]<< " Dmx["<< Index <<"]="<< 
			//			Dmx[Index]<< " Dmx["<< Index+1 <<"]="<< 
			//			Dmy[Index+1] << " Dmy["<< Index <<"]="<< 
			//			Dmy[Index]<< " Dmy["<< Index+Nx <<"]="<< 
			//			Dmy[Index+Nx]<< " P=" << P(epsilonL*sqrt(Tmx[Index]*Tmx[Index]+Tmy[Index]*Tmy[Index])) << 
			//			" FracSolSq= " << FracSolNew[Index]*FracSolNew[Index] << endl;
            Index++;
        }
        Index += 2;
    }
	//exit(-1);
	
    // 2. part solve phase field theta(orientation)
	n = 0;
	Index = (Nx+1);
	T_error = 0.0;
	for (i = 0; i < numY; i++) {
		for (j = 0; j < numX; j++) {
			ae[Index] = -Dmx[Index+1];
			ae[Index] *= XRQuadrat;
			aw[Index] = -Dmx[Index];
			aw[Index] *= XRQuadrat;
			as[Index] = -Dmy[Index];
			as[Index] *= XRQuadrat;
			an[Index] = -Dmy[Index+Nx];
			an[Index] *= XRQuadrat;
			ap[Index] = Dmbar[Index];
				
			b[Index] = tauTheta * TR * Pbar[Index] * CellInfo[Index].Theta;
			// added 19.02.14 additional corrections terms
			b[Index] += (Cmx[Index + 1]*Dmx[Index + 1] - Cmx[Index]*Dmx[Index] 
				- Cmy[Index] * Dmx[Index] + Cmy[Index + Nx] * Dmy[Index + Nx])*XRQuadrat;
				
			Index++;
		}
		Index += 2;
	}
	ExtrudTh();
	iccg( ae, aw, an, as, ap, b, Theta, Nx, Ny, maxit, 0, d_T, &NTheta, &T_error );	
	
	
	//fprintf(stdout, " T_error: %lf \n", T_error);
	//exit(-1);
    // now solve for theta field
    //this->OutCoeff();
    //exit(-1);
#ifndef DEBUG
    //NTheta = CGM(ae, aw, an, as, ap, b, Theta, Nx, Ny, maxit, GMRES);
#ifdef __NAG__	
    // using nag solver for the moment
    // NAGSOL(ae, aw, an, as, ap, b, Theta, Nx, Ny, "RGMRES", "CompletePiv", "UnModFact", maxit);
#endif
    /*
    for(i = 0; i < 25; i++) {
        ap[i] = 4.0;
        aw[i] = 1.0; ae[i] = 2.0;
        as[i] = -1.0; an[i] = 1.0;
        b[i] = 6.0;
    }
    */
    //this->OutMathematica(1);
    //exit(-1);

    //slapsol(ae, aw, an, as, ap, b, Theta, Nx, Ny, maxit);

#else
    //printout(ae, aw, an, as, ap, b, Nx, Ny);
    //NTheta = DNAG(ae, aw, an, as, ap, b, Theta, Nx, Ny, maxit, 0);
#endif

    // calculate concentration
    // 1. part calculate constants
    Index = (Nx + 1);
    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            // go with new FSnew or with middle value between old and new timestep
            oneminus[Index] = 1. - p * ( CellInfo[Index].FracSol + FracSolNew[Index] ) * 0.5;
			// correction 29.01.14 FLmy[Index - Nx]->FLmy[Index]
            Dbar[Index] = (FLmx[Index + 1] + FLmx[Index] + FLmy[Index + Nx] +
                    FLmy[Index]);
            Dbar[Index] *= Dls;
            Dbar[Index] += 4. * Dsk;
            Dbar[Index] *= XRQuadrat;
            Dbar[Index] *= TimeInfo->tWidth;

            D_FS[Index]  = FracSolNew[Index];
            D_FS[Index] -= CellInfo[Index].FracSol;
            D_FS[Index] *= p;

            C_over[Index] = 1.;
            C_over[Index] /= (oneminus[Index] + Dbar[Index] - 0.5 * D_FS[Index]);
            Index++;
        }
        // That was in x-direction
        Index += 2;
    }
    ExtrudDFS();
    // Now add anti-trapping currents
    Index = (Nx + 1);
    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            // first find gradients (right side first)
            dpxr = CellInfo[Index + 1].FracSol - CellInfo[Index].FracSol;
            dpyr = CellInfo[Index + Nx].FracSol + CellInfo[Index + Nx + 1].FracSol;
            dpyr -= CellInfo[Index - Nx].FracSol + CellInfo[Index - Nx + 1].FracSol;
            dpyr *= .25;
            float normr = sqrt(dpxr * dpxr + dpyr * dpyr);
            if (normr > normConst) {
                dpxr /= normr;
            }
            else {
                dpxr = 0.;
            }

            // now gradients for left side
            dpxl = CellInfo[Index].FracSol - CellInfo[Index - 1].FracSol;
            dpyl = CellInfo[Index + Nx].FracSol + CellInfo[Index + Nx - 1].FracSol;
            dpyl -= CellInfo[Index - Nx].FracSol + CellInfo[Index - Nx - 1].FracSol;
            dpyl *= 0.25;
            norml = sqrt(dpxl * dpxl + dpyl * dpyl);
            if (norml > normConst) {
                dpxl /= norml;
            } 
            else {
               dpxl = 0.;
            }

            // gradients for upper side
            dpxu = CellInfo[Index + Nx + 1].FracSol + CellInfo[Index + 1].FracSol;
            dpxu -= CellInfo[Index + Nx - 1].FracSol + CellInfo[Index - 1].FracSol;
            dpxu *= .25;
            dpyu = CellInfo[Index + Nx].FracSol - CellInfo[Index].FracSol;
            normu = sqrt(dpxu * dpxu + dpyu * dpyu);
            if (normu > normConst) {
                dpyu /= normu;
            }
            else {
                dpyu = 0.;
            }

            // gradients for down side
            dpxd = CellInfo[Index - Nx + 1].FracSol + CellInfo[Index + 1].FracSol;
            dpxd -= CellInfo[Index - Nx - 1].FracSol + CellInfo[Index - 1].FracSol;
            dpxd *= 0.25;
            dpyd = CellInfo[Index].FracSol - CellInfo[Index - Nx].FracSol;
            normd = sqrt(dpxd * dpxd + dpyd * dpyd);
            if (normd > normConst) {
                dpyd /= normd;
            }
            else {
                dpyd = 0.;
            }

            center = D_FS[Index] * CellInfo[Index].Conc;
            atc  = (D_FS[Index + 1] * CellInfo[Index + 1].Conc + center) * dpxr;
            atc -= (D_FS[Index - 1] * CellInfo[Index - 1].Conc + center) * dpxl;
            atc += (D_FS[Index + Nx] * CellInfo[Index + Nx].Conc + center) * dpyu;
            atc -= (D_FS[Index - Nx] * CellInfo[Index - Nx].Conc + center) * dpyd;
            atc *= atcm;
            KKonv[Index] = -atc;
            // now go to next cell in x-direction
            Index++;
        }
        Index += 2;
    }
    // 2. part solve concentration
    n = 0;
	do {
		Index = (Nx + 1);
		C_error = 0.;
		for (i = 0; i < numY; i++) {
			for (j = 0; j < numX; j++) {
				CNew = u[Index + 1] * (Dls * FLmx[Index + 1] + Dsk);
				CNew += u[Index - 1] * (Dls * FLmx[Index] + Dsk);
				CNew += u[Index + Nx] * (Dls * FLmy[Index + Nx] + Dsk);
				CNew += u[Index - Nx] * (Dls * FLmy[Index] + Dsk);
				CNew *= XRQuadrat;
			
				CNew -= KKonv[Index];
				CNew *= TimeInfo->tWidth;
			
				CNew += (oneminus[Index]+ 0.5 * D_FS[Index]) * CellInfo[Index].Conc;
				CNew *= C_over[Index];
			
			
				u[Index] = CNew;
				
				// go to next cell in x-direc
				Index++;
			}
			// go to next line in y-direction to the first cell
			Index += 2;
		}
		ExtrudC();
		
		// loop backwarts for error minimazing 
		Index = (Nxy - 2 - Nx);
		for (i = 0; i < numY; i++) {
			for (j = 0; j < numX; j++) {
				CNew = u[Index + 1] * (Dls * FLmx[Index + 1] + Dsk);
				CNew += u[Index - 1] * (Dls * FLmx[Index] + Dsk);
				CNew += u[Index + Nx] * (Dls * FLmy[Index + Nx] + Dsk);
				CNew += u[Index - Nx] * (Dls * FLmy[Index] + Dsk);
				CNew *= XRQuadrat;
			
				CNew -= KKonv[Index];
				CNew *= TimeInfo->tWidth;
			
				CNew += (oneminus[Index]+ 0.5 * D_FS[Index]) * CellInfo[Index].Conc;
				CNew *= C_over[Index];
			
			
				// Estimation of maximum relativ error againg last iteration step
				error = fabs(( CNew - u[Index] ) / CNew);
				//fprintf(stdout, " Error : %f \n", error);
				C_error = (C_error > error) ? C_error : error;
				u[Index] = CNew;
				
				// go to cell left in x-direc
				Index--;
			}
			
			// go to steps to left end of deeper line 
			Index -= 2;
		}
		ExtrudC();
		n++;
	} while (C_error > d_C && n < .5 * sqrt(Nxy));
	//fprintf(stdout, "#C_E %f %d %f\n", C_error, n, d_C);
	NC = n;
	n = 0;
	
	
	
    // solve equations using conjugate gradient methods (GMRES)
    //NC = CGM(ae, aw, an, as, ap, b, u, Nx, Ny, maxit, GMRES);
    // this->OutCoeff();
#ifdef __NAG__
	NAGSOL(ae, aw, an, as, ap, b, u, Nx, Ny, "RGMRES", "CompletePiv", "UnModFact", maxit);
#endif
    //slapsol(ae, aw, an, as, ap, b, Theta, Nx, Ny, maxit);

    //NC = CGM(ae, aw, an, as, ap, b, u, Nx, Ny, maxit, CGHS);
  
    // restoring of all calculated values as new values
    for (i = 0; i < DimInfo->Nxy; i++) {
        CellInfo[i].FracSol = FracSolNew[i];
        CellInfo[i].Theta = Theta[i];
		//cout<<"Theta["<<i<<"]="<<CellInfo[i].Theta << endl;
        FL[i] = 1.-FracSolNew[i];
        CellInfo[i].Conc = u[i];
        //cout<<"c["<<i<<"]="<<CellInfo[i].Conc << endl;
        CellInfo[i].CG = CellInfo[i].Conc * (1. - FracSolNew[i]*p);
    }
	
	// forgot to add cyclic boundary conditions 24.06.14 in y direction
	//for (i = 0; i < DimInfo->Nxy - Nx; i++) {
	//	CellInfo[i].FracSol = CellInfo[i + Nx].FracSol;
	//	CellInfo[i].Theta = CellInfo[i + Nx].Theta;
	//	FracSolNew[i]  = FracSolNew[i + Nx];
	//	FL[i] = FL[i + Nx];
	//	CellInfo[i].Conc = CellInfo[i + Nx].Conc;
	//	CellInfo[i].CG = CellInfo[i + Nx].CG;    
	//}
	//Extruder();
	//ExtrudC();
    return 0;
}


             

            

            

    
        
            
                
            
            
            
                
    
    
            
    
        


















    
