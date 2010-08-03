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

// Multiply Five diagonal matrix with vector v and
// returns as result vector w.
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

// Definition of Diagonal Matrix.
// It contains all neigh. coefficients including
// also extended laplacian coefficients.
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
        its = bicgsq(M,A,b,u,EPS);
        break;
        case BICGSTAB:
        its = bicgstab(M,A,b,u,EPS);
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
        its = bicgsq(M,A,b,u,EPS);
        break;
        case BICGSTAB:
        its = bicgstab(M,A,b,u,EPS);
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

// Allocates all necessary memory for
// arrays of coefficients. 
// Returns -1 error if not successed.
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

// Resets all vectors to zero.
void
Csim2d::resetv(void)
{
    // set all vectors to zero
    for (int i = 0; i < DimInfo->Nxy; i++) {
        ae[i] = aw[i] = an[i] = as[i] = ap[i] = 0.;
        ane[i] = anw[i] = ase[i] = asw[i] = b[i] = u[i] = 0.;
    }
}

// Deallocate all memory for vectors.
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

// Main routine it solves partial differential
// equations for phase-field and concentration.
int
Csim2d::solve2d(void)
{
    // constants
    double const normConst = 1.e-8;
    int const maxit = 20;
    int i,j;
    float x_maxPos = 0, y_maxPos = 0; 
    // coefficients of the "stiffness" matrix
    double atc = 0., center = 0.;
    double poly = 0.0,
           phi_x = 0.0,
           phi_y = 0.0,
           phi_xx = 0.0,
           phi_yy = 0.0,
           phi_xy = 0.0, co = 0., si = 0., thx = 0., thy = 0., anis = 0.;
    float FSold = 0.;

    // calculate values
    float deltaR = 1. / Delta;
    float deltaQR = 72. /(Delta * Delta);
    float dt_mobil = TimeInfo->tWidth * mobil;
    float a12 = 0.6267 * 5. / 24. ;
    a12 *= Delta * ML * (1. - k) / DC[1];
    float VMax = V_control;

    // calculate reciprocal quadrats
    float XR = 1. / LengthX;
    float twoXR = 0.5 * XR;
    float YR = 1. / LengthY;
    float TR = 1. / TimeInfo->tWidth;

    float XRQuadrat = XR * XR;
    float YRQuadrat = YR * YR;

    float fourXRQuadrat = XRQuadrat * 0.25;
    float atcm = 0.5 * ATC * Delta * XR * TR;

    int Nx = DimInfo->Nx;
    int Ny = DimInfo->Ny;
    int Nxy = DimInfo->Nxy;

    float Dsk = k * DC[0];
    float Dls = DC[1] - Dsk;
    float p = 1. - k;
    
    float KV1 = TimeInfo->tWidth;
    KV1 *= XRQuadrat;
    KV1 *= my;

    float kp1 = -LengthX * LengthX * Rho * TR;
    float kp1R = 1. / kp1;
    float kp2 = TimeInfo->tWidth / (Rho * LengthX);
    // summing up all liquid contributions for averaging liquid conc */
    // reset the old values
    float C_LA = 0.;
    float S_FL = 0.;
    float S_d_FS = 0.;
    int l      = 0;
    // reset all the coefficients
    resetv();
    
    int Index = Nx + 1;
    for (i = 0; i < numY; i++) {
        float C_LAR = 0.;
        float S_LAR = 0.;
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

    // start point for flmx: right down and then one row higher
    Index = 2 * Nx - 1;
    for (i = 0; i < numY; i++) {
        FLmx[Index] = .5 * (FL[Index] + FL[Index-1]);
        Index += Nx;
    }

    // start point for flmy: left up and then right up 
    Index = Nx * (numY + 1) + 1;
    for (i = 0; i < numX; i++) {
        FLmy[Index] = .5 * (FL[Index] + FL[Index - Nx]);
        Index++;
    }

    // friction term as function of frac-liquid in calculation of 
    // implicit solver for Velocity.
    for (i = 0; i < Nxy; i++) {
        // XXX leaving for new implementation
    }

    cpFL(DimInfo->Nx-2, 0); // on left (negative x axis)
    
    // iterative solution for velocity
    C_liquid = C_LA;

    l = 0;
    int m = 0;
    int n = 0;
    Index = (Nx + 1); 
    cout << "Mobil=" << mobil << " Gamma="<<Gamma
        <<" dtime=" << TimeInfo->tWidth <<" XRQuadrat="<<XRQuadrat<<" Delta="<<Delta<<endl;
    // calculate FS, Concentration 
    // calculate starting points for iteration of phase field
    // XXX implementation of moving frame
    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            FSold = CellInfo[Index].FracSol;
            float TM = Tmelt - ML * CellInfo[Index].Conc;
            float U  = (TM - (T_D + i * TGrad));
            float poly = FSold * (1. - FSold);
            // for now don't use laplace 8-star 
            if ( poly > (.01)) {
                /*
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
                si = 4. *(phi_x * phi_x * phi_x * phi_y - phi_y * phi_y * phi_x)
                    * norm2R * norm2R;
                co = 1. - 8.*phi_x*phi_x*phi_y*phi_y*norm2R*norm2R;
                thx = phi_x*phi_xy-phi_y*phi_xx;
                thx *= norm2R;
                thy = phi_x*phi_yy-phi_y*phi_yy;
                thy *= norm2R;
                anis = co*(2.+epssigma*co)*
                */  

                
            }
            else { // no anisotropy
                anis = 0;
                co   = 0;
            }
            b[Index] += anis;
            b[Index] += poly * (FSold - 0.5) * deltaQR;
            b[Index] *= Gamma;
            b[Index] += U * 30. * poly * poly * deltaR;
            dt_mobil  = 1. / (mobil * (epskin * co + 1.));
            dt_mobil += a12 * CellInfo[Index].Conc * (epssigma * co + 1.)
                *(epssigma * co + 1.);
            dt_mobil = TimeInfo->tWidth / dt_mobil;
            b[Index] *= dt_mobil;
            b[Index] += CellInfo[Index].FracSol;
            float coeff = mobil * Gamma * TimeInfo->tWidth;
            coeff *= 0.6666666*XRQuadrat; // 2/3
            b[Index] /= coeff;
            //cout << "b["<<Index<<"]="<<b[Index]<<endl;
            ap[Index] = 5.;
            ae[Index] = -1.; aw[Index] = -1.;
            an[Index] = -1.; as[Index] = -1.;
            // these are weights for extended laplacian
            ane[Index] = -.25; anw[Index] = -.25;
            ase[Index] = -.25; asw[Index] = -.25;
                
            Index++;
        }
        Index += 2;
    }
    // solve equations using conjugate gradient methods (GMRES)
    CGM(ae, aw, an, as, ane, anw, ase, asw, ap, b, u, Nx, Ny, maxit, GMRES);
    //CGM(ae, aw, an, as, ap, b, u, Nx, Ny, maxit, CGHS);
    // copy for now results in FracSolNew
    Index = (Nx+1);
    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            FracSolNew[Index] = u[Index];
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

    // calculate concentration
    // 1. part calculate constants
    Index = (Nx+1);
    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            // go with new FSnew or with middle value between old and new timestep
            oneminus[Index] = 1. - p * (CellInfo[Index].FracSol + FracSolNew[Index])*0.5;
            Dbar[Index] = (FLmx[Index+1] + FLmx[Index] + FLmy[Index+Nx] +
                    FLmy[Index-Nx]);
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
            float dpxr = CellInfo[Index + 1].FracSol - CellInfo[Index].FracSol;
            float dpyr = CellInfo[Index + Nx].FracSol + CellInfo[Index + Nx + 1].FracSol;
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
            float dpxl = CellInfo[Index].FracSol - CellInfo[Index - 1].FracSol;
            float dpyl = CellInfo[Index + Nx].FracSol + CellInfo[Index + Nx - 1].FracSol;
            dpyl -= CellInfo[Index - Nx].FracSol + CellInfo[Index - Nx - 1].FracSol;
            dpyl *= 0.25;
            float norml = sqrt(dpxl * dpxl + dpyl * dpyl);
            if (norml > normConst) {
                dpxl /= norml;
            } 
            else {
               dpxl = 0.;
            }

            // gradients for upper side
            float dpxu = CellInfo[Index + Nx + 1].FracSol + CellInfo[Index + 1].FracSol;
            dpxu -= CellInfo[Index + Nx - 1].FracSol + CellInfo[Index - 1].FracSol;
            dpxu *= .25;
            float dpyu = CellInfo[Index + Nx].FracSol - CellInfo[Index].FracSol;
            float normu = sqrt(dpxu * dpxu + dpyu * dpyu);
            if (normu > normConst) {
                dpyu /= normu;
            }
            else {
                dpyu = 0.;
            }

            // gradients for down side
            float dpxd = CellInfo[Index - Nx + 1].FracSol + CellInfo[Index + 1].FracSol;
            dpxd -= CellInfo[Index - Nx - 1].FracSol + CellInfo[Index - 1].FracSol;
            dpxd *= 0.25;
            float dpyd = CellInfo[Index].FracSol - CellInfo[Index - Nx].FracSol;
            float normd = sqrt(dpxd * dpxd + dpyd * dpyd);
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
    Index = (Nx + 1);
    C_error = 0.;
    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            ae[Index] = (Dls * FLmx[Index] + Dsk);
            ae[Index] *= XRQuadrat * TimeInfo->tWidth;
            aw[Index] = (Dls * FLmx[Index + 1] + Dsk);
            aw[Index] *= XRQuadrat * TimeInfo->tWidth;
            an[Index] = (Dls * FLmy[Index+Nx] + Dsk);
            an[Index] *= XRQuadrat * TimeInfo->tWidth;
            as[Index] = (Dls * FLmy[Index-Nx] + Dsk);
            as[Index] *= XRQuadrat * TimeInfo->tWidth;

            b[Index] = KKonv[Index] * TimeInfo->tWidth;
            ap[Index] = (oneminus[Index]+0.5*D_FS[Index]);
            ap[Index] *= C_over[Index];
            ae[Index] *= C_over[Index];
            aw[Index] *= C_over[Index];
            an[Index] *= C_over[Index];
            as[Index] *= C_over[Index];
            b[Index]  *= C_over[Index];

            // go to next cell in x-direc
            Index++;
        }
        // go to next line in y-direction to the first cell
        Index += 2;
    }
    ExtrudC();
    // solve equations using conjugate gradient methods (GMRES)
    NC = CGM(ae, aw, an, as, ap, b, u, Nx, Ny, maxit, GMRES);
    //NC = CGM(ae, aw, an, as, ap, b, u, Nx, Ny, maxit, CGHS);
    
    // restoring of all calculated values as new values
    for (i = 0; i < DimInfo->Nxy; i++) {
        CellInfo[i].FracSol = FracSolNew[i];
        FL[i] = 1.-FracSolNew[i];
        CellInfo[i].Conc = u[i];
        cout<<"c["<<i<<"]="<<CellInfo[i].Conc << endl;
        CellInfo[i].CG = CellInfo[i].Conc * (1. - FracSolNew[i]*p);
    }
    return 0;
}


             

            

            

    
        
            
                
            
            
            
                
    
    
            
    
        


















    
