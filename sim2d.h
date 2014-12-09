#ifndef __SIM2D__
#define __SIM2D__
#include<iostream>
using namespace std;

class Csim2d {
    float LengthX;    // dimension of one cell in microns
    float LengthY;    
    float LengthZ;

    int   numX;    // number of cells in x,y,z direction 
    int   numY;
    int   numZ;

    int   maxPos;  // maximal pos in y for solid for moving frame
    int   Xmax_Pos;  // maximal x solid position
    int   Ymax_Pos;  // maximal y solid position

    int   n_move;  // number of cell which droped out of frame.
    int   NC;      // iteration number for implicit concentration solver.
    int   NP;      // iteration number for implicit pressure solver.
    int   NV;      // iteration number for implicit velocity solver.
    int   ND;      // iteration number for divergence in velocity solver.
    int   NTheta;  // iteration number for implicit theta field solver.
	int   NF;      // iteration number for semi-implicit phase field solver.

    float NR0;     // relative number of particles.
    float NR;      // relative number of particles on start.
    float N0;
    float ML;      // slope of liquidus line.
    float DC[2];   // Diffusionconstant. concentration (Particles/M)
    float k;       // segregation coefficient.
    float C_0;     // global starting concentration.
    float C_control; // sum of C (controlling conservation of conc.)
    float C_liquid;  // sum of cl*fl/sum(fl).
    float C_error;   // found error for implicit concentration solver.
	float F_error;   // found error for semi-implicit phase solver.
	float T_error;   // found error for implicit theta solver.
    float P_error;   // found error for implicit pressure solver.
    float V_error;   // found error for implicit velocity solver.
    float D_error;   // found error for divergence in V-solver.
    float d_C;       // allowed error for conc.solver.
	float d_T;       // allowed error for theta solver.
    float d_F;       // allowed error for phase-field solver.
    double FS;       // sum of p-field for control.
    float AG;        //
    float Curv;      // control value for middle curvature p-field.
    float Curv1;     // control value for middle curvature n-vector.
    float crit;      // criteria for FL->0
    float lambda[2]; // conduction in W/cm/K
    float RhoCp[2];  // specific heat in J/ccm/K
    float D[2];      // temperature diffusion coeff. for matrix/ball
    float DR[2];     // reciprocal value.
    float L;         // lattent heat J/ccm
    float Tmelt;     // melting temperature
    float TPoint;    // deltaT what is given away per time-step.
    float T_0;       // starting temperature on the down border.
    float T_D;       // temperature at down border.
    float TGrad;     // temperature gradient over the cell.
    float Utip;      // tip-unterkuehlung = (T*-Tm(c_0))/(1-k)mc_0
    float Otip;      // tip-overconcentration 

    float y_bef;
    float SI1;
    float Delta;
    float ATC;       // prefactor for ATC 
    float mobil;     // mobility of phase-field interface
    float Gamma;     // Gibbs-Thomson coeff.
    float epskin;    
    float epssigma;  // Amplitudes of growing interface
    float Rho;       // density of liquid
    float d_Rho;     // density difference liquid-solid
    float my;        // viscosity
    float P_0;       // initial pressure for NS
    float P_G;       // pressure gradient.
    float V_0;       // 
    float V_control; // control velocity
    float KK;        // Kozeny-Karman factor
    float Div;       // allowed value for divergence.
    float DivMax;    // maximal divergence
    int Radius[200]; // radius of particles.
    int x0[200];     // center of particles.
    int y0[200];      
    float theta0[200]; // initial orientation of grain. -pi <= theta0 < pi

    float alfa;      // gradient coefficient for phi field
    float s;         // gradient coefficient for theta field
    float epsilonL;  // characteristic length of grain boundary
    float tauPhi;    // time-relaxation parameter for phi field
    float tauTheta;  // time-relaxation parameter for theta field
    float sigma;     // undercooling parameter in m(phi) function
    float beta;      // scaling coefficient in P function.
    float mu;        // scaling mobility in P function.
    float gamma;     // gamma big number (infinite diffusitivity)

    double *ap; 
    double *ae; 
    double *aw;
    double *an;
    double *as;
    double *ane;
    double *anw; 
    double *ase;
    double *asw;
    double *b;
    double *u;

    // structures 
    struct CTime* TimeInfo;
    struct Cell* CellInfo;
    struct Dimension* DimInfo;
public:
    Csim2d();
    ~Csim2d();
    Cell *getCellInfo() { return CellInfo; }
    CTime *getTimeInfo() { return TimeInfo; }
    Dimension *getDimension() {return DimInfo; }
    // functions
    void resetv(void);
    int alloc_c(void);
    void dealloc_c(void);
    void SetPar(int, int, int, float);
    int InitBound(void);
    int InitRegion(istream&);
    int InitTime(istream&);
    int InitTemperature();
    int InitCell(void);
    int InitDim(void);
    bool StabilCond(void);
    // inline functions
    inline double Igamma(double);
    inline double Agamma(double);
    inline double P(double);
    inline double Subt(double);
	inline double Corr(double);
    // sets boundary conditions;
    // from extruder.cc
    void Extruder(void);
    void CpX(int ysrc, int ydst);
    void CpY(int xsrc, int xdst);
    void cpFL(int xsrc, int xdst);
    void ExtrudV(void);
    void CpXV(int ysrc, int ydst);
    void CpYV(int xsrc, int xdst);
    void ExtrudD(void);
    void CpXD(int ysrc, int ydst);
    void CpYD(int xsrc, int xdst);
    void ExtrudP(void);
    void CpXP(int ysrc, int ydst);
    void CpYP(int xsrc, int xdst);
    void ExtrudPn(void);
    void CpXPn(int ysrc, int ydst);
    void CpYPn(int xsrc, int xdst);
    void ExtrudC(void);
    void CpXC(int ysrc, int ydst);
    void CpXC1(int ysrc, int ydst);
    void CpYC(int xsrc, int xdst);
    void ExtrudDFS(void);
    void CpXDFS(int ysrc, int ydst);
    void CpYDFS(int xsrc, int xdst);
    void ExtrudTh(void);
    void CpThX(int ysrc, int ydst);
    void CpThY(int xsrc, int ydst);
    // solver from solve2d.cc
    int solve2d();
    // outcolor.cc
    int OutContr(int);
    int OutTmpp(int);
    int OutSimgeo(void);
#ifdef __TECPLOT__
    int OutTec(int);
#endif
    int OutMathematica(int);
    void OutConv(void);
    int OutCoeff(void);

    // public fields
    float *FracSolNew;
    float *FSold;
    float *ConcNew;
    float *FL;
    float *FLmx;
    float *FLmy;
    float *FSmx;
    float *FSmy;
    double *Theta;
	double *T_over;
	double *Dmbar;
	double *Pbar;
    float *theta_over;
    double *Tmx;
    double *Tmy;
    double *Cmx;
    double *Cmy;
    double *Dmx;
    double *Dmy;
    float *RX;
    float *RY;
    float *RP;
    float *RPx;
    float *RPy;
    float *PNew;
    float *VxNew;
    float *VyNew;
    float *VXZ;
    float *VYZ;
    float *DIV;
    float *KKonv;
    float *Dbar;
    float *oneminus;
    float *C_over;
    float *D_FS;
	float ALPHA;
};

struct CTime {
    float tWidth;
    int   tSteps;
    float time;
};

struct Cell {
    float Conc;
    float ConcS;
    float CG;
    float T;
    float P;
    float VX;
    float VY;
    float FracSol;
    float Theta;
    int   state;
    int   material;
};

struct Dimension{
    int Nx;
    int Ny;
    int Nz;
    int Nxy;
    int Nyz;
    int Nxz;
    int Nxyz;
};

extern "C" void slapsol_(double *ae, double *aw, double *an, double *as, double *ap,
        double *b, double *u, int *ni, int *nj, int *maxit,char *method);

//Csim2d *csim2d;
// debug stuff
int DNAG (double *ae, double *aw, double *an, double *as, double *ap,
        double *b, double *u, int ni, int nj, int maxit, int method);
int
NAGSOL(double *ae, double *aw, double *an, double *as, double *ap, 
        double *bm,double *u, int ni, int nj, char* enum_meth, char* enum_piv,
        char* enum_fact, int maxiter);

#endif
    
        
