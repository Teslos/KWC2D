#include "sim2d.h"
#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;
#define M_PI 3.14

Csim2d::Csim2d()
{
    CellInfo = new Cell();
    TimeInfo = new CTime();
    DimInfo  = new Dimension();
}

void 
Csim2d::SetPar(int R, int x0, int y0, float theta0)
{
    int Index, i, j, k;
    float Summe, r, R2, r2, fsini, ftheta;
    R2 = (float)(R * R);
    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            r2  = ((float)(i - y0) + .5) * ((float)(i - y0) + .5);
            r2 += ((float)(j - x0) + .5) * ((float)(j - x0) + .5);
            r   = sqrt(r2);
			if (r < ((float)(R) + Delta / LengthX)) {
                Index  = DimInfo->Nx + 1;
                Index += j + i * (DimInfo->Nx);
                fsini  = r - (float)(R);
                //fsini *= 3. * LengthX / Delta;
				fsini *= LengthX / Delta;
                fsini  = tanh(fsini);
                fsini *= -0.5;
                fsini += 0.5;
                if (fsini > .99)
                    fsini = .99;
                CellInfo[Index].FracSol = fsini;
				FL[Index] = 1. - fsini;
            }
            // XXX set the orientation (theta) to a particle
            if (r < ((float)(R) + epsilonL / LengthX)) {
                Index = DimInfo->Nx + 1;
                Index += j + i * (DimInfo->Nx);
                ftheta = r - (float) (R);
                //ftheta *= 3. * LengthX / epsilonL;
				ftheta *= LengthX / epsilonL;
                ftheta = tanh(ftheta);
                ftheta *= -0.5*theta0;
                ftheta += 0.5*theta0;
                if (ftheta > .99*theta0)
                    ftheta = .99*theta0;
				//printf("Set particle index: %i, val: %lf, eps: %lf\n", Index, ftheta, epsilonL/LengthX);
                CellInfo[Index].Theta = ftheta;
                Theta[Index] = ftheta;
            }
        }
    }
}

int
Csim2d::InitBound(void)
{
    int sp;      // cell index
    int i,j,r,x,y,Xmax_xPos = 0, Ymax_xPos = 0, ij = 0;
    double FSR;  // sum of rows for frac. solid

    // call routine to set the round particles
    for (i = 1; i < Radius[0] + 1; i++) {
        SetPar(Radius[i], x0[i], y0[i], theta0[i]);
    }
    // planar front
    sp = DimInfo->Nx + 1;
    FS = 0.;
    for (j = 0; j < numY; j++) {
        FSR = 0.;
        for (i = 0; i < numX; i++) {
            CellInfo[sp].P = P_0;
            CellInfo[sp].P += (k - .5 * numX) * P_G;
        
            FSR += CellInfo[sp].FracSol;
            sp++;
        }
        FS += FSR;
        sp += 2;
    }
    FS /= (double)(numX * numY);

    // sets the starting concentration Lever/Scheil
    sp = DimInfo->Nx + 1;
    for (j = 0; j < numY; j++) {
        for (i = 0; i < numX; i++) {
            CellInfo[sp].Conc = C_0;
            CellInfo[sp].ConcS = CellInfo[sp].Conc * k;
            CellInfo[sp].CG = CellInfo[sp].Conc * FL[sp] + 
                CellInfo[sp].ConcS * CellInfo[sp].FracSol;
            CellInfo[sp].P = (k/(float)(numX) - .5) * P_G;
            if (CellInfo[sp].FracSol > 0.5) {
                Ymax_xPos = (Ymax_xPos > j) ? Ymax_xPos : j;
                Xmax_xPos = (Xmax_xPos > i) ? Xmax_xPos : i;
            }
            sp++;
        }
        sp += 2;
    }
    Xmax_Pos = Xmax_xPos;
    Ymax_Pos = Ymax_xPos;

    // and now set one more time particles to refresh concentration
    for (i = 1; i < Radius[0] + 1; i++) {
        SetPar(Radius[i], x0[i], y0[i], theta0[i]);
    }
	//exit(0);

    // XXX and finally setting of average values of fl at cell-boundries
    for (i = 0; i < DimInfo->Nxy; i++) {
        FLmx[i] = .5 * (FL[i] + FL[i-1]);
    }
    
    return 0;
}
    
                
    

int
Csim2d::InitRegion(istream &input)
{
    int i,j;
    float H[101];
    
    cout << "#####################" << endl;
    cout << "#" << endl;
    cout << "# Region informations" << endl;
    cout << "#" << endl;
    cout << "#####################" << endl;
    
    cout << "# How big is cell ? [cm] (float)" << endl;
    cout << "# for x ? " << endl;
    input >> LengthX;
    cout << "# for y ? " << endl;
    input >> LengthY;
    cout << "# for z ? " << endl;
    input >> LengthZ;

    cout << "# Number of cells...(int)" << endl;
    cout << "# in x direction ?" << endl;
    input >> numX;
    cout << "#" << numX << endl;
    cout << "# in y direction ?" << endl;
    input >> numY;
    cout << "#" << numY << endl;
    cout << "# in z direction ?" << endl;
    input >> numZ;
    cout << "#" << numZ << endl;

    // informations for moving frame
    cout << "# Maximal height of solid (moving frame)" << endl;
    input >> maxPos;
    cout << "#" << maxPos << endl;

    // informations for concentrations and diffusions
    cout << "# Slope of Liquidus in Phase-Diagram" << endl;
    input >> ML;
    cout << "#" << ML << endl;
    cout << "# Global starting concentration C0 [%] (float)" << endl;
    input >> C_0;
    cout << "#" << C_0 << endl;
    cout << "# Segregration coefficient k (float)" << endl;
    input >> k;
    cout << "#" << k << endl;
    cout << "# Diffusion coefficient liquid [cm^2/s] (float)" << endl;
    input >> DC[1];
    cout << "#" << DC[1] << endl;
    cout << "# Diffusion coefficient solid [cm^2/s] (float)" << endl;
    input >> DC[0];
    cout << "#" << DC[0] << endl;
    cout << "# End criteria for Fraction Liquid vs 1 " << endl;
    input >> crit;
    cout << "#" << crit << endl;
    
    // conductivity RhoCp from matrix and particles will be calculated
    // on diffusion coefficient 
    cout << "# Conduction coefficient of matrix ? [W/cm/K]" << endl;
    input >> lambda[1];
    cout << "#" << lambda[1] << endl;
    cout << "# Specific heat (rho_cp) of matrix ? [J/ccm/K]" << endl;
    input >> RhoCp[1];
    cout << "#" << RhoCp[1] << endl;
    cout << "# Conduction coefficient of particle ? [W/cm/K]" << endl;
    input >> lambda[0];
    cout << "#" << lambda[0] << endl;
    cout << "# Specific heat (rho_cp) of particle ? [J/ccm/K]" << endl;
    input >> RhoCp[0];
    cout << "#" << RhoCp[0] << endl;
    cout << "# Latent heat of matrix ? [J/ccm]" << endl;
    input >> L;
    cout << "#" << L << endl;
    cout << "# Latent heat of particle is set to 0" << endl;
    cout << "#" << endl;
    cout << "# Melting temperature of matrix [K] (float) " << endl;
    input >> Tmelt;
    cout << "#" << Tmelt << endl;
    cout << "# TPoint ( temperature derivation per time ) [K/s] (float)" << endl;
    input >> TPoint;
    cout << "#" << TPoint << endl;
    cout << "# Starting temperature at down border [K] (float)" << endl;
    input >> T_0;
    cout << "#" << T_0 << endl;
    cout << "# Temperature gradient over calculating region [K/cm] (float)" << endl;
    input >> TGrad;
    cout << "#" << TGrad << endl;
    cout << "# Delta - iface thickness [cm] (float)" << endl;
    input >> Delta;
    cout << "#" << Delta << endl;
    cout << "# Mobility of phase border [cm/s/K] (float)" << endl;
    input >> mobil;
    cout << "#" << mobil << endl;
    cout << "# Gibbs-Thomson Coef. [K cm] (float) " << endl;
    input >> Gamma;
    cout << "#" << Gamma << endl;
    cout << "# Kinetic anisotropy (float) " << endl;
    input >> epskin;
    cout << "#" << epskin << endl;
    cout << "# Border anisotropy (float) " << endl;
    input >> epssigma;
    cout << "#" << epssigma << endl;
    cout << "# Prefactor Anti-trapping (float)" << endl;
    input >> ATC;
    cout << "#" << ATC << endl;

    // Additional parameters for KWC model
    cout << "# Phase field gradient energy coefficient (float)" << endl;
    input >> alfa;
    cout << "#" << alfa << endl;
    cout << "# Phase field misorientation energy coefficient (float)" << endl;
    input >> s;
    cout << "#" << s << endl;
    cout << "# Characteristic length for grain boundary Epsilon(float)" << endl;
    input >> epsilonL;
    cout << "#" << epsilonL << endl;
    cout << "# Phase transient coefficient (float)" << endl;
    input >> tauPhi;
    cout << "#" << tauPhi << endl;
    cout << "# Orientation transient coefficient (float)" << endl;
    input >> tauTheta;
    cout << "#" << tauTheta << endl;
    cout << "# Undercooling parameter in M(phi) function (float)" << endl;
    input >> sigma;
    cout << "#" << sigma << endl;
    cout << "# scaling coefficient in P(grad phi) function Beta (float)" << endl;
    input >> beta;
    cout << "#" << beta << endl;
    cout << "# scaling mobility in P(grad phi) function Mu (float)" << endl;
    input >> mu;
    cout << "#" << mu << endl;
    cout << "# infinity diffusive constant (float)" << endl;
    input >> gamma;
    cout << "#" << gamma << endl;

    // particle informations and statistics of starting conditions.
    // number of particles is saved in 0 component.

    cout << "# Number of particles (< 100)" << endl;
    input >> Radius[0];
    cout << "#" << Radius[0] << endl;
    cout << "# Now " << Radius[0] << " x-coordinates and y-coordinates of particles" << endl;
    Curv = 0.;
    Curv1 = 0.;

    for (i = 1; i < Radius[0] + 1; i++) {
        cout << "#" << i << "Particle (Radius, x_M, y_M, theta_M)" << endl;
        input >> Radius[i] >> x0[i] >> y0[i] >> theta0[i];
        cout << "#" << Radius[i] <<"    "<< x0[i] 
            <<"    "<< y0[i] << "     " << theta0[i] << endl;
        Curv += (float) Radius[i];
        Curv1 += 1. / (float)(Radius[i]);
    }
    Curv /= (float)(Radius[0]);
    Curv1 /= LengthX;
    Curv1 /= (float)(Radius[0]);

    cout << "# Normal pressure (float) [N/m^2]" << endl;
    input >> P_0;
    // calculate from N/m^2 to g/cm/s^2
    P_0 *= 10.;
    cout << "#" << P_0 << endl;
    cout << "# Pressure gradient per cm (float) [N/m^2/cm]" << endl;
    input >> P_G;
    cout << "#" << P_G << endl;
    cout << "# Density in the liquid (float) [g/cm^3]" << endl;
    input >> Rho;
    cout << "#" << Rho << endl;
    cout << "# Variation of density with concentration (float)" << endl;
    input >> d_Rho;
    cout << "#" << d_Rho << endl;
    cout << "# Viscosity of melt (float) [cm^2/s]" << endl;
    input >> my;
    cout << "#" << my << endl;
    cout << "# velocity inside the melt (float)" << endl;
    input >> V_0;
    cout << "#" << V_0 << endl;
    cout << "# weight factor for Kozeny-Karman (float) " << endl;
    input >> KK;
    cout << "#" << KK << endl;
    cout << "# Accepted relative error in C-estimation (float) " << endl;
    input >> d_C;
    cout << "#" << d_C << endl;
    cout << "# Accepted relative error in V-estimation (float) " << endl;
    input >> Div;
    cout << "#" << Div << endl;

    // set the values
    AG = 0.;
    FS = 0.;
    NC = 0;
    V_control = 0.00001;
    n_move = 0;
    Ymax_Pos = 0;
    Xmax_Pos = 0;

    D[0] = lambda[0] / RhoCp[0];
    D[1] = lambda[1] / RhoCp[1];

    DR[0] = 1. / D[0];
    DR[1] = 1. / D[1];

    TGrad *= LengthX;
    T_D = T_0;

    return 0;
}

int 
Csim2d::InitDim(void)
{
    //Dimension *DimInfo = new Dimension();

    DimInfo->Nx = numX + 2;
    DimInfo->Ny = numY + 2;
    DimInfo->Nxy = DimInfo->Nx*DimInfo->Ny;
    return 0;
}
/// 
/// Initialize time
///
int
Csim2d::InitTime(istream &input)
{
    // initialize time informations
    // ask for time informations
    cout << "##########################" << endl;
    cout << "#" << endl;
    cout << "# Time informations" << endl;
    cout << "#" << endl;
    cout << "##########################" << endl;

    cout << "# Time step width ? [Seconds] float" << endl;
    input >> TimeInfo->tWidth;
    cout << "#" << TimeInfo->tWidth << endl;
    cout << "# How many time steps ? (int)" << endl;
    input >> TimeInfo->tSteps;
    cout << "#" << TimeInfo->tSteps <<endl;
   
    KK *= my;
    KK /= (Delta * Delta);
    KK *= 6. * 6. * 2.;
    KK *= TimeInfo->tWidth;
    
    // everything is ok
    return 0;
}
/// Initialize cell
///
///
int
Csim2d::InitCell(void)
{
    int i, j;

    CellInfo = new Cell[DimInfo->Nxy];
    if (CellInfo == NULL) {
        cerr << "Error allocating memory for: CellInfo" << endl;
        return -1;
    }
    // initialize cell informations with starting values
    for (j = 0; j < DimInfo->Nxy; j++) {
        CellInfo[j].Conc = C_0;
        CellInfo[j].ConcS = C_0 * k;
        CellInfo[j].CG   = 0.;
        CellInfo[j].FracSol = 0.;
        CellInfo[j].Theta = 0.;
        CellInfo[j].VX = 0.;
        CellInfo[j].VY = 0.;
        CellInfo[j].P  = P_0;
    }

    return 0;
}

///
/// initialize temperature
///
int
Csim2d::InitTemperature(void)
{
    int i,j;

    FracSolNew = new float[DimInfo->Nxy];
    if (FracSolNew == NULL) {
        cerr << "Error allocating memory for: FracSolNew" << endl;
        return -1;
    }

    D_FS = new float[DimInfo->Nxy];
    if (D_FS == NULL) {
        cerr << "Error allocating memory for: D_FS" << endl;
        return -1;
    }

    FL = new float[DimInfo->Nxy];
    if (FL == NULL) {
        cerr << "Error allocating memory for: FL Fraction Liquid" << endl;
        return -1;
    }

    FLmx = new float[DimInfo->Nxy];
    if (FLmx == NULL) {
        cerr << "Error allocating memory for: Fraction Liquid Middle value in x direction" << endl;
        return -1;
    }
    
    FLmy = new float[DimInfo->Nxy];
    if (FLmy == NULL) {
        cerr << "Error allocating memory for: Fraction Liquid Middle value in y direction" << endl;
        return -1;
    }

    FSmx = new float[DimInfo->Nxy];
    if (FSmx == NULL) {
        cerr << "Error allocating memory for: Fraction Solid Middle value in x direction" << endl;
        return -1;
    }
    
    FSmy = new float[DimInfo->Nxy];
    if (FSmy == NULL) {
        cerr << "Error allocating memory for: Fraction Solid Middle value in y direction" << endl;
        return -1;
    }
    
    Tmx = new double[DimInfo->Nxy];
    if (Tmx == NULL) {
        cerr << "Error allocating memory for: Gradient Theta field in the middle of x direction" << endl;
        return -1;
    }

    Tmy = new double[DimInfo->Nxy];
    if (Tmy == NULL) {
        cerr << "Error allocating memory for: Gradient Theta field in the middle of y direction" << endl;
        return -1;
    }
    
    Cmx = new double[DimInfo->Nxy];
    if (Cmx == NULL) {
        cerr << "Error allocating memory for: Correction theta in the middle of x direction" << endl;
        return -1;
    }

    Cmy = new double[DimInfo->Nxy];
    if (Cmy == NULL) {
        cerr << "Error allocating memory for: Correction theta in the middle of y direction" << endl;
        return -1;
    }

    Theta = new double[DimInfo->Nxy];
    if (Theta == NULL) {
        cerr << "Error allocating memory for: Theta field" << endl;
        return -1;
    }

    theta_over = new float[DimInfo->Nxy];
    if (theta_over == NULL) {
        cerr << "Error allocating memory for: Theta_over field" << endl;
        return -1;
    }

    Dmx = new double[DimInfo->Nxy];
    if (Dmx == NULL) {
        cerr << "Error allocating memory for: Dmx diffusion constant" << endl;
        return -1;
    }

    Dmy = new double[DimInfo->Nxy];
    if (Dmy == NULL) {
        cerr << "Error allocating memory for: Dmy diffusion constant" << endl;
        return -1;
    }

    RX = new float[DimInfo->Nxy];
    if (RX == NULL) {
        cerr << "Error allocating memory for: Friction coefficient in x direction" << endl;
        return -1;
    }
    
    RY = new float[DimInfo->Nxy];
    if (RY == NULL) {
        cerr << "Error allocating memory for: Friction coefficient in y direction" << endl;
        return -1;
    }

    RP = new float[DimInfo->Nxy];
    if (RP == NULL) {
        cerr << "Error allocating memory for: Friction coefficient over RP" << endl;
        return -1;
    }

    RPx = new float[DimInfo->Nxy];
    if (RPx == NULL) {
        cerr << "Error allocating memory for: Friction coefficient over RPx" << endl;
        return -1;
    }
    
    RPy = new float[DimInfo->Nxy];
    if (RPy == NULL) {
        cerr << "Error allocating memory for: Friction coefficient over RPy" << endl;
        return -1;
    }
    
    ConcNew = new float[DimInfo->Nxy];
    if (ConcNew == NULL) {
        cerr << "Error allocating memory for: Concentration new" << endl;
        return -1;
    }

    Dbar = new float[DimInfo->Nxy];
    if (Dbar == NULL) {
        cerr << "Error allocating memory for: Dbar" << endl;
        return -1;
    }

    oneminus = new float[DimInfo->Nxy];
    if (oneminus == NULL) {
        cerr << "Error allocating memory for: oneminus" << endl;
        return -1;
    }

    C_over = new float[DimInfo->Nxy];
    if (C_over == NULL) {
        cerr << "Error allocating memory for: C_over" << endl;
        return -1;
    }

    KKonv = new float[DimInfo->Nxy];
    if (KKonv == NULL) {
        cerr << "Error allocating memory for: KConv" << endl;
        return -1;
    }

    VXZ = new float[DimInfo->Nxy];
    if (VXZ == NULL) {
        cerr << "Error allocating memory for: VXZ" << endl;
        return -1;
    }

    VxNew = new float[DimInfo->Nxy];
    if (VxNew == NULL) {
        cerr << "Error allocating memory for: VxNew" << endl;
        return -1;
    }

    VyNew = new float[DimInfo->Nxy];
    if (VyNew == NULL) {
        cerr << "Error allocating memory for: VyNew" << endl;
        return -1;
    }

    DIV = new float[DimInfo->Nxy];
    if (DIV == NULL) {
        cerr << "Error allocating memory for: Divergence" << endl;
        return -1;
    }

    PNew = new float[DimInfo->Nxy];
    if (PNew == NULL) {
        cerr << "Error allocating memory for: Pressure new" << endl;
        return -1;
    }

    // initialize everything with 0.
    for (j = 0; j < DimInfo->Nxy; j++) {
        FracSolNew[j] = 0.;
        FL[j]         = 1.;
        FLmx[j]       = 1.;
        FLmy[j]       = 1.;
        FSmx[j]       = 0.;
        FSmy[j]       = 0.;
        Tmx[j]        = 0.;
        Tmy[j]        = 0.;
        Cmx[j]        = 0.;
        Cmy[j]        = 0.;
        Theta[j]      = M_PI;
        theta_over[j] = 0.;
        Dmx[j]        = 0.;
        Dmy[j]        = 0.;
        RX[j]         = 0.;
        RY[j]         = 0.;
        RP[j]         = 0.;
        RPx[j]        = 0.;
        RPy[j]        = 0.;
        VxNew[j]      = 0.;
        VyNew[j]      = 0.;
        DIV[j]        = 0.;
        ConcNew[j]    = 0.;
        D_FS[j]       = 0.;
        oneminus[j]   = 0.;
        Dbar[j]       = 0.;
        KKonv[j]      = 0.;
        C_over[j]     = 0.;
    }

    return 0;
}

//
// Returns true if stability condition
// is fullfilled
//
bool
Csim2d::StabilCond()
{
    float XQuadrat = LengthX * LengthX;
    return TimeInfo->tWidth < tauTheta / (alfa * alfa) * XQuadrat * 0.25; 
}
