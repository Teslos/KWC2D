#include<iostream>
#include "sim2d.h"
using namespace std;

void Csim2d::Extruder(void)
{
    CpX(1,0);  // go down (negative y-direction)
    CpX(DimInfo->Ny - 2, DimInfo->Ny - 1); // go up (positive y-direction)
    CpY(1,0);  // go right (positive x-direction)
    CpY(DimInfo->Nx - 2, DimInfo->Nx - 1); // go left (negative x-direction)
}

void Csim2d::CpX(int ysrc, int ydst)
{
    int i,     // counter in x-direction
        sp,    // index the source cell
        dp;    // destination cell
    // start points
    sp = ysrc * DimInfo->Nx;
    dp = ydst * DimInfo->Nx;
    // loop over all elements in x-line
    for (i = 0; i < DimInfo->Nx; i++) {
        CellInfo[dp].FracSol = CellInfo[sp].FracSol;
        FracSolNew[dp] = FracSolNew[sp];
        FL[dp] = FL[sp];
		// 24.06.14 added a new theta field which also goes into phase field equation
        //CellInfo[dp].Theta = CellInfo[sp].Theta;
        Theta[dp] = Theta[sp];
        // increase index
        sp++;
        dp++;
    }
}

void Csim2d::CpY(int xsrc, int xdst)
{
    int i,
        sp,
        dp;
    // start points
    sp = xsrc;
    dp = xdst;
    // loop over all elements in y-line
    for (i = 0; i < DimInfo->Ny; i++) {
        CellInfo[dp].FracSol = CellInfo[sp].FracSol;
        FracSolNew[dp] = FracSolNew[sp];
        FL[dp] = FL[sp];
		// 24.06.14 added a new theta field which also goes into phase field equation
        //CellInfo[dp].Theta = CellInfo[sp].Theta;
        Theta[dp] = Theta[sp];
        sp += DimInfo->Nx;
        dp += DimInfo->Nx;
    }
     
}

void Csim2d::cpFL(int xsrc, int xdst)
{
    int i,
        sp,
        dp;

    sp = xsrc;
    dp = xdst;

    for (i = 0; i < DimInfo->Ny; i++) {
        FLmy[dp] = FLmy[sp];
        RY[dp] = RY[sp];
        RPy[dp] = RPy[sp];
        sp += DimInfo->Nx;
        dp += DimInfo->Nx;
    }
}

void Csim2d::ExtrudV(void)
{
    CpXV(1, 0);
    CpXV(DimInfo->Ny - 2, DimInfo->Ny - 1);

    CpYV(1, 0);
    CpYV(DimInfo->Nx - 2, DimInfo->Nx - 1);
}

void Csim2d::CpXV(int ysrc, int ydst)
{
    int i,
        sp,
        dp;

    // calculate starting points
    sp = ysrc * DimInfo->Nx;
    dp = ydst * DimInfo->Nx;

    // loop over all elements in x-line
    for (i = 0; i < DimInfo->Nx; i++) {
        VxNew[dp] = VxNew[sp];
        VxNew[dp] = 0;
        sp++;
        dp++;
    }
}

void Csim2d::CpYV(int xsrc, int xdst)
{
    int i,
        sp,
        dp;

    sp = xsrc;
    dp = xdst;

    // loop over all elements in y-line
    for (i = 0; i < DimInfo->Ny; i++) {
        VxNew[dp] = 0.;
        VxNew[dp] = VxNew[sp];
        sp += DimInfo->Nx;
        dp += DimInfo->Nx;
    }
}

void Csim2d::CpXD(int ysrc, int ydst)
{
    int i,
        sp,
        dp;

    sp = ysrc * DimInfo->Nx;
    dp = ydst * DimInfo->Nx;

    for (i = 0; i < DimInfo->Nx; i++) {
        DIV[dp] = DIV[sp];
        sp++;
        dp++;
    }
}

void Csim2d::CpYD(int xsrc, int xdst)
{
    int i,
        sp,
        dp;

    sp = xsrc;
    dp = xdst;

    for (i = 0; i < DimInfo->Ny; i++) {
        DIV[dp] = DIV[sp];
        sp += DimInfo->Nx;
        dp += DimInfo->Nx;
    }
}

void Csim2d::ExtrudP(void)
{
    CpXP(1, 0);
    CpXP(DimInfo->Ny - 2, DimInfo->Ny - 1);

    CpYP(1, 0);
    CpYP(DimInfo->Nx - 2, DimInfo->Nx - 1);
}
void Csim2d::CpXP(int ysrc, int ydst)
{
    int i,
        sp,
        dp;

    sp = ysrc * DimInfo->Nx;
    dp = ydst * DimInfo->Nx;

    // loop over all elemente in x line
    for (i = 0; i < DimInfo->Nx; i++) {
        CellInfo[dp].P = CellInfo[sp].P;
        sp++;
        dp++;
    }
}

void Csim2d::CpYP(int xsrc, int xdst)
{
    int i,
        sp,
        dp;
    
    for (i = 0; i < DimInfo->Ny; i++) {
        // set pressure at boundary
        if (xdst == 0) 
            CellInfo[dp].P = CellInfo[sp].P - P_G;
        else
            CellInfo[dp].P = CellInfo[sp].P + P_G;
        sp += DimInfo->Nx;
        dp += DimInfo->Nx;
    }
}

void Csim2d::ExtrudPn(void)
{
    CpXPn(1, 0);
    CpXPn(DimInfo->Ny - 2, DimInfo->Ny - 1);

    CpYPn(1, 0);
    CpYPn(DimInfo->Nx - 2, DimInfo->Nx - 1);
    
}

void Csim2d::CpXPn(int ysrc, int ydst)
{
    int i,
        sp,
        dp;

    // calculate start points
    sp = ysrc * DimInfo->Nx;
    dp = ydst * DimInfo->Nx;

    // loop over all elements in x-line
    for (i = 0; i < DimInfo->Nx; i++) {
        PNew[dp] = PNew[sp];
        sp++;
        dp++;
    }
}

void Csim2d::CpYPn(int xsrc, int xdst)
{
    int i,
        sp,
        dp;

    sp = xsrc;
    dp = xdst;

    // loop over all elements in y-line
    for (i = 0; i < DimInfo->Ny; i++) {
        PNew[dp] = PNew[sp];
        sp += DimInfo->Nx;
        dp += DimInfo->Nx;
    }
}

void Csim2d::ExtrudC(void)
{
    CpXC(1, 0);
    CpXC1(DimInfo->Ny - 2, DimInfo->Ny - 1);

    CpYC(1, 0);
    CpYC(DimInfo->Nx - 2, DimInfo->Nx - 1);
}

void Csim2d::CpXC(int ysrc, int ydst)
{
    int i,
        sp,
        dp;

    sp = ysrc * DimInfo->Nx;
    dp = ydst * DimInfo->Nx;

    for (i = 0; i < DimInfo->Nx; i++) {
        u[dp] = u[sp];
        sp++;
        dp++;
    }
}

void Csim2d::CpXC1(int ysrc, int ydst)
{
    int i,
        sp,
        dp;

    sp = ysrc * DimInfo->Nx;
    dp = ydst * DimInfo->Nx;

    // loop over all elements in x-line
    for (i = 0; i < DimInfo->Nx; i++) {
        u[dp] = C_0;
        sp++;
        dp++;
    }
}

void Csim2d::CpYC(int xsrc, int xdst)
{
    int i,
        sp,
        dp;

    sp = xsrc;
    dp = xdst;
    
    for (i = 0; i < DimInfo->Ny; i++) {
        u[dp] = u[sp];
        sp += DimInfo->Nx;
        dp += DimInfo->Nx;
    }
}

void Csim2d::ExtrudDFS(void)
{
    CpXDFS(1, 0);
    CpXDFS(DimInfo->Ny - 2, DimInfo->Ny - 1);

    CpYDFS(1, 0);
    CpYDFS(DimInfo->Nx - 2, DimInfo->Nx - 1);
}

void Csim2d::CpXDFS(int ysrc, int ydst)
{
    int i,
        sp,
        dp;

    sp = ysrc * DimInfo->Nx;
    dp = ydst * DimInfo->Nx;

    // loop over all elements in x-line
    for (i = 0; i < DimInfo->Nx; i++) {
        D_FS[dp] = D_FS[sp];
        sp++;
        dp++;
    }
    
}

void Csim2d::CpYDFS(int xsrc, int xdst)
{
    int i,
        sp,
        dp;

    sp = xsrc;
    dp = xdst;

    for (i = 0; i < DimInfo->Ny; i++) {
        D_FS[dp] = D_FS[sp];
        sp += DimInfo->Nx;
        dp += DimInfo->Nx;
    }
}

void Csim2d::ExtrudTh(void)
{
    CpThX(1,0);  // go down (negative y-direction)
    CpThX(DimInfo->Ny - 2, DimInfo->Ny - 1); // go up (positive y-direction)
    CpThY(1,0);  // go right (positive x-direction)
    CpThY(DimInfo->Nx - 2, DimInfo->Nx - 1); // go left (negative x-direction)
}

void Csim2d::CpThX(int ysrc, int ydst)
{
    int i,     // counter in x-direction
        sp,    // index the source cell
        dp;    // destination cell
    // start points
    sp = ysrc * DimInfo->Nx;
    dp = ydst * DimInfo->Nx;
    // loop over all elements in x-line
    for (i = 0; i < DimInfo->Nx; i++) {
        CellInfo[dp].Theta = CellInfo[sp].Theta;
        Theta[dp] = Theta[sp];
		// added 24.06.14 liquid fraction necessary to calculate 
		// theta field.
		FL[dp] = FL[sp];
        // increase index
        sp++;
        dp++;
    }
}

void Csim2d::CpThY(int xsrc, int xdst)
{
    int i,
        sp,
        dp;
    // start points
    sp = xsrc;
    dp = xdst;
    // loop over all elements in y-line
    for (i = 0; i < DimInfo->Ny; i++) {
        CellInfo[dp].Theta = CellInfo[sp].Theta;
        Theta[dp] = Theta[sp];
		FL[dp] = FL[sp];
        sp += DimInfo->Nx;
        dp += DimInfo->Nx;
    }
     
}

