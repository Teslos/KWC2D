#include <iostream>
#include <math.h>

#include "sim1d.h"

// solvers
#include "cghs.h"

int 
Csim1d::solve1d(void)
{
    // constants
    // calculate reciprocal squares
    //
    float XR = 1. / LengthX;
    float twoXR = 0.5 * XR;

    float XRQuadrat = XR * XR;

    int Nx = DimInfo->Nx;

    // reset all coefficients



