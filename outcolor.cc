#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include "sim2d.h"

#define SIMGEONAME "sim.geoF"
#define TMPPNAME "sim.tmpp"
#define FSPNAME "sim.fsp"
#define KONZNAME "sim.konz"
#define PNAME "sim.druck"
#define VXNAME "sim.vx"
#define VYNAME "sim.vy"
#define TMPPSAVE "sav.tmpp"
#define FSPSAVE "sav.fsp"
#define CONCSAVE "sav.konz"
#define PSAVE "sav.druck"
#define VXSAVE "sav.vx"
#define VYSAVE "sav.vy"

#define HANAME "dis"

using namespace std;
// helper functions
void printarray(ofstream& out, int numPoints, double *c);

int Csim2d::OutContr(int tstep)
{
    char *preafix, // preafix for output file
        *suffix,   // suffix for output file
        *outname;  // outname = preafix.suffix
    int i,         // counter for two dimensions
        j,
        nxpoint,   // counter for statistics
        k=0, l=0, Nx, sp;

    FILE *fpout;   // file for distribution function

    float time,    // time in sec
          pk,      // 1-k
          K=0,     // Curvature
          R=0,     // Curvature radius
          GradX=0.,
          GradY=0.,
          N=0.,
          K_mp=0.,
          K_mN=0.,
          R_average=0.,
          H_sum=0.,
          Surf=0.,
          C_all=0.,
          C_lm=0.,
          DivMax=0.,
          DivAv=0.,
          v_superx=0.,
          v_supery=0.,
          v_supery1=0.,
          v_super1=0.,   // superficial velocity in one cell
          v_superxy=0.,  // superficial velocity in x+y direction
          d_0,sigma, Pe, k0,k1,kmean,kmean1,  // Curvature values in top dendrites
          phix0, phiy0,  // gradients at top dendrite
          phix1, phiy1,  // gradients at top dendrite
          x_max, y_max,  // positions of top-dendrite
          y_0, laplace, vel, vel1;
    float *Hau;          // vector for distribution-function
    int NumPoints,       // number of points to write
        RecLen = 0;      // length of data
    FILE *fpT, *fpF, *fpC, *fpP, *fpVX, *fpVY;
    float fs, temp, konz, p, vx, vy;

    time = tstep * (TimeInfo->tWidth);
    Nx = DimInfo->Nx;
    NumPoints = numX * numY;
    pk = 1. - k;
    // allocate for distribution function
    Hau = new float[DimInfo->Nx];
    if (Hau == NULL) {
        cout << "Error in allocation of memory for Hau" << endl;
        exit(-1);
    }

    char *praefix = new char[6];
    suffix = new char[25];

    outname = new char[32];
    // initialize control values
    DivAv = 0.;
    DivMax = 0.;
    k = 0;
    K_mp = 0.;
    R_average = 0.;

    // loop for checking the complete concentration, surface 
    sp = DimInfo->Nx + 1;
    // initialize
    AG = 0;
    C_control = 0.;
    C_liquid = 0.;
    k = 0;
    for (i = 0; i < numY; i++) {
        C_all = 0.;
        C_lm = 0.;
        Surf = 0.;
        for (j = 0; j < numX; j++) {
            // surface stress
            GradX = CellInfo[sp + 1].FracSol - CellInfo[sp - 1].FracSol;
            GradY = CellInfo[sp + Nx].FracSol - CellInfo[sp - Nx].FracSol;
            Surf += sqrt(GradX * GradX + GradY * GradY);

            // concentration control in all and liquid
            C_all += CellInfo[sp].CG;
            if (FL[sp] > .5) {
                C_lm += CellInfo[sp].Conc;
                k++;
            }
            // maximal divergence error
            DivMax = (DivMax > fabs(DIV[sp])) ? DivMax : fabs(DIV[sp]);
            // average divergence error
            DivAv += fabs(DIV[sp]);
            sp++;
        }
        AG += Surf;
        C_control += C_all;
        C_liquid += C_lm;
        sp += 2;
    }
    DivAv /= numX;
    DivAv /= numY;
    DivAv /= (LengthX * LengthX * Rho);
    DivAv *= TimeInfo->tWidth;

    DivMax = (LengthX * LengthX * Rho);
    DivMax *= TimeInfo->tWidth;

    cout << "# Dav " << DivAv << " Dmax " << DivMax << " ,ymax " << Ymax_Pos
         << " xmax " << Xmax_Pos << endl;

    AG *= LengthX * .5;

    // save concentration control values
    C_control /= numX;
    C_control /= numY;

    C_liquid /= (double)k;

    // calculate velocity
    // calculate superfical velocity from volume flow in y 
    // along x line in the middle
    // startpoint left above
    sp = DimInfo->Nx + 1;
    v_supery = 0.;
    v_supery1 = 0.;
    sp += (int) (0.5 * numX) * DimInfo->Nx;

    for (i = 0; i < numY; i++) {
        v_supery += CellInfo[sp].VY;
        v_supery1 += abs(CellInfo[sp].VY);
        // jump to next line
        sp += 1;
    }
    v_supery /= numY;
    v_supery1 /= numY;
    cout << TimeInfo->time << " " << 10000 * v_supery << " ";
    // calculate superfical velocity from volume flow in
    // x direction along y middle line
    sp = DimInfo->Nx + 1;
    sp += (int)(.5 * numX);

    v_superx = 0.;
    
    for (i = 0; i < numY; i++) {
        v_superx += CellInfo[sp].VX;
        // jump to next line
        sp += DimInfo->Nx;
    }
    v_superx /= (float) (numY);
    cout << 10000 * v_superx << " ";
    
    // and now average over all cells
    sp = DimInfo->Nx + 1;
    v_superxy = 0.;

    for (i = 0; i < numY; i++) {
        v_super1 = 0;
        for (j = 0; j < numX; j++) {
            v_super1 += .5 * (CellInfo[sp].VX + CellInfo[sp + 1].VX);
            v_super1 += .5 * (CellInfo[sp].VY + CellInfo[sp + Nx].VY);

            sp++;
        }
        v_superxy += v_super1;
        sp += 2;
    }
    v_superxy /= numY;
    v_superxy /= numX;
    cout << " " << 10000 * v_superxy << " " << 10000 * V_control << " ";
    
    // calculation of highest dendrite
    y_0 = 0.;
    y_max = y_0;
    x_max = j;
    k0 = 0.;
    k1 = 0.;
    kmean = 0.;
    vel = 0.;
    i = Ymax_Pos;

    nxpoint = Xmax_Pos + 10;
    if (nxpoint > numX)
        nxpoint = numX;
    for (j = Xmax_Pos; j < nxpoint; j++) {
        sp = Nx * (i + 1) + j + 1;
        y_0 = tanhf(2. * CellInfo[sp].FracSol - 1.);
        y_0 += Ymax_Pos;

        // when the highest point is reached, stop the position
        // and calculate radius and velocity (two point interpolation)
        if (y_0 > y_max) {
            y_max = y_0;
            x_max = j;
            phix0 = CellInfo[sp + 1].FracSol - CellInfo[sp - 1].FracSol;
            phiy0 = CellInfo[sp + Nx].FracSol - CellInfo[sp - Nx].FracSol;

            k0  = CellInfo[sp + 1].FracSol + CellInfo[sp - 1].FracSol;
            k0 += CellInfo[sp + Nx].FracSol + CellInfo[sp - Nx].FracSol;
            k0 -= 4. * CellInfo[sp].FracSol;

            laplace  = CellInfo[sp + 1 + Nx].FracSol + CellInfo[sp - 1 - Nx].FracSol;
            laplace += CellInfo[sp + 1 - Nx].FracSol + CellInfo[sp - 1 + Nx].FracSol;
            laplace -= 4. * CellInfo[sp].FracSol;
            laplace  = 2 * k0 + 0.5 * laplace;
            laplace /= 3.;

            k0  = laplace - CellInfo[sp + Nx].FracSol - CellInfo[sp - Nx].FracSol;
            k0 += 2. * CellInfo[sp].FracSol;
            k0 /= sqrt(phix0 * phix0 + phiy0 * phiy0);
            k0 *= -2. / LengthX;

            phix1  = CellInfo[sp + 1 + Nx].FracSol - CellInfo[sp - 1 + Nx].FracSol;
            phiy1  = CellInfo[sp + Nx + Nx].FracSol - CellInfo[sp].FracSol;

            k1  = CellInfo[sp + 1 + Nx].FracSol + CellInfo[sp - 1 + Nx].FracSol;
            k1 += CellInfo[sp + Nx + Nx].FracSol + CellInfo[sp].FracSol;
            k1 -= 4. * CellInfo[sp + Nx].FracSol;

            laplace  = CellInfo[sp + 1 + Nx + Nx].FracSol + CellInfo[sp - 1].FracSol;
            laplace += CellInfo[sp + 1].FracSol + CellInfo[sp - 1 + Nx + Nx].FracSol;
            laplace -= 4. * CellInfo[sp + Nx].FracSol;

            k1  = 2. * k1 + 0.5 * laplace;
            k1 /= 3.;
            k1 -= CellInfo[sp + Nx + Nx].FracSol + CellInfo[sp].FracSol;
            k1 += 2. * CellInfo[sp + Nx].FracSol;
            k1 /= sqrt(phix1 * phix1 + phiy1 * phiy1);
            k1 *= -2. / LengthX;

            // now interpolate dependent on vicinity to real y_0
            kmean  = k0 * (Ymax_Pos + 1. - y_0);
            kmean += k1 * (y_0 - Ymax_Pos);

            vel  = -D_FS[sp] / pk;
            vel /= phiy0;
            vel *= 2. * LengthX / TimeInfo->tWidth;
            vel *= (1. + Ymax_Pos - y_0);
            
            vel1  = -D_FS[sp + Nx] / pk;
            vel1 /= phiy1;
            vel1 *= 2. * LengthX / TimeInfo->tWidth;
            vel1 *= (y_0 - Ymax_Pos);
            vel1 += vel;
        }
    }
    y_max += n_move;
    y_max *= LengthX;
    x_max *= LengthX;
    vel  = (y_max - y_bef);
    vel /= TimeInfo->tWidth * SI1;
    y_bef = y_max;
    // output positions 6-9
    cout << " " << 10000 * y_max << " " << 10000 * x_max << " "
        << 10000. * vel << " " << 10000. / kmean;
    d_0  = Gamma;
    d_0 /= (ML * pk * C_0 * (1 + pk * Utip));
    sigma  = 2. * d_0 * DC[1];
    sigma /= vel;
    sigma *= kmean * kmean;
    Pe  = 0.5 * vel / kmean;
    Pe /= DC[1];
    // on positions 10-13
    cout << sigma << " " << Pe << " " << Utip << " " 
        << " " << Otip << " ";
    // print all the rest
    cout << FS << " " << AG << " " << C_control 
        << " " << C_liquid;
    cout << T_D << " " << CellInfo[DimInfo->Nx + 1 + Xmax_Pos].CG;
    cout << " " << Ymax_Pos << " " << 10000. * v_supery1 << endl;

    // calculate record length
    RecLen = sizeof(float) + sizeof(int) + NumPoints * sizeof(float);

    fpT = fopen(TMPPSAVE, "wb+");
    fpF = fopen(FSPSAVE, "wb+");
    fpC = fopen(CONCSAVE, "wb+");
    fpP = fopen(PSAVE, "wb+");
    fpVX = fopen(VXSAVE, "wb+");
    fpVY = fopen(VYSAVE, "wb+");

    if (fpT == NULL || fpF == NULL || fpC == NULL || fpP == NULL
         || fpVX == NULL || fpVY == NULL) {
        cerr << "Result files couldn't be open" << endl;
        return -1;
    }

    fwrite(&RecLen, sizeof(int), 1, fpT);
    fwrite(&time, sizeof(float), 1, fpT);
    fwrite(&NumPoints, sizeof(int), 1, fpT);
    
    fwrite(&RecLen, sizeof(int), 1, fpF);
    fwrite(&time, sizeof(float), 1, fpF);
    fwrite(&NumPoints, sizeof(int), 1, fpF);

    fwrite(&RecLen, sizeof(int), 1, fpC);
    fwrite(&time, sizeof(float), 1, fpC);
    fwrite(&NumPoints, sizeof(int), 1, fpC);

    fwrite(&RecLen, sizeof(int), 1, fpP);
    fwrite(&time, sizeof(float), 1, fpP);
    fwrite(&NumPoints, sizeof(int), 1, fpP);
    
    fwrite(&RecLen, sizeof(int), 1, fpVX);
    fwrite(&time, sizeof(float), 1, fpVX);
    fwrite(&NumPoints, sizeof(int), 1, fpVX);
      
    fwrite(&RecLen, sizeof(int), 1, fpVY);
    fwrite(&time, sizeof(float), 1, fpVY);
    fwrite(&NumPoints, sizeof(int), 1, fpVY);
        
    sp = DimInfo->Nx + 1;

    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            fs = CellInfo[sp].FracSol;
            konz = CellInfo[sp].CG;
            p = CellInfo[sp].P;
            temp = DIV[sp];
            vx = 10000 * CellInfo[sp].VX;
            vy = 10000 * CellInfo[sp].VY;
            fwrite(&fs, sizeof(float), 1, fpF);
            fwrite(&konz, sizeof(float), 1, fpC);
            fwrite(&temp, sizeof(float), 1, fpT);
            fwrite(&vx, sizeof(float), 1, fpVX);
            fwrite(&vy, sizeof(float), 1, fpVY);
            fwrite(&p, sizeof(float), 1, fpP);
            sp++;
        }
        sp += 2;
    }

    fwrite(&RecLen, sizeof(int), 1, fpT);
    fwrite(&RecLen, sizeof(int), 1, fpF);
    fwrite(&RecLen, sizeof(int), 1, fpC);
    fwrite(&RecLen, sizeof(int), 1, fpP);
    fwrite(&RecLen, sizeof(int), 1, fpVX);
    fwrite(&RecLen, sizeof(int), 1, fpVY);
    fclose(fpT);
    fclose(fpF);
    fclose(fpC);
    fclose(fpP);
    fclose(fpVX);
    fclose(fpVY);

    return 0;
}

//
// Output the results in text file suitable for
// plotting in Mathematica.
//
int Csim2d::OutMathematica(int tstep)
{
	int i = 0,                  // counter in three dimensions
        j = 0, Nx = 0, sp = 0,  // cell index
        NumPoints = 0,          // number of points to write
        RecLen    = 0;          // record length of data
    ofstream *fpout;            // file for distribution function
    
    float K = 0, time, fs, temp, konz, p, vx, vy, R = 0., // radius of curvature
          R_average = 0.,  // average radius
          H_sum = 0.;      // sum of distribution

    float *Hau;
    time = tstep * (TimeInfo->tWidth);
    Nx = DimInfo->Nx;

    NumPoints = numX * numY;

	std::cout << "Writing in the Mathematica file" << std::endl;
    // allocate memory for distribution function
    Hau = new float[ DimInfo->Nx ];
    if (Hau == NULL) {
        cerr << "Error in allocating memory for Hau" << endl;
        exit(-1);
    }

    RecLen = sizeof(float) + sizeof(int) + NumPoints * sizeof(float);

	ofstream fpT("sim.t", ios_base::out | ios_base::app);
	ofstream fpF("sim.f", ios_base::out | ios_base::app);
	ofstream fpC("sim.k", ios_base::out | ios_base::app);
	ofstream fpP("sim.p", ios_base::out | ios_base::app);

	if ( !(fpT.is_open() && fpF.is_open() && fpC.is_open() && fpP.is_open()) ) {
        cerr << "Error in opening files" << endl;
        return -1;
    }
    
	fpT << RecLen << " ";
    fpT << time << " ";
    fpT << NumPoints << " ";
    
    fpF << RecLen << " ";
    fpF << time << " ";
    fpF << NumPoints << " ";
    
    fpC << RecLen << " ";
    fpC << time << " ";
    fpC << NumPoints << " ";

    fpP << RecLen << " ";
    fpP << time << " ";
    fpP << NumPoints << " ";

	sp = DimInfo->Nx + 1;

    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            fs = CellInfo[sp].FracSol;
            konz = CellInfo[sp].CG;
            p = CellInfo[sp].P;
            // XXX output theta instead of temperature
            temp = CellInfo[sp].Theta;
            vx = 5000. * (CellInfo[sp].VX + CellInfo[sp + 1].VX);
            vy = 5000. * (CellInfo[sp].VY + CellInfo[sp + Nx].VY);
            fpF << fs << " ";
            fpC << konz << " ";
            fpP << p << " ";
            fpT << temp << " ";
            sp++;
        }
        sp += 2;
    }
    fpT << RecLen << endl;
    fpF << RecLen << endl;
    fpC << RecLen << endl;
    fpP << RecLen << endl;
    fpT.close();
    fpF.close();
    fpC.close();
    fpP.close();
	return 0;
}

//
// Output all coefficients of the A matrix 
// including also boundary.
// This is used for debugging purpose.
//
int Csim2d::OutCoeff()
{
	int i = 0,                  // counter in three dimensions
        j = 0,                  // cell index
        NumPoints = 0;          // number of points to write
       
    
	NumPoints = (numX+2) * (numY+2);
	std::cout << "Writing in the Coefficients file" << std::endl;
	ofstream fpout("sim.coeff", ios_base::out | ios_base::app);
	
	printarray(fpout, NumPoints, ae);
	printarray(fpout, NumPoints, aw);
	printarray(fpout, NumPoints, an);
	printarray(fpout, NumPoints, as);
	printarray(fpout, NumPoints, ap);
	printarray(fpout, NumPoints, b);
	
	return 0;
}

void printarray(ofstream& out, int numPoints, double *c)
{
	int i;
	for (i = 0; i < numPoints; i++) {
		out << c[i] << " ";
	}
	out << endl;
}

int Csim2d::OutTmpp(int tstep)
{
    int i = 0,                  // counter in three dimensions
        j = 0, Nx = 0, sp = 0,  // cell index
        NumPoints = 0,          // number of points to write
        RecLen    = 0;          // record length of data
    ofstream *fpout;                // file for distribution function
    
    float K = 0, time, fs, temp, konz, p, vx, vy, R = 0., // radius of curvature
          R_average = 0.,  // average radius
          H_sum = 0.;      // sum of distribution

    float *Hau;
    time = tstep * (TimeInfo->tWidth);
    Nx = DimInfo->Nx;

    NumPoints = numX * numY;

    // allocate memory for distribution function
    Hau = new float[ DimInfo->Nx ];
    if (Hau == NULL) {
        cerr << "Error in allocating memory for Hau" << endl;
        exit(-1);
    }

    RecLen = sizeof(float) + sizeof(int) + NumPoints * sizeof(float);

    ofstream fpT(TMPPNAME, ofstream::binary | ofstream::app);
    ofstream fpF(FSPNAME, ofstream::binary | ofstream::app);
    ofstream fpC(KONZNAME, ofstream::binary | ofstream::app);
    ofstream fpP(PNAME, ofstream::binary | ofstream::app);
    ofstream fpVX(VXNAME, ofstream::binary | ofstream::app);
    ofstream fpVY(VYNAME, ofstream::binary | ofstream::app);

    if ( !(fpT.is_open() && fpF.is_open() && fpC.is_open() && fpP.is_open() 
            && fpVX.is_open() && fpVY.is_open()) ) {
        cerr << "Error in opening files" << endl;
        return -1;
    }

    fpT.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpT.write(reinterpret_cast<const char*>(&time), sizeof(float));
    fpT.write(reinterpret_cast<const char*>(&NumPoints), sizeof(int));
    
    fpF.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpF.write(reinterpret_cast<const char*>(&time), sizeof(float));
    fpF.write(reinterpret_cast<const char*>(&NumPoints), sizeof(int));
    
    fpC.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpC.write(reinterpret_cast<const char*>(&time), sizeof(float));
    fpC.write(reinterpret_cast<const char*>(&NumPoints), sizeof(int));

    fpP.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpP.write(reinterpret_cast<const char*>(&time), sizeof(float));
    fpP.write(reinterpret_cast<const char*>(&NumPoints), sizeof(int));

    fpVX.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpVX.write(reinterpret_cast<const char*>(&time), sizeof(float));
    fpVX.write(reinterpret_cast<const char*>(&NumPoints), sizeof(int));

    fpVY.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpVY.write(reinterpret_cast<const char*>(&time), sizeof(float));
    fpVY.write(reinterpret_cast<const char*>(&NumPoints), sizeof(int));

    sp = DimInfo->Nx + 1;

    for (i = 0; i < numY; i++) {
        for (j = 0; j < numX; j++) {
            fs = CellInfo[sp].FracSol;
            konz = CellInfo[sp].CG;
            p = CellInfo[sp].P;
            // XXX output theta instead of temperature
            temp = CellInfo[sp].Theta;
            vx = 5000. * (CellInfo[sp].VX + CellInfo[sp + 1].VX);
            vy = 5000. * (CellInfo[sp].VY + CellInfo[sp + Nx].VY);
            fpF.write(reinterpret_cast<const char*>(&fs),sizeof(float));
            fpC.write(reinterpret_cast<const char*>(&konz),sizeof(float));
            fpP.write(reinterpret_cast<const char*>(&p), sizeof(float));
            fpT.write(reinterpret_cast<const char*>(&temp),sizeof(float));
            sp++;
        }
        sp += 2;
    }
    fpT.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpF.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpC.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpP.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpVX.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpVY.write(reinterpret_cast<const char*>(&RecLen), sizeof(int));
    fpT.close();
    fpF.close();
    fpC.close();
    fpP.close();
    fpVX.close();
    fpVY.close();
    return 0;
}

int Csim2d::OutSimgeo(void)
{
    int RecLen1 = 0;
    ofstream fp(SIMGEONAME, ios_base::out | ios_base::binary);
    if (!fp.is_open()) {
        cerr << "Error data can't be opened! " << endl;
        return -1;
    }
    RecLen1 = 3 * sizeof(int) + 3 * sizeof(float);
    
    // fwrite(&RecLen1, sizeof(int), 1, fp);
    fp.write(reinterpret_cast<const char*>(&numX), sizeof(int));
    fp.write(reinterpret_cast<const char*>(&numY), sizeof(int));
    fp.write(reinterpret_cast<const char*>(&numZ), sizeof(int));
    
    fp.write(reinterpret_cast<const char*>(&LengthX), sizeof(float));
    fp.write(reinterpret_cast<const char*>(&LengthY), sizeof(float));
    fp.write(reinterpret_cast<const char*>(&LengthZ), sizeof(float));

    // fwrite(&RecLen1, sizeof(int), 1, fp);
    fp.close();
    return 0;
}

void Csim2d::OutConv(void)
{
  cout << "#" << 10000*V_control << " " << V_error << " "
               << NV << " D " << D_error << " " << ND
               << " P " << P_error << " " << NP 
               << " C " << C_error << " " << NC << endl;
  cout << "#frame moved by " << n_move << " gridpoints" << endl;
}
    
    
    


             
