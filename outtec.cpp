// writing result to binary file to plot with tecplot.
// 
#include "../../sim2d.h"
#include "TECIO.h"


enum FileType {FULL=0, GRID=1, SOLUTION=2};


//int Csim2d::OutAsciiTec(int tstep)
//{
//
//}

/**
 * This function writes the solution to tecplot
 * \param tstep time step of the solution to write.
 *
 */
int 
Csim2d::OutTec(int tstep)
{
	char num[25];
	char name[255];
	double SolTime;
	int sp;
	float fs, temp, konz, p, vx, vy;
	INTEGER4 Debug,I,J,III,DIsDouble,VIsDouble,IMax,JMax,KMax,ZoneType,StrandID,ParentZn,IsBlock;
	INTEGER4 ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn, FileType;

	Debug     = 1;
	VIsDouble = 0;
	DIsDouble = 0;
	IMax      = numX;
	JMax      = numY;
	KMax      = 1;
	ZoneType  = 0;      /* Ordered */
	SolTime   = 360.0;
	StrandID  = 0;     /* StaticZone */
	ParentZn  = 0;      /* No Parent */
	IsBlock   = 1;      /* Block */
	ICellMax  = 0;
	JCellMax  = 0;
	KCellMax  = 0;
	NFConns   = 0;
	FNMode    = 0;
	ShrConn   = 0;
	FileType  = FULL;

	/*
	* Open the file and write the tecplot datafile 
	* header information 
	*/
	sprintf(name, "simt%i.plt", tstep);

	I = TECINI111("SIMPLE DATASET",
		"X Y FS C T P",
		name,
		".",
		&FileType,
		&Debug,
		&VIsDouble);

	// allocate data 
	
	float **X = (float**)malloc(sizeof(float*)*(numY+1));
	X[0] = (float*) malloc (sizeof(float)*numX*numY);
	float **Y = (float**)malloc(sizeof(float*)*(numY+1));
	Y[0] = (float*) malloc (sizeof(float)*numX*numY);
	float **FS = (float**)malloc(sizeof(float*)*(numY+1));
	FS[0] = (float*) malloc (sizeof(float)*numX*numY);
	float **CONC = (float**)malloc(sizeof(float*)*(numY+1));
	CONC[0] = (float*) malloc (sizeof(float)*numX*numY);
	float **T = (float**)malloc(sizeof(float*)*(numY+1));
	T[0] = (float*) malloc (sizeof(float)*numX*numY);
	float **P = (float**)malloc(sizeof(float*)*(numY+1));
	P[0] = (float*) malloc (sizeof(float)*numX*numY);

	for (int i=1; i<=numY; i++) {
		X[i] = X[i-1]+numX;
		Y[i] = Y[i-1]+numX;
		FS[i] = FS[i-1]+numX;
		CONC[i] = CONC[i-1]+numX;
		T[i] = T[i-1]+numX;
		P[i] = P[i-1]+numX;
	}

		
	
	sp = DimInfo->Nx + 1;
	int Nx = DimInfo->Nx;

	// write data to file
	for (I = 0; I < numY; I++) {
		for (J = 0; J < numX; J++) {
			fs = CellInfo[sp].FracSol;
			konz = CellInfo[sp].CG;
			p = CellInfo[sp].P;
			// XXX output theta instead of temperature
			temp = CellInfo[sp].Theta;
			vx = 5000. * (CellInfo[sp].VX + CellInfo[sp + 1].VX);
			vy = 5000. * (CellInfo[sp].VY + CellInfo[sp + Nx].VY);
			X[J][I] = J * LengthX;
			Y[J][I] = I * LengthY;
			FS[J][I] = fs;
			CONC[J][I] = konz;
			P[J][I] = p;
			T[J][I] = temp;
			
			sp++;
		}
		sp += 2;
	}
	
	/*
	* Write the zone header information.
	*/
	I = TECZNE111("Simple Zone",
		&ZoneType,
		&IMax,
		&JMax,
		&KMax,
		&ICellMax,
		&JCellMax,
		&KCellMax,
		&SolTime,
		&StrandID,
		&ParentZn,
		&IsBlock,
		&NFConns,
		&FNMode,
		0,              /* TotalNumFaceNodes */
		0,              /* NumConnectedBoundaryFaces */
		0,              /* TotalNumBoundaryConnections */
		NULL,           /* PassiveVarList */
		NULL,           /* ValueLocation = Nodal */
		NULL,           /* SharVarFromZone */
		&ShrConn);
	/*
	* Write out the field data.
	*/
	III = IMax*JMax;
	I   = TECDAT111(&III,&X[0][0],&DIsDouble);
	I   = TECDAT111(&III,&Y[0][0],&DIsDouble);
	I   = TECDAT111(&III,&FS[0][0],&DIsDouble);

	I   = TECDAT111(&III,&CONC[0][0],&DIsDouble);
	I   = TECDAT111(&III,&T[0][0],&DIsDouble);
	I   = TECDAT111(&III,&P[0][0],&DIsDouble);

	I = TECEND111();

	// deallocate mem
	free(FS[0]); 
	free (FS);
	free(CONC[0]); 
	free (CONC);
	free(T[0]); 
	free(T);
	free(P[0]); 
	free(P);
	free(X[0]); 
	free(X);
	free(Y[0]); 
	free(Y);
	return 0;
}
