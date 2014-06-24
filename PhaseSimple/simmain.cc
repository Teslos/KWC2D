#include <cstdlib>
#include "sim2d.h"
using namespace std;

Csim2d *csim2d;

int 
main(int argc, char **argv)
{
   // local variables
   int RetVal;                  // return value
   int i,j;
   float DRetVal;               // float return value
   int WriteInterval;           // when to write result
   int WriteInterval1;          // control values in logfile 
   int WriteInterval2;          // convergence in logfile 
   char *inname;                // input file name
   FILE *fd;
   struct CTime *TimeInfo;   
   
   // remove .tmpp, .fsp, .konz and .geop files
   system("rm sim.tmpp");
   system("rm sim.fsp");
   system("rm sim.druck");
   system("rm sim.konz");
   system("rm sim.vx sim.vy");
   system("rm dis.*");

   csim2d = new Csim2d();
   // read input and initialize RegInf
   // for now all input comes from std stream cin
   // but it possible also to read the stream from file
   RetVal = csim2d->InitRegion(cin);
   if (RetVal < 0) {
       cerr << "Error in initialization of: Region Info" << endl;
       exit(-1);
   }

   RetVal = csim2d->OutSimgeo();
   if (RetVal < 0) {
       cerr << "Error in writing the simgeo-files!" << endl;
   }
   cout << "############################" << endl;

   // initialize Dimension information and working region
   RetVal = csim2d->InitDim();
   if (RetVal < 0) {
       cerr << "Error in initialization of dimensions: DimInfo" << endl;
       exit(-1);
   }
   
   // read in time data and initialize the Time Info
   RetVal = csim2d->InitTime(cin);
   if (RetVal < 0) {
       cerr << "Error in initialization of time info: TimeInfo" << endl;
       exit(-1);
   }

   // see if system has a chance to converge
   //float DRet = CheckConv();
   //if (DRet > 0.0) {
   //    cerr << " Quotient: " << DRet << endl;
   //    exit(-1);
   //}
   // get the data from RegInf and initialize CellInfo
   RetVal = csim2d->InitCell();
   if (RetVal < 0) {
       cerr << "Error in initialization of Cell informations: CellInfo" << endl;
       exit(-1);
   }

   cout << "# Initialization of CellInfo is successfull" << endl;
   // set temperature field in the CellInfo, which is
   // needed for Temperature solver
   RetVal = csim2d->InitTemperature();
   if (RetVal < 0) {
       cerr << "Error in initialization of Temperature" << endl;
       exit(-1);
   }
   // sets the initial boundary
   csim2d->InitBound();
   cout << "# Initial boundary set" << endl;
   // read in for which timesteps is data written out
   cout << "# After how many timesteps are result data written out?" << endl;
   cin >> WriteInterval;
   cout << "#" << WriteInterval << endl;
   // read in for which timesteps are control values written out
   cout << "# After how many timesteps are control values written out" << endl;
   cin >> WriteInterval1;
   cout << "#" << WriteInterval1 << endl;
   //RegInf->SI1 = (float)(WriteInterval1);
   cout << "# After how many timesteps are convergence values written out" << endl;
   cin >> WriteInterval2;
   cout << "#" << WriteInterval2 << endl;
   // now we can start for t_steps we will solve fields in complete region
   cout << "# Output: Timestep (filter) and Control values:" << endl;
   cout << "#t,Dvy/Dvx/Dv/Vmax/ytip/xtip/v/rho/s*/Pe/Utp/Otp/FS/iface/Cav"
       << "/Cav_l/T_bot/C_bot/ymax,vy-mean" << endl;
   TimeInfo = csim2d->getTimeInfo();
   csim2d->alloc_c();
   
   // now check stability condition
   bool stable = csim2d->StabilCond();
   if (!stable) {
       cout << "# *** Stability condition not fullfilled! ***" << endl;
       return -1;
   }


   csim2d->OutTmpp(0);
   for (i = 0; i < TimeInfo->tSteps; i++) {
       TimeInfo->time = TimeInfo->tWidth * (float)i;
       // for every write interval results are written out
       if ((i % WriteInterval) == 0) {
           cout << "# *** Timestep " << i << " ***" << endl;
           RetVal = csim2d->OutTmpp(i);
#ifdef __TECPLOT__
 	   RetVal = csim2d->OutTec(i);
#endif
           if (RetVal < 0) {
               cerr << "Error in writing output data" << endl;
               exit(-1);
           }
       }
       // writing out control values
       if ((i % WriteInterval2) == 0) {
           RetVal = csim2d->OutContr(i); // write DivAV, DivMax, and 3*VSuper
           if (RetVal < 0) {
               cerr << "Error in calculating control values" << endl;
               exit(-1);
           }
           // on begin set N0
           if (i == 0) {
               //csim2d->NR0 = csim2d->NR;
           }
       }
       // number of iterations/convergence 
       if (((i-1) % ((int) WriteInterval2)) == 0) {
           csim2d->OutConv();
       }
       csim2d->Extruder();
       // solve the temperature, phases and diffusion for the region
       RetVal = csim2d->solve2d();
   }
   csim2d->dealloc_c();

   return 0;
}
   

