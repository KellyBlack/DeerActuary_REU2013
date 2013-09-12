/* ***************************************************************

   This code is subject to the MIT License (MIT)

   Copyright (c) 2013 Kelly Black, Candace Liu, Elizabeth Sweeney, Lindong Zhou

   Permission is hereby granted, free of charge, to any person
   obtaining a copy of this software and associated documentation
   files (the "Software"), to deal in the Software without
   restriction, including without limitation the rights to use, copy,
   modify, merge, publish, distribute, sublicense, and/or sell copies
   of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

   *****************************************************************

   Code to run a Monte Carlo simulation of the systems of equations
   representing the population of deer in one area and the associated
   funds in an account to offset insurance liabilities.

   ***************************************************************** */

#include <fstream>
#include <iostream>
#include <thread>
#include <mutex>
#include <cstdio>

#include <math.h>
#include <mpi.h>

/* Helper functions to set command line options. */
//#include <getopt.h>
#include <string.h>

#define DEFAULT_FILE "threaded_trial"
#define NUMBER_THREADS 4
#define DEBUG
#define VERBOSE

/* create a mutex that is used to protect the writing of the data to the file. */
std::mutex writeToFile;

/* Routine to calculate the step size for a given range and number of iterations. */
double calcDelta(double theMin,double theMax,int number)
{
  return((theMax-theMin)/((double)number));
}

/* Routine to write out the given data to the file. */
void printResultsCSV(double dt,int numberTimeSteps,
                     double P, double alpha,
                     double m1, double m2, 
                     double sumX, double sumX2,
                     double sumM, double sumM2,
                     int numberIters,
                     std::ofstream* dataFile)
{
  std::lock_guard<std::mutex> guard(writeToFile);  // Make sure that only this routine
                                                   // can access the file at any one time.
  *dataFile << dt*((double)numberTimeSteps) << "," 
            << P << "," 
            << alpha << "," 
            << m1 << "," << m2 << "," 
            << sumX << "," << sumX2 << "," 
            << sumM << "," << sumM2 << "," 
            << numberIters << std::endl;
  (*dataFile).flush();


}


/* Routine to write out the given data to the file. */
void printResultsMPI(double dt,int numberTimeSteps,
                     double P, double alpha,
                     double m1, double m2, 
                     double sumX, double sumX2,
                     double sumM, double sumM2,
                     int numberIters,
                     MPI_File* dataFile)
{
  std::lock_guard<std::mutex> guard(writeToFile);  // Make sure that only this routine
                                                   // can access the file at any one time.

  MPI_Status status;
  struct output 
  {
    double theTime;
    double P;
    double alpha;
    double m1;
    double m2;
    double sumX;
    double sumX2;
    double sumM;
    double sumM2;
    int iterations;
  };
  output basicStats;

  basicStats.theTime    = dt*((double)numberTimeSteps);
  basicStats.P          = P;
  basicStats.alpha      = alpha;
  basicStats.m1         = m1;
  basicStats.m2         = m2;
  basicStats.sumX       = sumX;
  basicStats.sumX2      = sumX2;
  basicStats.sumM       = sumM;
  basicStats.sumM2      = sumM2;
  basicStats.iterations = numberIters;
  MPI_File_write_ordered(*dataFile,&basicStats,sizeof(basicStats), MPI_INT, &status );
}




/* Routine to generate two normally distributed random numbers */
inline void randNormal(double nu[]) 
{
  double tmp = sqrt(-2.0*log(drand48()));
  double trig = 2.0*M_PI*drand48();
  nu[0] = tmp*sin(trig);
  nu[1] = tmp*cos(trig);
}



/* Routine to approximate one sample path. */
void samplePath(double P,double alpha,double beta,
                double r1,double h,double F,
                double rho,double g,
                int numberIters,
                double dt,
                double sdt,
                int numberTimeSteps,
                std::ofstream* dataFile
                //MPI_File* dataFile
                )
{

  /* Run time parameters */
  int timeLupe;
  double t;

  double m[2];
  double dW[2];

  double stochasticIntegral;  // The integral used for the sol. to the pop eqn.
  double z;                   // The transformed (linear) sol. to the pop eqn.
  double W;                   // The random walk.


  /* Determine and set the scaled parameters */
  double rtilde = r1-h;                     // scaled growth rate
  double ftilde = (rtilde/r1)*F;            // scaled carrying capacity
  double a      = rtilde-0.5*(alpha*alpha); // exp  exponent for sol. to deep eqn.
  double g0     = 0.5*alpha*alpha/a;        // int. constant for deer pop. solution.
  // todo - keep track of 1/a. 
  // keep track of exp(*) or factor it out appropriately?


  /* Start the loop for the multiple simulations. */
  double sumX  = 0.0;
  double sumX2 = 0.0;
  double sumM  = 0.0;
  double sumM2 = 0.0;
  int lupe;


  randNormal(dW); // calc. the initial set of random numbers.
  for(lupe=0;lupe<numberIters;++lupe)
    {
      /* set the initial conditions. */
      W    = 0.0;
      m[0] = ftilde;
      m[1] = (P-beta*ftilde)/(g-rho);
      stochasticIntegral = 0.0;

      for(timeLupe=0;timeLupe<numberTimeSteps;++timeLupe)
        {
          /* Set the time step. */
          t = ((double)timeLupe)*dt;

          /* Calc. two normally distributed random numbers */
          if(timeLupe%2==0)
            randNormal(dW); // calc. a new set of random numbers.
          else
            dW[0] = dW[1];    // shift the 2nd number into the first slot.
          dW[0] *= sdt;       // scale the change in W to have the proper variance.

          // Update the integral and then update the population and fund balance.
          stochasticIntegral += 
            exp(a*t+alpha*W)*(dW[0] + 0.5*alpha*(dW[0]*dW[0]-dt));
          z = rtilde/a - 
            exp(-a*t-alpha*W)*(g0 + ((alpha*rtilde)/a)*stochasticIntegral); 
          m[0] = ftilde/z;
          m[1] += (rho*m[1]+P-beta*m[0])*dt - beta*m[0]*dW[0] 
            - 0.5*alpha*beta*m[0]*(dW[0]*dW[0]-dt);

          W += dW[0];
        }

      // Update the tally used for the statistical ensemble
      sumX  += m[0];
      sumX2 += m[0]*m[0];
      sumM  += m[1]*1.0E-1;
      sumM2 += m[1]*m[1]*1.0E-2;
    }

  printResultsCSV(dt,numberTimeSteps,P,alpha,m[0],m[1],
                  sumX,sumX2,sumM,sumM2,numberIters,dataFile);

}




int main(int argc,char **argv)
{

  /* Define the basic run time variables. */
  double initialTime  =  0.0;
  double finalTime    = 10.0;
#ifdef DEBUG
  int numberIters     = 100;
  int numberTimeSteps = 5000;
#else
  int numberIters     = 100000;
  int numberTimeSteps = 500000;
#endif
  double dt;
  double sdt;


  /* define some book keeping variables. */
  /* Define the estimated parameters for the problem. */
  double r1   = log(1.702); // Deer max reproduction rate
  double h    = log(1.16);  // Harvest rate of the deer
  double F    = 28000.0;    // Carrying capacity of the deer.
  double rho  = 0.04;       // Bond fund rate of growth: log(1+rate); 
  double beta = 9.0;        // cost due to deer collisions .003*3000 */
  double g    = 0.05;       // Net target rate of growth of the fund.

  /* thread management */
  std::thread simulation[NUMBER_THREADS];
  int numberThreads = 0;


  /* define the parameters ranges*/
  double Pmin     = 430000.0;
  double alphaMin = 0.0;

  double Pmax     = 530000.0;
  double alphaMax = 0.15;

  double deltaP;
  double deltaAlpha;

#ifdef DEBUG
  int numP     = 100;
  int numAlpha = 100;
#else
  int numP     = 1000;
  int numAlpha = 1000;
#endif

   int boundsP[2];
   int lupeP,lupeAlpha;

  /* define the parameters */
  double P;
  double alpha;


   /* Define the output parameters. */
   char outFile[1024];

  /* File stuff */
  std::ofstream dataFile;
  //MPI_File mpiFileHandle;

   /* MPI related variables. */
   int  mpiResult;
   int  numMPITasks;
   int  mpiRank;
   char mpiHostname[MPI_MAX_PROCESSOR_NAME];
   int  mpiHostnameLen;
   MPI_Status mpiStatus;


   // initialize the mpi parameters
   mpiResult = MPI_Init (&argc,&argv);
   if(mpiResult!= MPI_SUCCESS)
     {
       std::cout << "MPI not started. Terminating the process." << std::endl;
       MPI_Abort(MPI_COMM_WORLD,mpiResult);
     }

   MPI_Comm_size(MPI_COMM_WORLD,&numMPITasks);
   MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
   MPI_Get_processor_name(mpiHostname, &mpiHostnameLen);
   MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
#ifdef VERBOSE
   std::cout << "Number of tasks= " <<  numMPITasks
                                         << " My rank= " << mpiRank
                                         << " Running on " << mpiHostname
                                         << std::endl;
#endif

  /* Set the step values for the parameters. */
  deltaP     = calcDelta(Pmin,Pmax,numP);
  deltaAlpha = calcDelta(alphaMin,alphaMax,numAlpha);

  /* Sort out which process is doing what. */
  if(mpiRank == 0)
    {
      // Figure out who does what.
      int totalP = numP/numMPITasks;
      for(lupeP=0;lupeP<(numMPITasks-1);++lupeP)
        {
          boundsP[0] = lupeP*totalP;
          boundsP[1] = (lupeP+1)*totalP;
          MPI_Send(boundsP,2,MPI_INT,lupeP+1,10,MPI_COMM_WORLD);
        }

      // Assign my own bounds to what is left.
      boundsP[0] = (numMPITasks-1)*totalP;
      boundsP[1] = numP;
    }

  else
    {
      // Get my bounds for the values of P to use.
      MPI_Recv(boundsP,2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&mpiStatus);
      boundsP[1] -= 1;
    }
#ifdef VERBOSE
  std::cout << "Process: " << mpiRank << " Got the bound: " << boundsP[0] << "," << boundsP[1] << std::endl;
  //MPI_Finalize();
  //return(0);
#endif


  /* Set the number of iterations used in the main loop. */
  dt  = ((finalTime-initialTime)/((double)numberTimeSteps));
  sdt = sqrt(dt);

#ifdef VERBOSE
  std::cout << "Starting iteration. " << numberTimeSteps << " iterations." << std::endl;
#endif

  /* Open the output file and print out the header. */
  sprintf(outFile,"%s-%d.dat",DEFAULT_FILE,mpiRank);
  dataFile.open(outFile,std::ios::out);
  dataFile << "time,P,alpha,x,m,sumx,sumx2,summ,summ2,N" << std::endl;
  std::cout << "opening " << outFile << std::endl;

  /*
  MPI_Status status;
  char err_buffer[MPI_MAX_ERROR_STRING];
  int resultlen;
  int ierr = MPI_File_open(MPI_COMM_WORLD,outFile, 
                           MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL, 
                           MPI_INFO_NULL, 
                           &mpiFileHandle);
  std::cout << "Open: " << ierr << "," << MPI_SUCCESS << std::endl;
  std::cout << status.MPI_ERROR << "," << status.MPI_SOURCE << "," << status.MPI_TAG << std::endl;
  MPI_Error_string(ierr,err_buffer,&resultlen);
  std::cout << "Error: " << err_buffer << std::endl;
  */

  //printResultsMPI(0.1,100,2.0,3.0,4.0,5.0,
  /*                6.0,7.0,8.0,9.0,10,&mpiFileHandle);
  lupeP = 1;
  //ierr = MPI_File_write_ordered(mpiFileHandle,&lupeP,sizeof(lupeP), MPI_INT, &status );
  std::cout << "Write: " << ierr << "," << MPI_SUCCESS << std::endl;
  std::cout << status.MPI_ERROR << "," << status.MPI_SOURCE << "," << status.MPI_TAG << std::endl;
  MPI_Error_string(ierr,err_buffer,&resultlen);
  std::cout << "Error: " << err_buffer << std::endl;
  MPI_File_close(&mpiFileHandle);
  MPI_Finalize();
  exit(0);
  */

  /* Set the seed for the random number generator. */
  srand48(time(NULL));


  /* Get the details of the mpi processes. */


  /* Go through and run the simulations for all possible values of the parameters. */
  for(lupeP=boundsP[0];lupeP<=boundsP[1];++lupeP)
    {
      P = Pmin + deltaP*((double)lupeP);

      for(lupeAlpha=0;lupeAlpha<=numAlpha;++lupeAlpha) 
        {
          alpha = alphaMin + deltaAlpha*((double)lupeAlpha);

          if(numberThreads >= NUMBER_THREADS)
            {
              // There are too many threads. Wait for each run to end.
              while(numberThreads>0)
                {
#ifdef DEBUG
                  std::cout << "Waiting on thread " << simulation[numberThreads-1].get_id() 
                            << std::endl;
#endif
                  simulation[--numberThreads].join();
                }
            }

          // Make this run a separate thread.
          simulation[numberThreads++] = std::thread(samplePath,
                                                    P,alpha,beta,
                                                    r1,h,F,rho,g,
                                                    numberIters,dt,sdt,
                                                    numberTimeSteps,
                                                    &dataFile); //mpiFileHandle);
#ifdef VERBOSE
          /* print a notice */
          std::cout << "Simulation: " 
                    << dt*((double)numberTimeSteps) << "," 
                    << P << "," << alpha << "," 
                    << simulation[numberThreads-1].get_id() << std::endl;
#endif
        }



    }

  // Wait until all threads are done.
  while(numberThreads>0)
    {
#ifdef DEBUG
      std::cout << "Waiting on thread " << simulation[numberThreads-1].get_id() 
                << std::endl;
#endif
      simulation[--numberThreads].join();
    }


  dataFile.close();
  //MPI_File_close(&mpiFileHandle);
  MPI_Finalize();
  return(0);
}

