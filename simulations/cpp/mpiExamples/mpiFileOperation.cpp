#include <fstream>
#include <iostream>
#include <mpi.h>

// mpic++ -o  mpiManagementExample mpiManagementExample.cpp 
// mpirun -np 4 --host localhost mpiManagementExample

#define NUMBER 10

int main(int argc,char **argv)
{
	int mpiResult;
	int numProcesses;
	int numtasks;
	int rank;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	int len;

	double val[NUMBER];
	int lupe;
	MPI_Status theStatus;

  // File related stuffs
  char outFile[1024];
  MPI_File mpiFileHandle;


	mpiResult = MPI_Init (&argc,&argv);
	if(mpiResult!= MPI_SUCCESS)
		{
			std::cout << "MPI not started. Terminating the process." << std::endl;
			MPI_Abort(MPI_COMM_WORLD,mpiResult);
		}

	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Get_processor_name(hostname, &len);
	//std::cout << "Number of tasks= " <<  numtasks
	//					<< " My rank= " << rank
	//					<< " Running on " << hostname
	//					<< std::endl;

	// Initialize the buffer
	for(lupe=0;lupe<NUMBER;++lupe)
		{
			val[lupe] = (double)(lupe+rank);
		}

  // open the file
  sprintf(outFile,"fileExample-%d.dat",rank);
  std::cout << "opening " << outFile << std::endl;
  MPI_Status status;
  char err_buffer[MPI_MAX_ERROR_STRING];
  int resultlen;
  int ierr = MPI_File_open(MPI_COMM_WORLD,outFile, 
                           MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL, 
                           MPI_INFO_NULL, 
                           &mpiFileHandle);
  std::cout << "Open: " << ierr << "," << MPI_SUCCESS << std::endl;
  std::cout << "Status: " << status.MPI_ERROR << "," << status.MPI_SOURCE << "," << status.MPI_TAG << std::endl;
  MPI_Error_string(ierr,err_buffer,&resultlen);
  std::cout << "Error: " << err_buffer << std::endl;


	// print out what was passed
	std::cout << "Process " << rank << " writing: ";
	for(lupe=0;lupe<NUMBER;++lupe)
		std::cout << val[lupe] << ", ";
	std::cout << std::endl;

  MPI_File_close(&mpiFileHandle);
	MPI_Finalize();
}
