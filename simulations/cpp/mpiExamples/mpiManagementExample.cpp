#include <fstream>
#include <iostream>
#include <mpi.h>

// mpic++ -o  mpiManagementExample mpiManagementExample.cpp 
// mpirun -np 4 --host localhost mpiManagementExample

int main(int argc,char **argv)
{
	int mpiResult;
	int numProcesses;
	int numtasks;
	int rank;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	int len;

	mpiResult = MPI_Init (&argc,&argv);
	if(mpiResult!= MPI_SUCCESS)
		{
			std::cout << "MPI not started. Terminating the process." << std::endl;
			MPI_Abort(MPI_COMM_WORLD,mpiResult);
		}

	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Get_processor_name(hostname, &len);
	std::cout << "Number of tasks= " <<  numtasks
						<< " My rank= " << rank
						<< " Running on " << hostname
						<< std::endl;

	MPI_Finalize();
}
