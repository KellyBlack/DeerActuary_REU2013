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
			val[lupe] = (double)lupe;
		}

	if(rank == 0)
		{
			// This is the first process.
			std::cout << "Sending message from the zero process" << std::endl;
			MPI_Send(val,NUMBER,MPI_DOUBLE,1,10,MPI_COMM_WORLD);
			std::cout << "Zero process waiting to hear the message from " << numtasks-1 << std::endl;
			MPI_Recv(val,NUMBER,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&theStatus);
			std::cout << "Zero process heard from  " << theStatus.MPI_SOURCE << " with tag " 
								<< theStatus.MPI_TAG << std::endl;
		}
	else
		{
			// Wait to hear the message.
			std::cout << "Process " << rank << " waiting to hear the message" << std::endl;
			MPI_Recv(val,NUMBER,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&theStatus);
			std::cout << "Process " << rank << " heard from  " << theStatus.MPI_SOURCE << " with tag " 
								<< theStatus.MPI_TAG << std::endl;

			// update what was sent.
			for(lupe=0;lupe<NUMBER;++lupe)
				val[lupe] += 1.0;

			// pass it along to the next process
			std::cout << "Process " << rank << " sending message to " << (rank+1)%numtasks << " process." << std::endl;
			MPI_Send(val,NUMBER,MPI_DOUBLE,(rank+1)%numtasks,10,MPI_COMM_WORLD);
		}

	// print out what was passed
	std::cout << "Process " << rank << " heard: ";
	for(lupe=0;lupe<NUMBER;++lupe)
		std::cout << val[lupe] << ", ";
	std::cout << std::endl;

	MPI_Finalize();
}
