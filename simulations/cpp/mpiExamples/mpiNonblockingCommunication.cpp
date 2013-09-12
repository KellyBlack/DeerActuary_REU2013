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
	MPI_Status  theStatus[2];
	MPI_Request theRequests[2];

	int test1;
	int test2;

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
			val[lupe] = (double)lupe + rank;
		}

	// pass the info along to the next process
	std::cout << "Process " << rank << " sending message to " << (rank+1)%numtasks << " process." << std::endl;
	MPI_Isend(val,NUMBER,MPI_DOUBLE,(rank+1)%numtasks,10,MPI_COMM_WORLD,&theRequests[0]);

	// Wait to hear the message.
	std::cout << "Process " << rank << " waiting to hear the message" << std::endl;
	MPI_Irecv(val,NUMBER,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&theRequests[1]);

	// keep checking to see when everything is done.
	// This is not a good way to do it, and it is only for demonstration purposes.
	test1 = 0;
	test2 = 0;
	while(!test1 && !test2)
		{
			MPI_Test (&theRequests[0],&test1,&theStatus[0]);
			MPI_Test (&theRequests[1],&test2,&theStatus[1]);
		}

	// better yet just wait for them all at once....
	MPI_Waitall(2, theRequests, theStatus);

	// print out what was passed
	std::cout << "Process " << rank << " heard: ";
	for(lupe=0;lupe<NUMBER;++lupe)
		std::cout << val[lupe] << ", ";
	std::cout << std::endl;

	MPI_Finalize();
}
