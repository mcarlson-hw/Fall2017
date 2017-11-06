#include <iostream>
#include <mpi.h>
using namespace std;

#include "tests.h"
#include "SpectralClusterer.h";

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	MPI_Comm world = MPI_COMM_WORLD;

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	test1(world);
		
	MPI_Finalize();
	return 0;
}