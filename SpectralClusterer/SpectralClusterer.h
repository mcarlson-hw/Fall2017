#pragma once
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>
#include <fstream>
#include <tuple>
#include <algorithm>
class SClusterer
{
	class Vertex;
	class Partition;
	class Subpartition;

	std::string filename, ext;
	std::vector<std::tuple<int, int, double>> V;
	int K, L, N, nonzeros;
	int rank, size;
	MPI_Comm comm;
	
	std::vector<Vertex> W;
	std::vector<int> lid;

	Partition *rp;
	Partition *cp;

	void parse();
	void prepare_local_vertices();
	void init_partition();
	void compute_spectrum();
	void applyW(Subpartition, std::vector<double>&, std::vector<double>&);

public:
	SClusterer(int, std::string, MPI_Comm);
	std::string getExt();
	void printW();
};

class SClusterer::Vertex
{
	int id, deg, vol;
	std::vector<int> e;
	std::vector<double> w, u;

public:
	Vertex(int, std::vector<int>, std::vector<double>);
	int get_deg();
};

class SClusterer::Partition
{
	int K;
	double lQ, Q;
	std::vector<Subpartition> P;

	MPI_Comm comm;

public:
	Partition(int, std::vector<std::vector<int>>, MPI_Comm);
};

class SClusterer::Subpartition
{
	int L;
	double lvol, vol, sQ, lsQ;
	std::vector<int> v;

	MPI_Comm comm;
	int rank, size;

public:
	Subpartition(std::vector<int>, MPI_Comm);
};