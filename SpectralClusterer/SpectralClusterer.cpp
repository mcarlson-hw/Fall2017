#include "SpectralClusterer.h"
SClusterer::SClusterer(int _K, std::string __filename, MPI_Comm _comm)
{
	comm = _comm;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
	filename = __filename;

	K = _K;

	parse();
	if (rank == 0) { std::cout << "Weight Matrix:\n"; printW(); std::cout << "\n"; }
	prepare_local_vertices();
	//init_partition();
}
SClusterer::Vertex::Vertex(int _id, std::vector<int> _e, std::vector<double> _w)
{
	id = _id;
	vol = 0;
	for (int d : _e)   e.push_back(d);
	for (double x : _w) { w.push_back(x); vol += x; }
	deg = e.size();
}
SClusterer::Partition::Partition(int _K, std::vector<std::vector<int>> Vk, MPI_Comm _comm)
{
	K = 1;
	lQ = 0.0;
	Q = 0.0;
	comm = _comm;
	P.push_back(Subpartition(Vk[0], _comm));
}
SClusterer::Subpartition::Subpartition(std::vector<int> _v, MPI_Comm _comm)
{
	lvol = 0.0;
	lsQ = 0.0;
	// On v, push back elements of _v that are in the local set of vertices
	// Compute local volume; Every time an element of _v is added, lvol += _v.weight
	// Gather local volumes at root, sum them, broadcast sum as global volume
	// Compute local sQ, gather at root, broadcast sum as global sQ
}

bool indexcompare(const std::tuple<int, int, double> &lhs, const std::tuple<int, int, double> &rhs)
{
	if (std::get<0>(lhs) != std::get<0>(rhs)) return std::get<0>(lhs) < std::get<0>(rhs);
	else return std::get<1>(lhs) < std::get<1>(rhs);
}
void SClusterer::parse()
{
	if (rank != 0) return;

	size_t ind = filename.find_last_of('.');
	ext = filename.substr(ind+1);

	std::ifstream inputfile(filename);
	std::string line;
	int a, b, c;
	double x;

	// Read first line:
	//	rows cols nonzeros
	inputfile >> a >> b >> c;
	N = a;
	N = b;
	nonzeros = c;

	// Read next 'nonezeros' liens
	//	i j value
	for (int line = 0; line < nonzeros; line++)
	{
		inputfile >> a >> b >> x;
		V.push_back(std::make_tuple(a, b, x));
		
		// matrix is symmetric, A(i,j) = A(j,i)
		if (a != b) V.push_back(std::make_tuple(b, a, x));
	}

	std::sort(V.begin(), V.end(), indexcompare);

	// Prepare CSR arrays ia, ja, a
	/*int globalcount = 0;
	int localcount = 0;
	for (int row = 1; row <= M; row++)
	{
		localcount = 0;
		while (globalcount + localcount < nonzeros && std::get<0>(W[globalcount + localcount]) == row)
		{
			jw.push_back(std::get<1>(W[globalcount + localcount]) - 1);
			w.push_back(std::get<2>(W[globalcount + localcount]));
			localcount += 1;
		}
		iw.push_back(globalcount);
		globalcount += localcount;
	}*/
}
void SClusterer::prepare_local_vertices()
{
	std::vector<int> I, J, lI, lJ;
	std::vector<double> w, lw;
	std::vector<int> VertsPerP, NZPerP;
	std::vector<int> Verts0, NZ0;
	int lnz;
	double curr_w;
	if (rank == 0)
	{
		for (int i = 0; i < N; i++)
		{
			I.push_back(std::get<0>(V[i]) - 1);
			J.push_back(std::get<1>(V[i]) - 1);
			w.push_back(std::get<2>(V[i]));
		}

		for (int i = 0; i < size; i++)
			VertsPerP.push_back((i < N % size) ? N / size + 1 : N / size);
		Verts0.push_back(0);
		for (int i = 1; i < size; i++)
			Verts0.push_back(Verts0[i - 1] - VertsPerP[i - 1]);
		

		std::cout << "Size = " << size << ", VertsPerP[0] = " << VertsPerP[0] << ", Verts0[0] = " << Verts0[0] << "\n";
		for (int i = 0; i < size; i++)
		{
			curr_w = 0.0;
			for (int j = Verts0[i]; j < Verts0[i] + VertsPerP[i]; j++)
				curr_w += W[j].get_deg();
			NZPerP.push_back(curr_w);
		}
		/*NZ0.push_back(0);
		for (int i = 1; i < size; i++)
			NZ0.push_back(NZ0[i - 1] - NZPerP[i - 1]);*/
	}

	//lI.resize(lnz);
	//lJ.resize(lnz);
	//lw.resize(lnz);
	//MPI_Scatter(NZPerP.data(), 1, MPI_INT, &lnz, 1, MPI_INT, 0, comm);
	/*MPI_Scatterv(I.data(), NZPerP.data(), NZ0.data(), MPI_INT, lI.data(), lnz, MPI_INT, 0, comm);
	MPI_Scatterv(J.data(), NZPerP.data(), NZ0.data(), MPI_INT, lJ.data(), lnz, MPI_INT, 0, comm);
	MPI_Scatterv(w.data(), NZPerP.data(), NZ0.data(), MPI_INT, lw.data(), lnz, MPI_INT, 0, comm);*/

	/*int count = 0;
	int vcount = -1;
	std::vector<int> e;
	std::vector<double> wp;
	for (int i = 0; i < lnz; i++)
	{
		if (I[i] != vcount)
		{
			if (count != 0) lid.push_back(count);
			count = 0;
			vcount = I[i];
			W.push_back(Vertex(vcount, e, wp));
			e = std::vector<int>();
			wp = std::vector<double>();
			e.push_back(lJ[i]);
			wp.push_back(lw[i]);
			count += 1;
		}
		else
		{
			e.push_back(lJ[i]);
			wp.push_back(lw[i]);
			count += 1;
		}
	}
	lid.push_back(count);
	W.push_back(Vertex(vcount, e, wp));*/
}
void SClusterer::init_partition()
{
	std::vector<std::vector<int>> v;
	std::vector<int> v0;
	for (int i = 0; i < N; i++)
		v0.push_back(i);
	v.push_back(v0);
	rp = &Partition(1, v, comm);
}
void SClusterer::compute_spectrum()
{

}
void SClusterer::applyW(Subpartition p, std::vector<double> &x, std::vector<double> &y)
{

}


// Misc. Utilities
std::string SClusterer::getExt() { return ext; }
void SClusterer::printW()
{
	if (nonzeros < 20)
		for (int i = 0; i < nonzeros; i++)
			std::cout << "Entry " << i << ": " << std::get<0>(V[i]) << ", " << std::get<1>(V[i]) << ", " << std::get<2>(V[i]) << std::endl;
	else
	{
		for (int i = 0; i < 10; i++)
			std::cout << "Entry " << i << ": " << std::get<0>(V[i]) << ", " << std::get<1>(V[i]) << ", " << std::get<2>(V[i]) << std::endl;
		std::cout << ".\n";
		std::cout << ".\n";
		std::cout << ".\n";
		for (int i = nonzeros-10; i < nonzeros; i++)
			std::cout << "Entry " << i << ": " << std::get<0>(V[i]) << ", " << std::get<1>(V[i]) << ", " << std::get<2>(V[i]) << std::endl;
	}
}
int SClusterer::Vertex::get_deg() { return deg; }

