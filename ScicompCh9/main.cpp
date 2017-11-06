#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <mkl.h>
using namespace std;

// File I/O
void read_csv(string filename, vector<double> *A, int* rows, int* cols)
{
	ifstream infile(filename);
	if (infile.is_open())
	{
		*rows = 0;
		*cols = 0;
		string line;
		string val;
		cout << "Starting read.\n";
		while (getline(infile, line, '\n'))
		{
			istringstream iline(line);
			*cols = 0;
			while (getline(iline, val, ','))
			{
				(*A).push_back(stod(val));
				*cols += 1;
			}
			*rows += 1;
		}
	}	
}
void write_binary(string filename, string ext, vector<double> *A, int rows, int cols)
{
	string temp = filename + "_size." + ext;
	string temp2 = filename + "." + ext;
	
	FILE *file1 = fopen(temp.c_str(), "wb");
	fwrite(&rows, sizeof(int), 1, file1);
	fwrite(&cols, sizeof(int), 1, file1);
	fclose(file1);

	FILE *file2 = fopen(temp2.c_str(), "wb");
	fwrite((*A).data(), rows*cols * sizeof(double), 1, file2);
	fclose(file2);
}
void read_binary(string filename, string ext, vector<double> *A, int* rows, int* cols)
{
	string temp = filename + "_size." + ext;
	string temp2 = filename + "." + ext;

	FILE *file1 = fopen(temp.c_str(), "rb");
	fread(rows, sizeof(int), 1, file1);
	fread(cols, sizeof(int), 1, file1);
	fclose(file1);

	FILE *file2 = fopen(temp2.c_str(), "rb");
	(*A).resize((*rows)*(*cols));
	fread((*A).data(), sizeof(double)*(*rows)*(*cols), 1, file2);
	fclose(file2);
}
void print_mat(vector<double> *A, int rows, int cols)
{
	if (rows > 10 && cols > 10)
	{
		for (int row = 0; row < 4; row++)
		{
			for (int col = 0; col < 4; col++)
			{
				cout << (*A)[row*cols + col] << " ";
			}
			cout << ". . . ";
			for (int col = cols - 5; col < cols; col++)
			{
				cout << (*A)[row*cols + col] << " ";
			}
			cout << "\n";
		}
		cout << "\n.\n.\n.\n\n";
		for (int row = rows - 5; row < rows; row++)
		{
			for (int col = 0; col < 4; col++)
			{
				cout << (*A)[row*cols + col] << " ";
			}
			cout << ". . . ";
			for (int col = cols - 5; col < cols; col++)
			{
				cout << (*A)[row*cols + col] << " ";
			}
			cout << "\n";
		}
	}
	else
	{
		for (int row = 0; row < rows; row++)
		{
			for (int col = 0; col < cols; col++)
			{
				cout << (*A)[row*cols + col] << " ";
			}
			cout << "\n";
		}
	}
}
void prepare_data(string filename, string binary_ext, string text_ext, vector<double> *A, int *rows, int *cols)
{
	ifstream infile(filename + "." + binary_ext);
	if (infile.good())
	{
		cout << "Reading binary file.\n";
		read_binary(filename, binary_ext, A, rows, cols);
	}
	else
	{
		cout << "Binary file does not exist! Attempting to read CSV and produce binary file.\n";
		read_csv(filename + "." + text_ext, A, rows, cols);
		write_binary(filename, binary_ext, A, *rows, *cols);
	}
}

// Data Preperation
void center_data(vector<double> *A, int rows, int cols, int S_type, vector<double> *means, vector<double> *s)
{
	// S_type:
	//	1: 1/n * one_norm
	//	2: 1/sqrt(n) * two_norm	(default)
	//	3: 1/sqrt(n-1) * two_norm
	//	4: infinity_norm
	double oldmean, mean, var, S, x;
	int imax;
	for (int row = 0; row < rows; row++)
	{
		mean = 0;
		for (int col = 0; col < cols; col++)
		{
			x = (*A)[row*cols + col];
			mean = mean + x;
		}
		mean = mean / cols;
		(*means).push_back(mean);
		for (int col = 0; col < cols; col++)
			(*A)[row*cols + col] = (*A)[row*cols + col] - mean;

		switch (S_type)
		{
		case 0:
			S = 1.0;
			break;
		case 1:
			S = 1.0 / cols * cblas_dasum(cols, &(*A)[row*cols], 1);
			break;
		case 2:
			S = 1.0 / sqrt(cols) * cblas_dnrm2(cols, &(*A)[row*cols], 1);
			break;
		case 3:
			S = 1.0 / sqrt(cols - 1) * cblas_dnrm2(cols, &(*A)[row*cols], 1);
			break;
		case 4:
			imax = cblas_idamax(cols, &(*A)[row*cols], 1);
			S = fabs((*A)[row*cols + imax]);
			break;
		default:
			S = 1.0 / sqrt(cols) * cblas_dnrm2(cols, &(*A)[row*cols], 1);
			break;
		}
		(*s).push_back(S);
		cblas_dscal(cols, 1.0/S, &(*A)[row*cols], 1);
	}
}


int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	vector<double> A;
	vector<double> lA;
	vector<double> s, u, vt, superb;
	int rows = 0;
	int cols = 0;
	int info = 0;
	vector<double> means, S;
	// Read File
	if (rank == 0)
	{
		double start = MPI_Wtime();
		prepare_data("open", "b", "csv", &A, &rows, &cols);

		// Compute SVD of A:
		s.resize(max(1, min(rows,cols)));
		u.resize(rows*rows);
		vt.resize(cols*cols);
		superb.resize(min(rows, cols) - 1);
		info = LAPACKE_dgesvd(CblasRowMajor, 'A', 'A', rows, cols, A.data(), cols, s.data(), u.data(), rows, vt.data(), cols, superb.data());
		double end = MPI_Wtime();

		cout << "\nInfo = " << info << "\n";
		cout << "\nU = \n";
		print_mat(&u, rows, rows);
		cout << "\nSingular Values: ";
		print_mat(&s, 1, min(rows, cols));
		cout << "\n";
		cout << "Number of singular values: " << min(rows, cols) << "\n";
		cout << "Time Elapsed: " << end - start << "\n\n";
	}

	// Split data among processors
	/*MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int* NPerP = new int[size];
	int* N0 = new int[size]();
	for (int i = 0; i < size; i++)
		NPerP[i] = (i < rows % size) ? cols*(rows / size + 1) : (cols*(rows / size));
	for (int i = 1; i < size; i++)
		N0[i] = N0[i - 1] + NPerP[i - 1];
	lA.resize(NPerP[rank]);
	MPI_Scatterv(A.data(), NPerP, N0, MPI_DOUBLE, lA.data(), NPerP[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);*/

	// Center Data

	

	MPI_Finalize();
	return 0;
}