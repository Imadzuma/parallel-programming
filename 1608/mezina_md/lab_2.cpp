#include <iostream>
#include <ctime>
#include <algorithm>
#include "mpi.h"

using namespace std;

typedef unsigned char TPixel;

void EnterInfo(int& size_x, int& size_y, int& intensity, bool& writing) {
	cout << "I'm proccess 0! Enter length: ";
	cin >> size_x;
	cout << "I'm proccess 0! Enter width: ";
	cin >> size_y;
	cout << "I'm proccess 0! Enter intensity: ";
	cin >> intensity;
	cout << "I'm proccess 0! Enter 1, if information sould be displayed, else enter any value: ";
	int res_write;
	cin >> res_write;
	if (res_write == 1)
		writing = true;
	else
		writing = false;
}

TPixel* GenerateMatix(const int size_x, const int size_y, const bool writing) {
	TPixel* matrix = new TPixel[size_x * size_y];
	for (int y = 0; y < size_y; ++y) {
		for (int x = 0; x < size_x; ++x)
			matrix[y*size_x + x] = rand() % 256;
	}
	if (writing) {
		cout << "I'm proccess 0! I generate matrix:";
		for (int y = 0; y < size_y; ++y) {
			cout << "\n\t";
			for (int x = 0; x < size_x; ++x)
				cout << (int)matrix[y*size_x + x] << " ";
		}
		cout << endl;
	}
	return matrix;
}

void TransferEnterInfo(int& size_x, int& size_y, int& intensity, bool& writing) {
	MPI_Bcast(&size_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&size_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&intensity, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&writing, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
}

int* FindSizeArea(const int num_proc, const int size_x, const int size_y) {
	int* size_area = new int[num_proc];
	int whole = (size_x * size_y) / num_proc;
	int frac = (size_x * size_y) % num_proc;
	for (int i = 0; i < num_proc; ++i) 
		size_area[i] = whole + ((i < frac) ? 1 : 0);
	return size_area;
}

int* FindBeginArea(const int num_proc, const int size_x, const int size_y) {
	int* begin_area = new int[num_proc];
	int whole = (size_x * size_y) / num_proc;
	int frac = (size_x * size_y) % num_proc;
	for (int i = 0; i < num_proc; ++i) 
		begin_area[i] = i * whole + min(i, frac);
	return begin_area;
}

TPixel* TransferMatrix(int& work_size, const int num_proc, const int rank, const int size_x, const int size_y, const bool writing, const TPixel* matrix, const int* size_area, const int* begin_area) {
	work_size = ((size_x * size_y) / num_proc) + ((rank < (size_x * size_y) % num_proc) ? 1 : 0);
	TPixel* work_area = new TPixel[work_size];
	MPI_Scatterv(matrix, size_area, begin_area, MPI_UNSIGNED_CHAR, work_area, work_size, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	if (writing) {
		cout << "I'm proccess " << rank << "! ";
		if (work_size == 0)
			cout << "My work area is empty" << endl;
		else {
			cout << "My work area is: ";
			for (int i = 0; i < work_size; ++i)
				cout << (int)work_area[i] << " ";
			cout << endl;
		}
	}
	return work_area;
}

double FindAverage(const int rank, const int work_size, const TPixel* work_area, const bool writing, const int size_x, const int size_y) {
	long long sum = 0;
	for (int i = 0; i < work_size; ++i)
		sum += work_area[i];
	long long res_sum;
	MPI_Allreduce(&sum, &res_sum, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	double average = (double)res_sum / (size_x * size_y);
	if (rank == 0 && writing)
		cout << "I'm proccess " << rank << "! Average is " << average << endl;
	return average;
}

TPixel* FindResult(const int rank, const int work_size, const TPixel* work_area, const double average, const int intensity, const bool writing) {
	TPixel* res_area = new TPixel[work_size];
	int table[256];
	double correct = 1.0 + intensity / 100.0;
	for (int i = 0; i <= 255; ++i) {
		double del = (double)i - average;
		table[i] = average + del * correct;
		table[i] = (table[i] < 0) ? 0 : table[i];
		table[i] = (table[i] > 255) ? 255 : table[i];
	}
	for (int i = 0; i < work_size; ++i) 
		res_area[i] = table[work_area[i]];
	if (writing) {
		cout << "I'm proccess " << rank << "! ";
		if (work_size == 0)
			cout << "I don't have work" << endl;
		else {
			cout << "My result is ";
			for (int i = 0; i < work_size; ++i)
				cout << (int)res_area[i] << " ";
			cout << endl;
		}
	}
	return res_area;
}

TPixel* CollectResult(const int rank, const int work_size, const TPixel* res_area, const int* size_area, const int* begin_area, const int size_x, const int size_y, const bool writing) {
	TPixel* res_matrix = new TPixel[size_x * size_y];
	MPI_Gatherv(res_area, work_size, MPI_UNSIGNED_CHAR, res_matrix, size_area, begin_area, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	if (rank == 0 && writing) {
		cout << "I'm proccess " << rank << "! We find result matrix:";
		for (int y = 0; y < size_y; ++y) {
			cout << "\n\t";
			for (int x = 0; x < size_x; ++x)
				cout << (int)res_matrix[y*size_x + x] << " ";
		}
		cout << endl;
	}
	return res_matrix;
}

void Check(const TPixel* matrix, const TPixel* res_matrix, const int size_x, const int size_y, const int intensity, const bool writing) {
	long long sum = 0;
	double start = MPI_Wtime();
	for (int y = 0; y < size_y; ++y) {
		for (int x = 0; x < size_x; ++x) 
			sum += matrix[size_x * y + x];
	}
	double average = (double)sum / (size_x * size_y);
	TPixel* check_matrix = FindResult(0, size_x * size_y, matrix, average, intensity, false);
	if (writing) {
		cout << "Checking matrix: ";
		for (int y = 0; y < size_y; ++y) {
			cout << "\n\t";
			for (int x = 0; x < size_x; ++x)
				cout << (int)check_matrix[y*size_x + x] << " ";
		}
		cout << endl;
	}
	bool check = true;
	for (int y = 0; y < size_y; ++y) {
		for (int x = 0; x < size_x; ++x) {
			if (check_matrix[size_x * y + x] != res_matrix[size_x * y + x]) {
				check = false;
				cout << "Error in x = " << x << ", y = " << y << endl;
			}
		}
	}
	double finish = MPI_Wtime();
	cout << "Time: " << finish - start << endl;
	if (check)
		cout << "Check finished with success!" << endl;
	else
		cout << "Check finished with failure!" << endl;
	delete[] check_matrix;
}

int main(int argc, char** argv) {
	cout << fixed;
	cout.precision(10);

	//time
	double start_time, end_time;

	srand(time(0));

	//initialization MPI
	int num_proc, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//enter sizes and intensity
	int size_x = -1, size_y = -1;
	int intensity;
	bool writing;
	if (rank == 0)
		EnterInfo(size_x, size_y, intensity, writing);

	//generate matrix
	TPixel* matrix = NULL;
	if (rank == 0) 
		matrix = GenerateMatix(size_x, size_y, writing);

	//transfer sizes and intensity
	if (rank==0)
		start_time = MPI_Wtime();
	TransferEnterInfo(size_x, size_y, intensity, writing);

	//find separation boundary
	int* size_area = NULL;
	int* begin_area = NULL;
	if (rank == 0) {
		size_area = FindSizeArea(num_proc, size_x, size_y);
		begin_area = FindBeginArea(num_proc, size_x, size_y);
	}

	//transfer matrix
	int work_size;
	TPixel* work_area = TransferMatrix(work_size, num_proc, rank, size_x, size_y, writing, matrix, size_area, begin_area);

	//find average
	double average = FindAverage(rank, work_size, work_area, writing, size_x, size_y);

	//find result in proccess
	TPixel* res_area = FindResult(rank, work_size, work_area, average, intensity, writing);

	//collect results
	TPixel* res_matrix = CollectResult(rank, work_size, res_area, size_area, begin_area, size_x, size_y, writing);

	//print time
	if (rank == 0) {
		end_time = MPI_Wtime();
		cout << "I'm proccess " << rank << "! Time: " << end_time - start_time << endl;
	}

	if (rank == 0)
		Check(matrix, res_matrix, size_x, size_y, intensity, writing);

	//free memory
	if (size_area != NULL) delete[] size_area;
	if (begin_area != NULL) delete[] begin_area;
	if (work_area != NULL) delete[] work_area;
	if (res_area != NULL) delete[] res_area;
	if (matrix != NULL) delete[] matrix;
	if (res_matrix != NULL) delete[] res_matrix;

	MPI_Finalize();
}