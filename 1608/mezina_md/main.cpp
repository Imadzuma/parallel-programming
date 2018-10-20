#define _CRT_SECURE_NO_WARNINGS

#include <algorithm>
#include <cstdio>
#include <mpi.h>
#include <iostream>
#include <vector>
#include <ctime>

int main(int argc, char **argv)
{
	int rank, num_proc;
	int size_v;
	int* v = NULL;
	srand(time(0));

	//Èíèöèàëèçàöèÿ MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//Îïðåäåëåíèå è ïåðåäà÷à ðàçìåðîâ âåêòîðà
	if (rank == 0) {
		std::cout << "I'm process " << rank << "! Please, enter vector size: ";
		std::cin >> size_v;
		//std::cout << std::endl;
	}
	MPI_Bcast(&size_v, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Îïðåäåëåíèå è ïåðåäà÷à âåêòîðà
	v = new int[size_v];
	if (rank == 0) {
		std::cout << "I'm process " << rank << "! I generate vector:";
		for (int i = 0; i < size_v; ++i) {
			v[i] = rand() % 200 - 100;
			std::cout << " " << v[i];
		}
		std::cout << std::endl;
	}
	MPI_Bcast(v, size_v, MPI_INT, 0, MPI_COMM_WORLD);

	//Ïîäñ÷åò ðåçóëüòàòà â ïðîöåññàõ
	int size_area = (size_v - 1) / num_proc + ((rank < (size_v - 1) % num_proc) ? 1 : 0);
	int begin = ((size_v - 1) / num_proc) * rank + std::min((size_v - 1) % num_proc, rank);
	int res_proc = 0;
	for (int i = 0; i < size_area; ++i)
		res_proc += ((v[begin + i] >= 0) != (v[begin + i + 1] >= 0)) ? 1 : 0;
	std::cout << "I'm process " << rank << "! My area is " << begin << " to " << begin + size_area << ", my result is " << res_proc << std::endl;

	//Ïîäñ÷åò îáùåãî ðåçóëüòàòà
	int res;
	MPI_Reduce(&res_proc, &res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		std::cout << "I'm process " << rank << "! We find result! It is " << res << std::endl;
	}
	delete[] v;
	MPI_Finalize();
	return 0;
}
