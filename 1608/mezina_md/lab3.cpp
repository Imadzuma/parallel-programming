#include <iostream>
#include <ctime>
#include <cstring>
#include <algorithm>
#include <climits>

#include "mpi.h"

using namespace std;

//Ввод данных
void EnterInfo(int num_proc, int rank, int& num_vert, int& num_edge, int& start, bool& writing) {
	if (rank != 0)
		return;
	//Ввод числа вершин
	cout << "I'm proccess " << rank << "! Enter number of vertexes: ";
	cin >> num_vert;
	//Ввод числа ребер
	cout << "I'm proccess " << rank << "! Enter number of edges: ";
	cin >> num_edge;
	//Ввод стартовой вершины
	cout << "I'm proccess " << rank << "! Enter start vertex: ";
	cin >> start;
	start--;
	//Ввод инфорцмации о необходимоти вывода
	cout << "I'm proccess " << rank << "! Enter 1, if information sould be displayed, else enter any value: ";
	int tmp;
	cin >> tmp;
	if (tmp == 1)
		writing = true;
	else
		writing = false;
}

//Генерация матрицы смежности графа
unsigned char* GenerateMatrix(int num_proc, int rank, int num_vert, int num_edge, bool writing) {
	if (rank != 0)
		return NULL;
	//Выделение памяти для матрицы
	unsigned char* matrix = new unsigned char[num_vert*num_vert];
	int full_edge = num_proc * (num_proc - 1) / 2;
	//Если ребер не очень много, то будем случайно добавлять ребра, иначе - удалять
	if (num_edge <= full_edge / 2) {
		memset(matrix, 0, num_vert*num_vert * sizeof(unsigned char));
		//Выполняем для каждого ребра
		for (int i = 0; i < num_edge; ++i) {
			//Случайным образом подбираем подходящие вершины
			int x = 0, y = 0;
			while (matrix[y*num_vert + x] != 0 || x == y) {
				x = rand() % num_vert;
				y = rand() % num_vert;
			}
			//Генерируем вес ребра и заносим в матрицу
			unsigned char val = rand() % 200 + 1;
			matrix[y*num_vert + x] = val;
			matrix[x*num_vert + y] = val;
		}
	}
	else {
		//Звполняем полный граф
		for (int i = 0; i < num_vert; ++i)
			matrix[i*num_vert + i] = 0;
		for (int y = 0; y < num_vert; ++y) {
			for (int x = y+1; x < num_vert; ++x) {
				unsigned char val = rand() % 200 + 1;
				matrix[y*num_vert + x] = val;
				matrix[x*num_vert + y] = val;
			}
		}
		//Выполняем для каждого ребра
		for (int i = 0; i < full_edge - num_edge; ++i) {
			//Случайным образом подбираем подходящие вершины
			int x = 0, y = 0;
			while (matrix[y*num_vert + x] != 0 || x == y) {
				x = rand() % num_vert;
				y = rand() % num_vert;
			}
			//Генерируем вес ребра и заносим в матрицу
			unsigned char val = 0;
			matrix[y*num_vert + x] = val;
			matrix[x*num_vert + y] = val;
		}
	}
	//Вывод матрицы смежности при необходимости
	if (writing) {
		cout << "I'm procces " << rank << "! I generate matrix:";
		for (int y = 0; y < num_vert; ++y) {
			cout << "\n\t";
			for (int x = 0; x < num_vert; ++x)
				cout << (int)matrix[y*num_vert + x] << "\t";
		}
		cout << endl;
	}
	return matrix;
}

//Пересылка данных другим процессам
unsigned char* Transfer(int num_proc, int rank, unsigned char* matrix, int& num_vert, bool& writing, int& start, int& left_border, int& size_border) {
	//Пересылаем число вершин, стартовую вершину и необходимость вывода информации
	MPI_Bcast(&num_vert, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&writing, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
	MPI_Bcast(&start, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//Каждый процесс находит начало своей области и ее размер
	left_border = (num_vert / num_proc)*rank + min(rank, num_vert % num_proc);
	size_border = num_vert / num_proc + ((rank < num_vert % num_proc) ? 1 : 0);
	//Генерируем массив размеров и смещений для функции Scatterv путем сбора полученных величин от каждого процесса
	int* begin_area = NULL;
	int* size_area = NULL;
	if (rank == 0) {
		begin_area = new int[num_proc];
		size_area = new int[num_proc];
	}
	left_border *= num_vert;
	size_border *= num_vert;
	MPI_Gather(&left_border, 1, MPI_INT, begin_area, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&size_border, 1, MPI_INT, size_area, 1, MPI_INT, 0, MPI_COMM_WORLD);
	left_border /= num_vert;
	size_border /= num_vert;
	//Раздаем рабочую область каждому процессу
	unsigned char* work_area = new unsigned char[size_border*num_vert];
	MPI_Scatterv(matrix, size_area, begin_area, MPI_UNSIGNED_CHAR, 
		         work_area, size_border*num_vert, MPI_UNSIGNED_CHAR, 
		         0, MPI_COMM_WORLD);
	//Каждый процесс информирует о диапазоне своей области при необходимости
	if (writing) {
		if (size_border>0)
			cout << "I'm proccess " << rank << "! My borders is from " << left_border + 1 << " to " << left_border + size_border << endl;
		else
			cout << "I'm proccess " << rank << "! I don't have work area" << endl;
	}
	//Очищаем ненужную память
	if (begin_area != NULL)
		delete[] begin_area;
	if (size_area != NULL)
		delete[] size_area;
	return work_area;
}

//Алгоритм Дейкстры
unsigned int* DijkstraAlgorithm(int num_proc, int rank, int num_vert, int start, bool writing, int left_border, int size_border, unsigned char* work_area) {
	//Выделяем память под массив расстояний до стартовой вершины и заполняем их максимальными значениями (пути нет)
	unsigned int* distance = new unsigned int[size_border];
	memset(distance, UCHAR_MAX, size_border * sizeof(unsigned int));
	//Выделяем память под множество пройденных вершин и заполняем их нулями (пройденных вершин нет)
	bool* was = new bool[size_border];
	memset(was, false, size_border * sizeof(bool));
	//Объявляем текущую вершину как стартовую и расстояние до нее как 0
	int cur_vert = start;
	unsigned int cur_dist = 0;
	//Процесс, в рабочей области которого находится стартовая вершина, помечает у себя 0 расстояние до нее
	if (left_border <= cur_vert && cur_vert < left_border + size_border)
		distance[cur_vert - left_border] = 0;
	//Выделение памяти для масива минимальных дистанций, полученных от каждого процесса в 0
	unsigned int* glob_min = NULL;
	if (rank == 0)
		glob_min = new unsigned int[num_proc];
	if (writing)
		MPI_Barrier(MPI_COMM_WORLD);
	//Выполняем цикл, пока есть текущая вершина
	while (cur_vert != -1) {
		//При еобходимости 0 процесс объявляет, какую вершину смотрят на данный момент
		if (writing && rank == 0)
			cout << "I'm proccess " << rank << "! Current vertex is " << cur_vert + 1 << ", distance to it is " << cur_dist << endl;
		//Процесс, в рабочей области которого находится стартовая вершина, помечает у себя, что она уже пройдена
		if (left_border <= cur_vert && cur_vert < left_border + size_border)
			was[cur_vert - left_border] = true;
		//Обновляем дистанции относительно текущей вершины
		for (int i = 0; i < size_border; ++i) {
			if (!was[i] && work_area[i*num_vert + cur_vert] != 0 && distance[i] > cur_dist + work_area[i*num_vert + cur_vert])
				distance[i] = cur_dist + work_area[i*num_vert + cur_vert];
		}
		//Найдем непройденную вершину с минимальным расстоянием от стартовой
		unsigned int min_dist = UINT_MAX;
		int min_vert = -1;
		for (int i = 0; i < size_border; ++i) {
			if (!was[i] && min_dist > distance[i]) {
				min_vert = i;
				min_dist = distance[i];
			}
		}
		//Передаем найденные расстояния 0 процессу
		MPI_Gather(&min_dist, 1, MPI_UNSIGNED_LONG, glob_min, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
		//Нулевой процесс находит минимальное расстояние и сообщает всем процесс, от которого его получил (процесс-источник)
		int finder = 0;
		if (rank == 0) {
			for (int i = 0; i < num_proc; ++i) {
				if (glob_min[i] < glob_min[finder])
					finder = i;
			}
		}
		MPI_Bcast(&finder, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//Процесс-источник объявляет всем новую текущую вершину и расстояние до нее
		if (rank == finder) {
			cur_vert = min_vert + left_border;
			cur_dist = min_dist;
		}
		MPI_Bcast(&cur_vert, 1, MPI_INT, finder, MPI_COMM_WORLD);
		MPI_Bcast(&cur_dist, 1, MPI_UNSIGNED_LONG, finder, MPI_COMM_WORLD);
	}
	if (was != NULL)
		delete[] was;
	return distance;
}

//Складываем результаты работы в одном процессе
unsigned int* FindResults(int num_proc, int rank, int num_vert, int start, bool writing, int left_border, int size_border, unsigned int* distance) {
	//Выделяем память под массив результатов, а также массивы размеров и смещений для функции Gatherv
	unsigned int* results = NULL;
	int* size_area = NULL;
	int* begin_area = NULL;
	if (rank == 0) {
		results = new unsigned int[num_vert];
		size_area = new int[num_proc];
		begin_area = new int[num_proc];
	}
	//Собираем информацию о размерах и смещениях от каждого процесса
	MPI_Gather(&left_border, 1, MPI_INT, begin_area, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&size_border, 1, MPI_INT, size_area, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//Собираем результаты
	MPI_Gatherv(distance, size_border, MPI_UNSIGNED_LONG, results, size_area, begin_area, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	//Выводим информацию о реультатах
	if (rank == 0 && writing) {
		cout << "I'm proccess " << rank << "! We find distances from " << start + 1 << " vertex:";
		for (int i = 0; i < num_vert; ++i) {
			if (results[i] != UINT_MAX)
				cout << "\n\tTo " << i + 1 << ": " << results[i];
			else
				cout << "\n\tTo " << i + 1 << ": " << "No distance";
		}
		cout << endl;
	}
	//Очищаем ненужную память
	if (begin_area != NULL)
		delete[] begin_area;
	if (size_area != NULL)
		delete[] size_area;
	return results;
}

//Проврка корректности (алгоритм тот же)
void Check(int num_vect, int start, bool writing, unsigned char* matrix, unsigned int* results) {
	double start_time = MPI_Wtime();
	unsigned int* distance = new unsigned int[num_vect];
	memset(distance, UCHAR_MAX, num_vect * sizeof(unsigned int));
	bool* was = new bool[num_vect];
	memset(was, false, num_vect*sizeof(bool));
	int cur_vect = start;
	distance[start] = 0;
	while (cur_vect != -1) {
		was[cur_vect] = 1;
		for (int i = 0; i < num_vect; ++i) {
			if (!was[i] && matrix[cur_vect*num_vect + i]!=0 && distance[i] > distance[cur_vect] + matrix[cur_vect*num_vect + i])
				distance[i] = distance[cur_vect] + matrix[cur_vect * num_vect + i];
		}
		cur_vect = -1;
		unsigned int min_dist = UINT_MAX;
		for (int i = 0; i < num_vect; ++i) {
			if (!was[i] && min_dist > distance[i]) {
				cur_vect = i;
				min_dist = distance[i];
			}
		}
	}
	if (writing) {
		cout << "Check correctness:";
		for (int i = 0; i < num_vect; ++i)
			cout << "\n\tTo " << i + 1 << ": " << distance[i];
		cout << endl;
	}
	double finish_time = MPI_Wtime();
	bool check = true;
	for (int i = 0; i < num_vect; ++i) {
		if (distance[i] != results[i])
			check = false;
	}
	if (check)
		cout << "Check finished with success" << endl;
	else
		cout << "Check finished with failure" << endl;
	cout << "Time: " << finish_time - start_time << endl;
}

int main(int argc, char** argv) {
	cout << fixed;
	cout.precision(10);
	srand(time(0));
	int num_proc, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	int num_vert, num_edge;
	bool writing;
	int start = -1;
	EnterInfo(num_proc, rank, num_vert, num_edge, start, writing);
	unsigned char* matrix = GenerateMatrix(num_proc, rank, num_vert, num_edge, writing);
	double start_time = MPI_Wtime();
	int left_border;
	int size_border;
	unsigned char* work_area = Transfer(num_proc, rank, matrix, num_vert, writing, start, left_border, size_border);
	unsigned int* distance = DijkstraAlgorithm(num_proc, rank, num_vert, start, writing, left_border, size_border, work_area);
	unsigned int* results = FindResults(num_proc, rank, num_vert, start, writing, left_border, size_border, distance);
	double finish_time = MPI_Wtime();
	if (rank == 0) 
		cout << "I'm proccess " << rank << "! Time: " << finish_time - start_time<<endl;
	if (rank == 0)
		Check(num_vert, start, writing, matrix, results);
	if (matrix != NULL)
		delete[] matrix;
	if (distance!=NULL)
		delete[] distance;
	if (results != NULL)
		delete[] results;
	MPI_Finalize();
	return 0;
}