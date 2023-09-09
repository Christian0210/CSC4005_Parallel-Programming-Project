#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <iomanip>
using namespace std;

const int N = 20;
const int MAX_NUM = 1000;

int OddEvenSort(int local[], const int m, const int my_rank, int p, MPI_Comm comm)
{
	int i, j, send_temp = 0, recv_temp = 10001, rrank = (my_rank + 1) % p, lrank = (my_rank + p - 1) % p;

	for (i = 0; i < p * m; i++)
	{
		if (i % 2)
		{
			for (j = m - 2; j > 0; j -= 2)
				if (local[j] < local[j - 1])
					swap(local[j], local[j - 1]);

			if (my_rank != 0)
			{
				send_temp = local[0];
				MPI_Send(&send_temp, 1, MPI_INT, lrank, 0, comm);
				MPI_Recv(&recv_temp, 1, MPI_INT, lrank, 0, comm, MPI_STATUS_IGNORE);
				if (recv_temp > local[0])
					local[0] = recv_temp;
			}
			if (my_rank != p - 1) {
				send_temp = local[m - 1];
				MPI_Recv(&recv_temp, 1, MPI_INT, rrank, 0, comm, MPI_STATUS_IGNORE);
				MPI_Send(&send_temp, 1, MPI_INT, rrank, 0, comm);
				if (recv_temp < local[m - 1])
					local[m - 1] = recv_temp;
			}
		}
		else
		{
			for (j = m - 1; j > 0; j -= 2)
				if (local[j] < local[j - 1])
					swap(local[j], local[j - 1]);
		}
	}
	return 0;
}

int main(int argc, char* argv[])
{
	MPI_Comm comm;
	int i, n = N, m, rank, proc;

	int* Arr = 0;
	clock_t start, end;
	double time;

	start = clock();
	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc);

	m = n / proc;

	int* num = (int*)malloc(sizeof(int) * n);
	int* local = (int*)malloc(sizeof(int) * m);

	if (!rank)
	{
		srand(clock());
		for (i = 0; i < n; i++)
			num[i] = rand() % MAX_NUM;
		printf("Name: Yuncong Cui\nStudent ID: 118010045\nHomework 1, Odd Even Sort, MPI Implementation\n\n---------------Before Odd Even Sort---------------\nArray: ");
		for (i = 0; i < n; i++)
			printf("%d ", num[i]);
	}

	MPI_Scatter(num, m, MPI_INT, local, m, MPI_INT, 0, comm);
	OddEvenSort(local, m, rank, proc, comm);
	MPI_Gather(local, m, MPI_INT, num, m, MPI_INT, 0, comm);

	end = clock();
	time = (double)(end - start) / CLOCKS_PER_SEC;

	if (!rank)
	{
		printf("\n\n---------------After Odd Even Sort----------------\nArray: ");
		for (i = 0; i < n; i++)
			printf("%d ", num[i]);
		printf("\n\n---------------------Analysis---------------------\nProcessor: %d\nThe number of elements: %d\nRuntime: %lf seconds\n\n", proc, n, time);
	}

	MPI_Finalize();

	return 0;
}