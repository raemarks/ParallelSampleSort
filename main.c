#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include "mpi_constants.h"
#include "sample_sort.h"
#include "random.h"

#define RED "\x1B[31m"
#define GREEN "\x1B[32m"
#define BLUE "\x1B[34m"
#define YELLOW "\x1B[33m"
#define WHITE "\x1B[37m"

int test_sample_sort();
int test_prefix_operation();
int test_all_to_all();
void test_partition();

int main(int argc, char *argv[])
{
	P=9973;
	B=19;
	A=7;
	seed=3557;
#if 0
	struct timeval t1, t2;
	int *result, i, it, err, iterations = 10;
	double myruntime, *runtimes, iteration_sum;
	double sum = 0.0;
	FILE *outfile;
	char *filename;
#endif

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&p);

	assert(p>=1);

	test_sample_sort();
	//test_prefix_operation();
	//test_partition();

	MPI_Finalize();

	return 0;
}

int test_all_to_all()
{
	return 0;
}

int test_sample_sort()
{
	int expected[256], *arr;
	int err, i;
	int *sorted, sorted_size;
	int *recvcounts = NULL, *rdispls = NULL, *recvbuf = NULL, rsize;

	n = 256;

	err = generate_random_series_gathered(n, &arr);

	if (rank == 0) {
		memcpy(expected, arr, n*sizeof(int));
		qsort(expected, n, sizeof(int), compare);
		free(arr);
	}

	err = run_sample_sort(&sorted, &sorted_size);
	if (err) {
		printf("Got an error running sample sort: %d\n", err);
		goto out;
	}

	rdispls = (int *) malloc(sizeof(int)*p);
	recvcounts = (int *) malloc(sizeof(int)*p);
	recvbuf = (int *) malloc(sizeof(int)*n);

	err = MPI_Gather(&sorted_size, 1, MPI_INT, recvcounts, 1, MPI_INT, 0,
	    MPI_COMM_WORLD);

	if (rank == 0) {
		rdispls[0]=0;
		rsize = recvcounts[0];
		for(i = 1; i < p; i++) {
			rdispls[i] = recvcounts[i-1] + rdispls[i-1];
			rsize += recvcounts[i];
		}
	}
	if (rank == 0 && rsize != n)
		printf("Not going to receive n! Receiving this instead: %d\n",
		    rsize);

	err = MPI_Gatherv(sorted, sorted_size, MPI_INT, recvbuf, recvcounts,
	    rdispls, MPI_INT, 0, MPI_COMM_WORLD);
	if (err) {
		printf("Error calling gatherv: %d \n", err);
		goto out;
	}

	if (rank == 0) {
		if (memcmp(expected, recvbuf, n) == 0) {
			printf("%sSuccess!\n%s\n", GREEN, WHITE);
		}
		else {
			printf("%sFAILURE! \n%sExpected:\n", RED, WHITE);
			for (i = 0; i < n; i++) {
				printf("%d ", expected[i]);
			}
			printf("\nActual:\n");
			for (i = 0; i < n; i++) {
				printf("%d ", recvbuf[i]);
			}
			putchar('\n');
		}
	}

out:
	if (rdispls)
		free(rdispls);
	if (recvcounts)
		free(recvcounts);
	if (recvbuf)
		free(recvbuf);
	if (sorted)
		free(sorted);

	return 0;
}

int
test_prefix_operation()
{
	int i, err, prev, res;
	int expected[2048];
	int *result;
	int64_t t;

	n = 2048;
	prev = seed;
	expected[0] = seed;
	for (i = 1; i < n; i++) {
		t = (((int64_t)A)*((int64_t)prev) + (uint64_t) B) % (uint64_t) P;
		expected[i] = (int) t;
		prev = expected[i];
	}

	set_constants(n, p, rank, A, B, P, seed);
	//err = generate_fake_random(n, &result);
	err = generate_random_series_gathered(n, &result);

	if (rank == 0) {
		if ((res = memcmp(expected, result, n)) == 0) {
			printf("%sSuccess!\n%s\n", GREEN, WHITE);
		}
		else {
			printf("%sFAILURE! \n%sExpected:\n", RED, WHITE);
			for (i = 0; i < n; i++) {
				printf("%d ", expected[i]);
			}
			printf("\nActual:\n");
			for (i = 0; i < n; i++) {
				printf("%d ", result[i]);
			}
			putchar('\n');
		}
		free(result);
	}

	return res;
}

void
test_partition()
{
	int i;
	int a1[12] = {0, 3, 6, 9, 11, 13, 15, 15, 17, 20, 25, 92};
	int pivots1[3] = {9, 16, 21};
	int edispls1[4] = {0, 4, 8, 10};
	int ecounts1[4] = {4, 4, 2, 2};
	int *displs1, *counts1;

	int a2[12] = {1, 2, 2, 2, 3, 3, 4, 16, 17, 24, 24, 15};
	int pivots2[3] = {0, 2, 16};
	int edispls2[4] = {0, 0, 4, 8};
	int ecounts2[4] = {0, 4, 4, 4};
	int *displs2, *counts2;

	int a3[12] = {1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 4, 5};
	int pivots3[3] = {1, 6, 7};
	int edispls3[4] = {0, 3, 11, 11};
	int ecounts3[4] = {3, 9, 0, 0};
	int *displs3, *counts3;

	if (rank != 0)
		return;

	partition(a1, 12, pivots1, 3, &displs1, &counts1);
	if (memcmp(displs1, edispls1, 4*sizeof(int)) == 0) {
		printf("\n%sDISPLS1 SUCCESS%s\n", GREEN, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", edispls1[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", displs1[i]);
		}
		putchar('\n');
	}
	else {
		printf("\n%sDISPLS1 FAILURE! \n%sExpected:\n", RED, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", edispls1[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", displs1[i]);
		}
		putchar('\n');
	}
	if (memcmp(counts1, ecounts1, 4*sizeof(int)) == 0) {
		printf("\n%sCOUNTS1 SUCCESS%s\n", GREEN, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", ecounts1[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", counts1[i]);
		}
		putchar('\n');
	}
	else {
		printf("\n%sCOUNTS1 FAILURE! \n%sExpected:\n", RED, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", ecounts1[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", counts1[i]);
		}
		putchar('\n');
	}
	free(displs1);
	free(counts1);

	partition(a2, 12, pivots2, 3, &displs2, &counts2);
	if (memcmp(displs2, edispls2, 4*sizeof(int)) == 0) {
		printf("\n%sDISPLS2 SUCCESS%s\n", GREEN, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", edispls2[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", displs2[i]);
		}
		putchar('\n');
	}
	else {
		printf("\n%sDISPLS2 FAILURE! \n%sExpected:\n", RED, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", edispls2[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", displs2[i]);
		}
		putchar('\n');
	}
	if (memcmp(counts2, ecounts2, 4*sizeof(int)) == 0) {
		printf("\n%sCOUNTS2 SUCCESS%s\n", GREEN, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", ecounts2[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", counts2[i]);
		}
		putchar('\n');
	}
	else {
		printf("\n%sCOUNTS2 FAILURE! \n%sExpected:\n", RED, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", ecounts2[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", counts2[i]);
		}
		putchar('\n');
	}
	free(displs2);
	free(counts2);

	partition(a3, 12, pivots3, 3, &displs3, &counts3);
	if (memcmp(displs3, edispls3, 4*sizeof(int)) == 0) {
		printf("\n%sDISPLS3 SUCCESS%s\n", GREEN, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", edispls3[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", displs3[i]);
		}
		putchar('\n');
	}
	else {
		printf("\n%sDISPLS3 FAILURE! \n%sExpected:\n", RED, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", edispls3[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", displs3[i]);
		}
		putchar('\n');
	}
	if (memcmp(counts3, ecounts3, 4*sizeof(int)) == 0) {
		printf("\n%sCOUNTS3 SUCCESS%s\n", GREEN, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", ecounts3[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", counts3[i]);
		}
		putchar('\n');
	}
	else {
		printf("\n%sCOUNTS3 FAILURE! \n%sExpected:\n", RED, WHITE);
		for (i = 0; i < 4; i++) {
			printf("%d ", ecounts3[i]);
		}
		printf("\nActual:\n");
		for (i = 0; i < 4; i++) {
			printf("%d ", counts3[i]);
		}
		putchar('\n');
	}
	free(displs3);
	free(counts3);
}
