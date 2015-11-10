#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include "random.h"
#include "mpi_constants.h"
#include "sample_sort.h"

double comm_time;

int compare (const void * a, const void * b)
{
	const int *ia = (const int *)a; // casting pointer types
	const int *ib = (const int *)b;
	return *ia  - *ib;
}

int *
pick_pivots(int *to_sort, int size)
{
	int i, *pivots, chunk_size;

	chunk_size = size/p;

	pivots = (int *) malloc(sizeof(int)*(p-1));

	for (i = 1; i < p; i++) {
		pivots[i-1] = to_sort[chunk_size*i];
	}

	return pivots;
}

void
partition(int *to_partition, int size, int *pivots, int npivots, int **displs,
    int **counts)
{
	int i, count, displ;

	*displs = (int *) malloc(sizeof(int)*(npivots+1));
	*counts = (int *) malloc(sizeof(int)*(npivots+1));

	displ = 0;
	count = 0;
	for (i = 0; (i < npivots && displ + count < size); i++) {
		(*displs)[i] = displ;
		while (displ + count < size &&
		    to_partition[displ + count] <= pivots[i])
			count++;
		(*counts)[i] = count;
		displ += count;
		count = 0;
	}
	if (displ != size) {
		/* Room left at the end */
		(*displs)[npivots] = displ;
		(*counts)[npivots] = size - displ;
	}
	else {
		/* Not room left at end, all point at last element and send 0 */
		while (i <= npivots + 1) {
			(*displs)[i] = size-1;
			(*counts)[i] = 0;
			i++;
		}
	}
}

int
sample_sort(int n, int *to_sort, int **sorted, int *sorted_size)
{
	int err, i;
	int *lpivots = NULL, n_lpivots, *all_lpivots = NULL, n_all_lpivots,
	    *gpivots = NULL, n_gpivots;
	int *sdispls = NULL, *sendcounts = NULL, *rdispls = NULL, *recvcounts,
	    *recvbuf = NULL, rsize;
	struct timeval t1, t2;

	comm_time = 0.0;

	/* Local sort */
	qsort(to_sort, n/p, sizeof(int), compare);

	/* Pick local pivots */
	n_lpivots = p - 1;
	lpivots = pick_pivots(to_sort, n/p);


	/* Gather local pivots from all processes for sorting */
	n_all_lpivots = n_lpivots * p;
	all_lpivots = (int *) malloc(sizeof(int)*n_all_lpivots);

	gettimeofday(&t1, NULL);
	err = MPI_Allgather(lpivots, n_lpivots, MPI_INT, all_lpivots,
	    n_lpivots, MPI_INT, MPI_COMM_WORLD);
	gettimeofday(&t2, NULL);
	if (err)
		goto out;

	comm_time += (t2.tv_sec-t1.tv_sec)*1000 +
	    ((double) t2.tv_usec-t1.tv_usec)/1000;

	/* Sort local pivots from all processes */
	qsort(all_lpivots, n_all_lpivots, sizeof(int), compare);


	/* Pick global pivots */
	n_gpivots = p - 1;
	gpivots = pick_pivots(all_lpivots, n_all_lpivots);

	/* Partition with global pivots */
	gettimeofday(&t1, NULL);
	partition(to_sort, n/p, gpivots, n_gpivots, &sdispls, &sendcounts);

	/* Tell other processes how much data to expect */
	rdispls = (int *) malloc(sizeof(int)*p);
	recvcounts = (int *) malloc(sizeof(int)*p);
	err = MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT,
	    MPI_COMM_WORLD);
	if (err)
		goto out;


	/* Calculate receiving displacements */
	rdispls[0]=0;
	rsize = recvcounts[0];
	for(i = 1; i < p; i++) {
		rdispls[i] = recvcounts[i-1] + rdispls[i-1];
		rsize += recvcounts[i];
	}


	/* Shuffle chunks of data to appropriate processes */
	recvbuf = (int *) malloc(sizeof(int)*rsize);
	err = MPI_Alltoallv(to_sort, sendcounts, sdispls, MPI_INT, recvbuf,
	    recvcounts, rdispls, MPI_INT, MPI_COMM_WORLD);
	gettimeofday(&t2, NULL);
	if (err)
		goto out;

	comm_time += (t2.tv_sec-t1.tv_sec)*1000 +
	    ((double) t2.tv_usec-t1.tv_usec)/1000;


	/* Final local sort */
	qsort(recvbuf, rsize, sizeof(int), compare);

	*sorted = recvbuf;
	*sorted_size = rsize;


out:
	/* Cleanup */
	if (lpivots)
		free(lpivots);
	if (all_lpivots)
		free(all_lpivots);
	if (gpivots)
		free(gpivots);
	if (sdispls)
		free(sdispls);
	if (rdispls)
		free(rdispls);
	if (sendcounts)
		free(sendcounts);
	if (recvcounts)
		free(recvcounts);

	return err;
}

int
run_sample_sort(int **sorted, int *sorted_size)
{
	int err;
	int *arr = NULL;

	err = generate_random_series(n, &arr);
	if (err)
		goto out;

	err = sample_sort(n, arr, sorted, sorted_size);
	if (err)
		goto out;

out:
	if (arr)
		free(arr);

	return err;
}
