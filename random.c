#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "random.h"
#include "matrix.h"

/* For parallel prefix */
int64_t local_product[4];
int64_t global_product[4];

/* Input size */
int n;

/* MPI constants */
int p;
int rank;

/* Big primes */
int64_t A;
int64_t B;
int64_t P;
int64_t seed;

void
set_constants(int _n, int _p, int _rank, int _A, int _B, int _P, int _seed)
{
	n = _n;
	p = _p;
	rank = _rank;
	A = _A;
	B = _B;
	P = _P;
	seed = _seed;
}

void
compute_matrix_power(int power, int64_t matrix[4])
{
	int64_t matrix1[4] = {A, 0, B, 1};
	int64_t matrix2[4] = {A, 0, B, 1};
	int64_t matrix3[4] = {A, 0, B, 1};
	int64_t *m1, *m2, *m3, *temp;
	int i;

	if (power == 0) {
		matrix[0] = 1;
		matrix[1] = 0;
		matrix[2] = 0;
		matrix[3] = 1;
		return;
	}

	m1 = matrix1;
	m2 = matrix2;
	m3 = matrix3;
	for (i = 0; i < power - 1; i++) {
		mult_matrices(m1, m2, m3);
		matrix_mod_p(m3, P);
		temp = m3;
		m3 = m2;
		m2 = temp;
	}
	memcpy(matrix, m2, sizeof(int64_t)*4);
}

void
compute_initial_global_product()
{
	/* Shift power down by one on the first matrix. If we don't do this,
	 * then every rank will need to get the last element from the previous
	 * rank. Shifting the first element down by one will shift evey
	 * subsequent element down by one because of the prefix. This is a
	 * little perforamnce cheat. */
	if (rank == 0)
		compute_matrix_power(n/p-1, global_product);
	else
		compute_matrix_power(n/p, global_product);

	local_product[0] = A;
	local_product[1] = 0;
	local_product[2] = B;
	local_product[3] = 1;
}

int
parallel_prefix()
{
	int mate, i, err;
	//Remote global product.
	int64_t global_product_remote[4], matrix[4];

	for (i = 0; i < log2((double) p); i++) {
		mate = rank ^ (1 << i);

		err = swap_matrices(global_product, global_product_remote, mate);
		if (err) {
			return err;
		}

		mult_matrices(global_product, global_product_remote, matrix);
		matrix_mod_p(matrix, P);
		memcpy(global_product, matrix, sizeof(int64_t)*4);

		if (mate < rank) {
			mult_matrices(local_product, global_product_remote,
			    matrix);
			matrix_mod_p(matrix, P);
			memcpy(local_product, matrix, sizeof(int64_t)*4);
		}
	}

	return 0;
}

int
generate_random_series_gathered(int n, int **series_out) {
	int err = 0;
	int *local_series;

	if (rank == 0)
		*series_out = (int *) malloc(sizeof(int)*n);

	err = generate_random_series(n, &local_series);
	if (err) return err;

	err = MPI_Gather(local_series, n/p, MPI_INT, *series_out, n/p, MPI_INT,
	    0, MPI_COMM_WORLD);

	free(local_series);

	return err;
}

int
generate_random_series(int n, int **series_out)
{
	int err = 0, i, prev;
	int64_t init_vector[2] = {seed, 1}, result_vector[2] = {0};
	int *series;
	/*
	 * Algorithm:
	 *
	 * 1. Distribute job amongst procs. i.e. decide which range each rank
	 * will be responsible for. Broadcast the seed.
	 *
	 * 2. Compute local product of matrices.
	 *
	 * 3. Compute parallel prefix.
	 *
	 * 4. Compute local random numbers.
	 *
	 * 5. Reduce.
	 */

	/* 1. */
	if (err)
		return err;

	/* 2. */
	compute_initial_global_product();

	/* 3. */
	err = parallel_prefix();
	if (err)
		return err;

	/* 4. */
	*series_out = (int *) malloc(sizeof(int)*n/p);
	series = *series_out;

	if (rank == 0)
		series[0] = seed;
	else {
		mult_vector_matrix(init_vector, local_product, result_vector);
		series[0] = result_vector[0] % P;
	}

	prev = series[0];
	for (i = 1; i < n/p; i++) {
		series[i] = (A*prev + B) % P;
		prev = series[i];
	}

	return 0;
}

int
generate_fake_random(int n, int **series_out)
{
	int64_t m[4];

	int i;
	int *series;

	MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&A, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&B, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&P, 1, MPI_INT, 0, MPI_COMM_WORLD);

	*series_out = (int *) malloc(sizeof(int)*n/p);
	series = *series_out;

	for (i = 0; i < n/p; i++) {
		int64_t init_vector[2] = {(int64_t)seed, (int64_t)1},
			result_vector[2] = {0};
		if (rank == 0 && i == 0) {
			series[0] = init_vector[0];
			continue;
		}
		compute_matrix_power((rank*n/p) + i, m);
		mult_vector_matrix(init_vector, m, result_vector);
		series[i] = result_vector[0] % P;
	}

	/* 4. */
	return 0;
}

void
test_matrix_power()
{
	int64_t m1[4] = {A, 0, B, 1}, m2[4] = {A*A, 0, B*A + B, 1}, m3[4] = {0,0,1,1};
	int64_t m[4];
	m3[0] = m2[0]*m1[0];
	m3[2] = m1[2]*m2[0] + m2[2];

	printf("Testing matrix power\n");
	compute_matrix_power(1, m);
	if (memcmp(m, m1, 4*sizeof(int64_t)) != 0) {
		printf("Power 1 wrong!\n");
	}
	compute_matrix_power(2, m);
	if (memcmp(m, m2, 4*sizeof(int64_t)) != 0) {
		printf("Power 2 wrong!\n");
	}
	compute_matrix_power(3, m);
	if (memcmp(m, m3, 4*sizeof(int64_t)) != 0) {
		printf("Power 3 wrong!\n");
	}

	compute_matrix_power(2, m);
	printf("original:\n%ld %ld\n%ld %ld\npower 2:\n%ld %ld\n%ld %ld\n", m1[0], m1[1],
	    m1[2], m1[3], m[0], m[1], m[2], m[3]);
}

void
test_parallel_prefix() {
	compute_initial_global_product();
	parallel_prefix();
}
