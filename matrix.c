#include "matrix.h"
#include <mpi.h>
#include <stdio.h>
#include <stdint.h>

int64_t *
get_cell(int64_t *matrix, int x, int y)
{
	return &(matrix[x*2 + y]);
}

void
mult_matrices(int64_t *matrix1, int64_t *matrix2, int64_t *out_matrix)
{
	*(get_cell(out_matrix, 0, 0)) =
	    (*get_cell(matrix1, 0, 0))*(*get_cell(matrix2, 0, 0)) +
	    (*get_cell(matrix1, 0, 1))*(*get_cell(matrix2, 1, 0));
	*(get_cell(out_matrix, 0, 1)) =
	    (*get_cell(matrix1, 0, 0))*(*get_cell(matrix2, 0, 1)) +
	    (*get_cell(matrix1, 0, 1))*(*get_cell(matrix2, 1, 1));
	*(get_cell(out_matrix, 1, 0)) =
	    (*get_cell(matrix1, 1, 0))*(*get_cell(matrix2, 0, 0)) +
	    (*get_cell(matrix1, 1, 1))*(*get_cell(matrix2, 1, 0));
	*(get_cell(out_matrix, 1, 1)) =
	    (*get_cell(matrix1, 1, 0))*(*get_cell(matrix2, 0, 1)) +
	    (*get_cell(matrix1, 1, 1))*(*get_cell(matrix2, 1, 1));
}

void
matrix_mod_p(int64_t matrix[4], int p)
{
	matrix[0] %= p;
	matrix[1] %= p;
	matrix[2] %= p;
	matrix[3] %= p;
}

void
mult_vector_matrix(int64_t *vector, int64_t *matrix, int64_t *out_vector)
{
	*(get_cell(out_vector, 0, 0)) =
	    (*get_cell(vector, 0, 0))*(*get_cell(matrix, 0, 0)) +
	    (*get_cell(vector, 0, 1))*(*get_cell(matrix, 1, 0));
	*(get_cell(out_vector, 0, 1)) =
	    (*get_cell(vector, 0, 0))*(*get_cell(matrix, 0, 1)) +
	    (*get_cell(vector, 0, 1))*(*get_cell(matrix, 1, 1));

}

void test_mult_matrices() {
	int64_t matrix1[4] = {1,2,3,4}, matrix2[4]={5,6,7,8}, matrix3[4]={19,22,43,50},
	    matrix4[4];


	mult_matrices(matrix1, matrix2, matrix4);

	printf("Expected:\n%3ld %3ld\n%3ld %3ld\n\n", matrix3[0], matrix3[1],
	    matrix3[2], matrix3[3]);

	printf("Result:\n%3ld %3ld\n%3ld %3ld\n\n", matrix4[0], matrix4[1],
	    matrix4[2], matrix4[3]);

}

void test_mult_vector_matrix() {
	int64_t vector1[2] = {2,1}, matrix1[4]={5,0,7,1}, vector2[2]={17,1}, vector3[2];


	mult_vector_matrix(vector1, matrix1, vector3);

	printf("Expected:\n%3ld %3ld\n\n", vector2[0], vector2[1]);
	printf("Result:\n%3ld %3ld\n\n", vector3[0], vector3[1]);
}

int
swap_matrices(int64_t send_matrix[4], int64_t recv_matrix[4], int mate)
{
	MPI_Status status;

	return MPI_Sendrecv(send_matrix, 4, MPI_LONG_LONG_INT, mate, 0, recv_matrix, 4,
	    MPI_LONG_LONG_INT, mate, 0, MPI_COMM_WORLD, &status);
}
