#ifndef _MATRIX_H_
#define _MATRIX_H_
#include <stdint.h>

void
mult_matrices(int64_t *matrix1, int64_t *matrix2, int64_t *out_matrix);

void
mult_vector_matrix(int64_t *vector, int64_t *matrix, int64_t *out_vector);

int
swap_matrices(int64_t send_matrix[4], int64_t recv_matrix[4], int mate);

void
mult_vector_matrix(int64_t *vector, int64_t *matrix, int64_t *out_vector);

void
matrix_mod_p(int64_t matrix[4], int p);

#endif
