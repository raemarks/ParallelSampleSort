#ifndef _RANDOM_H_
#define _RANDOM_H_
#include <stdint.h>

int
generate_random_series(int n, int **series_out);

int
generate_random_series_gathered(int n, int **series_out);

int
generate_fake_random(int n, int **result);

void
set_constants(int _n, int _p, int _rank, int _A, int _B, int _P, int _seed);

void
test_matrix_power();

#endif
