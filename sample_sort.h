#ifndef _SAMPLE_SORT_H_
#define _SAMPLE_SORT_H_

int compare (const void * a, const void * b);
void
partition(int *to_partition, int size, int *pivots, int npivots, int **displs,
    int **counts);

int
sample_sort(int n, int *to_sort, int **sorted, int *sorted_size);

int
run_sample_sort(int **sorted, int *sorted_size);

#endif
