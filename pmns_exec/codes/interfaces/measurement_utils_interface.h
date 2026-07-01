#ifndef MEASUREMENT_UTILS_H
#define MEASUREMENT_UTILS_H

#include <stdint.h>
#include <stdlib.h>

#define N_BENCH_TESTS 501
#define N_BENCH_SAMPLES 1001

#ifdef __cplusplus
extern "C" {
#endif

void quicksort(uint64_t* t, int n);

uint64_t *quartiles(uint64_t *tab, uint64_t size);

uint64_t cpucyclesStart(void);

uint64_t cpucyclesStop(void);

uint64_t cycles_diff(uint64_t start, uint64_t stop);

#ifdef __cplusplus
}
#endif

#endif
