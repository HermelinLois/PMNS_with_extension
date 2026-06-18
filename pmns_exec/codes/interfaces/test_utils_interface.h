#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>
#include "../params/pmns_params.h"

void rand_field_element(mp_limb_t out[EXTENSION_DEGREE][N_LIMBS], gmp_randstate_t state);

void gen_random_pmns(int64_t out[DEGREE]);

void print_pol(int64_t *P);

void check_equality(int64_t P[DEGREE], int64_t C[DEGREE], const char* method);

void reset_polynomial(int64_t polynomial[DEGREE]);

#endif
