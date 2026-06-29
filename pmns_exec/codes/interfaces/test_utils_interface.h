#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>
#include "../params/pmns_params.h"

void rand_field_element(mp_limb_t out[EXTENSION_DEGREE][N_LIMBS], gmp_randstate_t state);

void print_pol(int degree, int64_t P[degree]);

void check_equality(int degree, int64_t P[degree], int64_t C[degree], const char* method);

void reset_polynomial(int degree, int64_t polynomial[degree]);

#endif
