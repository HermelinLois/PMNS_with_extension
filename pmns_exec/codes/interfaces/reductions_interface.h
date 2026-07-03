#ifndef REDUCTION_API
#define REDUCTION_API

#include <stdint.h>
#include "../params/pmns_params.h"
#include "../params/reductions_params.h"

# include "./interfaces/reductions_interface.h"
# include "./interfaces/operations_interface.h"
# include <stdio.h>
# include <string.h>
# include <gmp.h>

/*=================================================================
                        BABAI REDUCTION FUNCTION 
=================================================================*/
# if IS_BABAI_USABLE
void reduction_babai(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]);
# endif

/*=================================================================
                    MONTGOMERY REDUCTION FUNCTION 
=================================================================*/
void reduction_montgomery_lattice(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]);

# if IS_TOEPLITZ_USABLE
void reduction_montgomery_toeplitz(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[2 * DEGREE - 1], const uint64_t sublattice_inv[2 * DEGREE - 1]);

void reduction_montgomery_toeplitz_recursive(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[2 * DEGREE - 1], const uint64_t sublattice_inv[2 * DEGREE - 1]);
# endif

# if IS_DOUBLE_SPARSE
void reduction_montgomery_linear(int64_t out[DEGREE], __int128 polynomial[DEGREE]);
# endif

/*=================================================================
                MONTGOMERY REDUCTION FOR CONVERSION 
=================================================================*/
void reduction_montgomery_mpz(mpz_t out[DEGREE], mpz_t polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]);

void reduction_montgomery_mpn(int n_limbs, mp_limb_t out[DEGREE][n_limbs], mp_limb_t pol[DEGREE][n_limbs], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]);
#endif