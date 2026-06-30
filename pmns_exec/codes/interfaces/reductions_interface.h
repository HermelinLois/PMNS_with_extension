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
void reduction_babai_int128(int degree, int64_t out[degree], __int128 polynomial[degree], const int64_t sublattice[degree][degree], const int64_t sublattice_inv[degree][degree]);
# endif

/*=================================================================
                    MONTGOMERY REDUCTION FUNCTION 
=================================================================*/
void reduction_montgomery_int128(int degree, int64_t out[degree], __int128 polynomial[degree], const int64_t sublattice[degree][degree], const int64_t sublattice_inv[degree][degree]);

# if IS_TOEPLITZ_USABLE
void reduction_montgomery_toeplitz(int degree, int64_t out[degree], __int128 polynomial[degree], const int64_t sublattice[2*degree - 1], const uint64_t sublattice_inv[2*degree - 1]);

void reduction_montgomery_toeplitz_recursive(int degree, int64_t out[degree], __int128 polynomial[degree], const int64_t sublattice[2 * degree - 1], const uint64_t sublattice_inv[2 * degree - 1]);
# endif

# if IS_DOUBLE_SPARSE
void reduction_montgomery_linear(int extension_degree, int degree, int64_t out[degree], __int128 polynomial[degree]);
# endif

/*=================================================================
                MONTGOMERY REDUCTION FOR CONVERSION 
=================================================================*/
void reduction_montgomery_mpz(int degree, mpz_t out[degree], mpz_t polynomial[degree], const int64_t sublattice[degree][degree], const int64_t sublattice_inv[degree][degree]);

void reduction_montgomery_mpn(int n_limbs, int degree, mp_limb_t out[degree][n_limbs], mp_limb_t pol[degree][n_limbs], const int64_t sublattice[degree][degree], const int64_t sublattice_inv[degree][degree]);
#endif