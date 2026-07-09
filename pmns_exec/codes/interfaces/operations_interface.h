#ifndef OPERATIONS_H
#define OPERATIONS_H

#include <stdint.h>
#include "../params/pmns_params.h"
#include "./interfaces/operations_interface.h"
#include "./interfaces/reductions_interface.h"
#include "../params/pmns_params.h"
#include <string.h>

#include "./interfaces/operations_interface.h"
#include "./interfaces/reductions_interface.h"
#include "../params/pmns_params.h"
#include <string.h>

#define TOEPLITZ_THRESHOLD 50
#define KARATSUBA_THRESHOLD 30

/*=================================================================
                FUNCTIONS OVER POLYNOMIALS COEFFICIENTS
=================================================================*/
# if BABAI_IS_USABLE
void coeff_shift_i64(__int128 out[DEGREE], __int128 polynomial[DEGREE], int n_shift);

void coeff_shift_i128(__int128 out[DEGREE], __int128 polynomial[DEGREE], int n_shift);
#endif

/*=================================================================
        FUNCTIONS FOR PRODUCTS WITH POLYNOMIALS AND INTEGERS 
=================================================================*/
void addmul_pol64_int64(__int128 out[DEGREE], const int64_t polynomial[DEGREE], int64_t coeff);

void addmul_pol64_int64_mpn(int n_limbs, mp_limb_t out[DEGREE][n_limbs], const int64_t polynomial[DEGREE], int64_t coeff);

/*=================================================================
        FUNCTIONS FOR PRODUCTS WITH POLYNOMIALS AND POLYNOMIALS 
=================================================================*/

void addmul_polmpn_pol64(int n_limbs, mp_limb_t out[DEGREE][n_limbs], mp_limb_t pol_mpn[DEGREE][n_limbs], const int64_t pol64[DEGREE]);

void recursive_karatsuba_product(__int128 out[2 * DEGREE - 1], int64_t pol_A[DEGREE], int64_t pol_B[DEGREE]);

void polynomials_product(__int128 out[DEGREE], int64_t PolA[DEGREE], int64_t PolB[DEGREE]);

# if ELEMENTS_ARE_IN_GAMMA_BASIS
void addmul_polmpn_Xpow_modE(int n_limbs, mp_limb_t out[DEGREE][n_limbs], mp_limb_t polmpn[DEGREE][n_limbs], unsigned int pow);
#endif

/*=================================================================
        FUNCTIONS FOR PRODUCTS WITH POLYNOMIALS AND MATRICES 
=================================================================*/
void prod_polmpn_mat_i64(int n_limbs, int64_t out[DEGREE], mp_limb_t pol[DEGREE][n_limbs], const int64_t sublattice[DEGREE][DEGREE]);

void prod_pol_mat_mpn(int n_limbs, mp_limb_t out[DEGREE][n_limbs+1], int64_t pol[DEGREE], const int64_t sublattice[DEGREE][DEGREE]);

void prod_pol_mat_mpz(mpz_t out[DEGREE], int64_t pol[DEGREE], const int64_t sublattice[DEGREE][DEGREE]);

void prod_polmpz_mat_i64(int64_t out[DEGREE], mpz_t pol[DEGREE], const int64_t sublattice[DEGREE][DEGREE], mpz_t phi);

void prod_pol_mat_i64(__int128 out[DEGREE], __int128 polynomial[DEGREE], const int64_t matrix[DEGREE][DEGREE]);

void prod_pol_mat_i128(__int128 out[DEGREE], __int128 polynomial[DEGREE], const int64_t matrix[DEGREE][DEGREE]);

# if LATTICE_IS_DOUBLE_SPARSE
void prod_pol_mat_toeplitz_i64(__int128 out[DEGREE], const __int128 polynomial[DEGREE], const uint64_t matrix_toeplitz[2*DEGREE - 1]);

void prod_pol_mat_toeplitz_i128(__int128 out[DEGREE], const __int128 polynomial[DEGREE], const int64_t matrix_toeplitz[2*DEGREE - 1]);

void prod_pol_mat_toeplitz_recursive_i64(__int128 out[DEGREE], const __int128 vector[DEGREE], const uint64_t toeplitz_matrix[2*DEGREE - 1]);

void prod_pol_mat_toeplitz_recursive_i128(__int128 out[DEGREE], const __int128 vector[DEGREE], const int64_t toeplitz_matrix[2*DEGREE - 1]);
# endif

# if IS_DOUBLE_SPARSE
void prod_pol_mat_linear_i64(int64_t out[DEGREE], __int128 polynomial[DEGREE]);

void prod_pol_mat_linear_i128(__int128 out[DEGREE], int64_t polynomial[DEGREE]);
#endif
#endif