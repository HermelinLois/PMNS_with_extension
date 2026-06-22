#ifndef OPERATIONS_H
#define OPERATIONS_H

#include <stdint.h>
#include "../params/pmns_params.h"

/* Polynomial product over extension field */
void polynomials_product(__int128 out[DEGREE], int64_t PolA[DEGREE], int64_t PolB[DEGREE]);

void addmul_pol64_int64(__int128 out[DEGREE], const int64_t polynomial[DEGREE], int64_t coeff);

void addmul_pol64_int64_mpn(int n_limbs, mp_limb_t out[DEGREE][n_limbs], const int64_t polynomial[DEGREE], int64_t coeff);

void addmul_polmpn_pol64(int n_limbs, mp_limb_t out[DEGREE][n_limbs], mp_limb_t pol_mpn[DEGREE][n_limbs], const int64_t pol64[DEGREE]);

void prod_pol_mat_i64(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t matrix[DEGREE][DEGREE]);

void prod_pol_mat_i128(__int128 out[DEGREE], int64_t polynomial[DEGREE], const int64_t matrix[DEGREE][DEGREE]);

#if ELEMENTS_EXPRESS_IN_GAMMA_BASIS
void addmul_polmpn_Xpow_modE(int n_limbs, mp_limb_t out[DEGREE][n_limbs], mp_limb_t polmpn[DEGREE][n_limbs], unsigned int pow);
#endif

#if IS_TOEPLITZ_USABLE
void prod_pol_mat_toeplitz_i64(int64_t out[DEGREE], __int128 polynomial[DEGREE], const uint64_t matrix_toeplitz[2*DEGREE - 1]);

void prod_pol_mat_toeplitz_i128(__int128 out[DEGREE], int64_t polynomial[DEGREE], const int64_t matrix_toeplitz[2*DEGREE - 1]);

void toeplitz_recursive_vector_matrix(int n, __int128 *out, const __int128 *vector, const int64_t *toeplitz_matrix, int with_mod_64bits);
#endif

#if IS_BABAI_USABLE
void coeff_shift_i64(int64_t out[DEGREE], __int128 in[DEGREE], int shift);

void coeff_shift_i128(__int128 out[DEGREE], __int128 in[DEGREE], int shift);
#endif

#if IS_DOUBLE_SPARSE
void linear_prod_pol_lattice_inv_i64(int64_t out[DEGREE], __int128 polynomial[DEGREE]);

void linear_prod_pol_lattice_i128(__int128 out[DEGREE], int64_t polynomial[DEGREE]);
#endif
#endif
