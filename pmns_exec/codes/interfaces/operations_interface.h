#ifndef OPERATIONS_H
#define OPERATIONS_H

#include <stdint.h>
#include "../params/pmns_params.h"
#include "./interfaces/operations_interface.h"
#include "./interfaces/reductions_interface.h"
#include "../params/pmns_params.h"
#include <string.h>


/*=================================================================
                FUNCTIONS OVER POLYNOMIALS COEFFICIENTS
=================================================================*/
#if IS_BABAI_USABLE
void coeff_shift_i64(int degree, __int128 out[degree], __int128 polynomial[degree], int n_shift);

void coeff_shift_i128(int degree, __int128 out[degree], __int128 polynomial[degree], int n_shift);
#endif

/*=================================================================
        FUNCTIONS FOR PRODUCTS WITH POLYNOMIALS AND INTEGERS 
=================================================================*/
void addmul_pol64_int64(int degree, __int128 out[degree], const int64_t polynomial[degree], int64_t coeff);

void addmul_pol64_int64_mpn(int degree, int n_limbs, mp_limb_t out[degree][n_limbs], const int64_t polynomial[degree], int64_t coeff);

/*=================================================================
        FUNCTIONS FOR PRODUCTS WITH POLYNOMIALS AND POLYNOMIALS 
=================================================================*/
void addmul_polmpn_pol64(int degree, int n_limbs, mp_limb_t out[degree][n_limbs], mp_limb_t pol_mpn[degree][n_limbs], const int64_t pol64[degree]);

void recursive_karatsuba_product(int degree, __int128 out[2 * degree - 1], int64_t pol_A[degree], int64_t pol_B[degree]);

void polynomials_product(int degree, __int128 out[2 * degree - 1], int64_t PolA[degree], int64_t PolB[degree]);

#if IS_TOEPLITZ_USABLE
void addmul_polmpn_Xpow_modE(int degree, int n_limbs, mp_limb_t out[degree][n_limbs], mp_limb_t polmpn[degree][n_limbs], unsigned int pow);
#endif

/*=================================================================
        FUNCTIONS FOR PRODUCTS WITH POLYNOMIALS AND MATRICES 
=================================================================*/
void prod_pol_mat_i64(int degree, __int128 out[degree], __int128 polynomial[degree], const int64_t matrix[degree][degree]);

void prod_pol_mat_i128(int degree, __int128 out[degree], __int128 polynomial[degree], const int64_t matrix[degree][degree]);

# if IS_TOEPLITZ_USABLE
void prod_pol_mat_toeplitz_i64(int degree, __int128 out[degree], const __int128 polynomial[degree], const uint64_t matrix_toeplitz[2*degree - 1]);

void prod_pol_mat_toeplitz_i128(int degree, __int128 out[degree], const __int128 polynomial[degree], const int64_t matrix_toeplitz[2*degree - 1]);

void prod_pol_mat_toeplitz_recursive_i64(int degree, __int128 out[degree], const __int128 vector[degree], const uint64_t toeplitz_matrix[2*degree - 1]);

void prod_pol_mat_toeplitz_recursive_i128(int degree, __int128 out[degree], const __int128 vector[degree], const int64_t toeplitz_matrix[2*degree - 1]);
# endif


#if IS_DOUBLE_SPARSE
void prod_pol_mat_linear_i64(int64_t out[DEGREE], __int128 polynomial[DEGREE]);

void prod_pol_mat_linear_i128(__int128 out[DEGREE], int64_t polynomial[DEGREE]);
#endif

#endif