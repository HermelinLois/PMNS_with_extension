# include "./interfaces/reductions_interface.h"
# include "./interfaces/operations_interface.h"
# include <stdio.h>
# include <string.h>
# include <gmp.h>

/*=================================================================
                        BABAI REDUCTION FUNCTION 
=================================================================*/
# if IS_BABAI_USABLE
void reduction_babai_int128(int degree, int64_t out[degree], __int128 polynomial[degree], const int64_t sublattice[degree][degree], const int64_t sublattice_inv[degree][degree]) {
    /*Define the Babai internal reduction approach*/
    __int128 S[degree];
    __int128 SL[degree];
    __int128 PH2[degree];

    coeff_shift_i64(degree, PH2, polynomial, H2);
    prod_pol_mat_i128(degree, S, PH2, sublattice_inv);
    coeff_shift_i128(degree, S, S, (H1 - H2));
    prod_pol_mat_i64(degree, SL, S, sublattice);

    for (int i = 0; i < degree; i++)
        out[i] = polynomial[i] - SL[i];
}
# endif


/*=================================================================
                    MONTGOMERY REDUCTION FUNCTION 
=================================================================*/
// define a macro to define a generic structure for the montgomery reduction function
#define MONTGOMERY_REDUCTION_CORE(Q_TYPE, T_TYPE, DEGREE, PRODUCT_FUNC64, PRODUCT_FUNC128, OUT, POLYNOMIAL, LATTICE, LATTICE_INV)   \
    Q_TYPE Q[DEGREE];                                                                                                               \
    T_TYPE T[DEGREE];                                                                                                               \
                                                                                                                                    \
    PRODUCT_FUNC64(DEGREE, Q, POLYNOMIAL, LATTICE_INV);                                                                             \
    PRODUCT_FUNC128(DEGREE, T, Q, LATTICE);                                                                                         \
                                                                                                                                    \
    for (int deg = 0; deg < DEGREE; deg++)                                                                                          \
        OUT[deg] = (T[deg] + POLYNOMIAL[deg]) >> 64;                                                


void reduction_montgomery_int128(int degree, int64_t out[degree], __int128 polynomial[degree], const int64_t sublattice[degree][degree], const int64_t sublattice_inv[degree][degree]) {
    /*Define the Montgomery internal reduction approach with the multiplication of the polynomial and sublattices*/
    MONTGOMERY_REDUCTION_CORE(__int128, __int128, degree, prod_pol_mat_i64, prod_pol_mat_i128, out, polynomial, sublattice, sublattice_inv);
}


# if IS_TOEPLITZ_USABLE
void reduction_montgomery_toeplitz(int degree, int64_t out[degree], __int128 polynomial[degree], const int64_t sublattice[2*degree - 1], const uint64_t sublattice_inv[2*degree - 1]) {
    /*Define the Montgomery Toeplitz internal reduction approach using the Toeplitz structure*/
    MONTGOMERY_REDUCTION_CORE(__int128, __int128, degree, prod_pol_mat_toeplitz_i64, prod_pol_mat_toeplitz_i128, out, polynomial, sublattice, sublattice_inv);
}


void reduction_montgomery_toeplitz_recursive(int degree, int64_t out[degree], __int128 polynomial[degree], const int64_t sublattice[2 * degree - 1], const uint64_t sublattice_inv[2 * degree - 1]){
    /*Define the Montgomery Toeplitz recursive internal reduction approach using the recursive Toeplitz structure*/
    MONTGOMERY_REDUCTION_CORE(__int128, __int128, degree, prod_pol_mat_toeplitz_recursive_i64, prod_pol_mat_toeplitz_recursive_i128, out, polynomial, sublattice, sublattice_inv);
}
# endif


# if IS_DOUBLE_SPARSE
void reduction_montgomery_linear(int extension_degree, int degree, int64_t out[degree], __int128 polynomial[degree]) {
    /*Define the Montgomery internal reduction approach by using the double sparse structure of the sublattice and its inverse.
    This approach is more efficient than the previous ones, as it takes advantage of the sparsity of the sublattice and its inverse to 
    have a linear complexity.*/
    int64_t Q[degree];
    __int128 T[degree];

    prod_pol_mat_linear_i64(extension_degree, degree, Q, polynomial);
    prod_pol_mat_linear_i128(extension_degree, degree, T, Q);

    for (int deg = 0; deg < degree; deg++)
        out[deg] = (T[deg] + polynomial[deg]) >> 64;
}
# endif


/*=================================================================
                MONTGOMERY REDUCTION FOR CONVERSION 
=================================================================*/
static inline void init_phi(mpz_t phi) {
    //Construct the value of phi as 2^PHI_POW
    mpz_init_set_ui(phi, 1);
    mpz_mul_2exp(phi, phi, PHI_POW);
}

void reduction_montgomery_mpz(int degree, mpz_t out[degree], mpz_t polynomial[degree], const int64_t sublattice[degree][degree], const int64_t sublattice_inv[degree][degree]){
    /*Define the Montgomery internal reduction approach with the multiplication of the polynomial and sublattices. It use mpz types to 
    handle multiple precision integers during the computation.*/
    int64_t Q[degree];
    mpz_t T[degree];

    mpz_t tmp, phi;
    mpz_init(tmp);
    init_phi(phi);

    prod_polmpz_mat_i64(degree, Q, polynomial, sublattice_inv, phi);
    prod_pol_mat_mpz(degree, T, Q, sublattice);

    // add the polynomials and divide by phi to get the final result
    for (int deg = 0; deg < degree; deg++) {
        mpz_init_set_ui(tmp, 0);
        mpz_add(tmp, T[deg], polynomial[deg]);
        mpz_fdiv_q_2exp(out[deg], tmp, PHI_POW);
    }

    mpz_clears(tmp, phi, NULL);
    for (int j=0; j<degree; j++) mpz_clear(T[j]);
}


void reduction_montgomery_mpn(int n_limbs, int degree, mp_limb_t out[degree][n_limbs], mp_limb_t pol[degree][n_limbs], const int64_t sublattice[degree][degree], const int64_t sublattice_inv[degree][degree]){
    /*Define the Montgomery internal reduction approach with the multiplication of the polynomial and sublattices. It use mpn types to 
    have efficient computation over multiple precision integers.*/
    int SIZE = n_limbs + 1;
    int64_t Q[degree];
    mp_limb_t T[degree][SIZE];

    prod_polmpn_mat_i64(degree, n_limbs, Q, pol, sublattice_inv);
    prod_pol_mat_mpn(degree, n_limbs, T, Q, sublattice);

    // Do a register shift after the add as PHI_POW is equal to the size of a limb in the current implementation
    mp_limb_t tmp[SIZE];
    for (int deg = 0; deg < degree; deg++) {
        mpn_zero(tmp, SIZE);
        int sP = (pol[deg][n_limbs - 1] >> (GMP_NUMB_BITS - 1)) & 1;

        mp_limb_t coeff_extended[SIZE];
        mpn_copyi(coeff_extended, pol[deg], n_limbs);
        coeff_extended[n_limbs] = sP ? ~0UL : 0UL;
        
        mpn_add_n(tmp, coeff_extended, T[deg], SIZE);
        mpn_copyi(out[deg], tmp + 1, n_limbs);
    }
}