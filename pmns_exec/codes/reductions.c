# include "./interfaces/reductions_interface.h"
# include "./interfaces/operations_interface.h"
# include <stdio.h>
# include <string.h>
# include <gmp.h>

/*=================================================================
                        BABAI REDUCTION FUNCTION 
=================================================================*/
# if IS_BABAI_USABLE
void reduction_babai_int128(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]) {
    __int128 S[DEGREE];
    __int128 SL[DEGREE];
    __int128 PH2[DEGREE];

    coeff_shift_i64(DEGREE, PH2, polynomial, H2);
    prod_pol_mat_i128(DEGREE, S, PH2, sublattice_inv);
    coeff_shift_i128(DEGREE, S, S, (H1 - H2));
    prod_pol_mat_i64(DEGREE, SL, S, sublattice);

    for (int i = 0; i < DEGREE; i++)
        out[i] = polynomial[i] - SL[i];
}
# endif


/*=================================================================
                    MONTGOMERY REDUCTION FUNCTION 
=================================================================*/
#define MONTGOMERY_REDUCTION_CORE(Q_TYPE, T_TYPE, DEGREE, PRODUCT_FUNC64, PRODUCT_FUNC128, OUT, POLYNOMIAL, LATTICE, LATTICE_INV)   \
    Q_TYPE Q[DEGREE];                                                                                                               \
    T_TYPE T[DEGREE];                                                                                                               \
                                                                                                                                    \
    PRODUCT_FUNC64(DEGREE, Q, POLYNOMIAL, LATTICE_INV);                                                                             \
    PRODUCT_FUNC128(DEGREE, T, Q, LATTICE);                                                                                         \
                                                                                                                                    \
    for (int deg = 0; deg < DEGREE; deg++)                                                                                          \
        OUT[deg] = (T[deg] + POLYNOMIAL[deg]) >> 64;                                                


void reduction_montgomery_int128(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]) {
    MONTGOMERY_REDUCTION_CORE(__int128, __int128, DEGREE, prod_pol_mat_i64, prod_pol_mat_i128, out, polynomial, sublattice, sublattice_inv);
}


# if IS_TOEPLITZ_USABLE
void reduction_montgomery_toeplitz(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[2*DEGREE - 1], const uint64_t sublattice_inv[2*DEGREE - 1]) {
    MONTGOMERY_REDUCTION_CORE(__int128, __int128, DEGREE, prod_pol_mat_toeplitz_i64, prod_pol_mat_toeplitz_i128, out, polynomial, sublattice, sublattice_inv);
}


void reduction_montgomery_toeplitz_recursive(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[2 * DEGREE - 1], const uint64_t sublattice_inv[2 * DEGREE - 1]){
    MONTGOMERY_REDUCTION_CORE(__int128, __int128, DEGREE, prod_pol_mat_toeplitz_recursive_i64, prod_pol_mat_toeplitz_recursive_i128, out, polynomial, sublattice, sublattice_inv);
}
# endif


# if IS_DOUBLE_SPARSE
void reduction_montgomery_linear(int64_t out[DEGREE], __int128 polynomial[DEGREE]) {
    int64_t Q[DEGREE] = {0};
    __int128 T[DEGREE] = {0};

    prod_pol_mat_linear_i64(Q, polynomial);
    prod_pol_mat_linear_i128(T, Q);

    for (int deg = 0; deg < DEGREE; deg++)
        out[deg] = (T[deg] + polynomial[deg]) >> 64;
}
# endif


/*=================================================================
                MONTGOMERY REDUCTION FOR CONVERSION 
=================================================================*/
static inline void init_phi(mpz_t phi) {
    mpz_init_set_ui(phi, 1);
    mpz_mul_2exp(phi, phi, PHI_POW);
}

void reduction_montgomery_mpz(mpz_t out[DEGREE], mpz_t polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]){
    int64_t Q[DEGREE];
    mpz_t T[DEGREE];

    mpz_t tmp, phi;
    mpz_init(tmp);
    init_phi(phi);

    prod_polmpz_mat_i64(DEGREE, Q, polynomial, sublattice_inv, phi);
    prod_pol_mat_mpz(DEGREE, T, Q, sublattice);

    for (int deg = 0; deg < DEGREE; deg++) {
        mpz_init_set_ui(tmp, 0);
        mpz_add(tmp, T[deg], polynomial[deg]);
        mpz_fdiv_q_2exp(out[deg], tmp, PHI_POW);
    }

    mpz_clears(tmp, phi, NULL);
    for (int j=0; j<DEGREE; j++) mpz_clear(T[j]);
}


void reduction_montgomery_mpn(int n_limbs, mp_limb_t out[DEGREE][n_limbs], mp_limb_t pol[DEGREE][n_limbs], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]){
    int SIZE = n_limbs + 1;
    int64_t Q[DEGREE];
    mp_limb_t T[DEGREE][SIZE];
    
    prod_polmpn_mat_i64(DEGREE, n_limbs, Q, pol, sublattice_inv);
    prod_pol_mat_mpn(DEGREE, n_limbs, T, Q, sublattice);

    mp_limb_t tmp[SIZE];
    for (int deg = 0; deg < DEGREE; deg++) {
        mpn_zero(tmp, SIZE);
        int sP = (pol[deg][n_limbs - 1] >> (GMP_NUMB_BITS - 1)) & 1;

        mp_limb_t coeff_extended[SIZE];
        mpn_copyi(coeff_extended, pol[deg], n_limbs);
        coeff_extended[n_limbs] = sP ? ~0UL : 0UL;
        
        mpn_add_n(tmp, coeff_extended, T[deg], SIZE);
        mpn_copyi(out[deg], tmp + 1, n_limbs);
    }
}