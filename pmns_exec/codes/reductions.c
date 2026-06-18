# include "./interfaces/reductions_interface.h"
# include "./interfaces/operations_interface.h"
# include <stdio.h>
# include <string.h>
# include <gmp.h>


# ifndef IS_SPARSE
void reduction_babai_int128(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]) {
    __int128 S[DEGREE];
    int64_t SL[DEGREE];
    int64_t PH2[DEGREE];

    coeff_shift_i64(PH2, polynomial, H2);
    prod_pol_mat_i128(S, PH2, sublattice_inv);
    coeff_shift_i128(S, S, (H1 - H2));
    prod_pol_mat_i64(SL, S, sublattice);

    for (int i = 0; i < DEGREE; i++)
        out[i] = polynomial[i] - SL[i];
}

# else
void reduction_montgomery_linear(int64_t out[DEGREE], __int128 polynomial[DEGREE]) {
    int64_t Q[DEGREE] = {0};
    __int128 T[DEGREE] = {0};

    linear_prod_pol_lattice_inv_i64(Q, polynomial);
    linear_prod_pol_lattice_i128(T, Q);

    for (int deg = 0; deg < DEGREE; deg++)
        out[deg] = (T[deg] + polynomial[deg]) >> 64;
}
# endif

void reduction_montgomery_int128(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]) {
    int64_t Q[DEGREE] = {0};
    __int128 T[DEGREE] = {0};

    prod_pol_mat_i64(Q, polynomial, sublattice_inv);
    prod_pol_mat_i128(T, Q, sublattice);

    for (int deg = 0; deg < DEGREE; deg++)
        out[deg] = (T[deg] + polynomial[deg]) >> 64;
}


void reduction_montgomery_toeplitz(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[2*DEGREE - 1], const uint64_t sublattice_inv[2*DEGREE - 1]) {
    int64_t Q[DEGREE] = {0};
    __int128 T[DEGREE] = {0};

    prod_pol_mat_toeplitz_i64(Q, polynomial, sublattice_inv);
    prod_pol_mat_toeplitz_i128(T, Q, sublattice);

    for (int deg = 0; deg < DEGREE; deg++)
        out[deg] = (T[deg] + polynomial[deg]) >> 64;
}




void reduction_montgomery_mpn(int n_limbs, mp_limb_t out[DEGREE][n_limbs], mp_limb_t pol[DEGREE][n_limbs], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]){
    int SIZE = n_limbs + 1;

    int64_t Q[DEGREE];
    mp_limb_t T[DEGREE][SIZE];
    mp_limb_t tmp[SIZE];

    mp_limb_t acc[SIZE];
    mp_limb_t prod[SIZE];

    for (int j = 0; j < DEGREE; j++) {
        mpn_zero(acc, SIZE);

        for (int i = 0; i < DEGREE; i++) {
            int64_t mat_coeff = sublattice_inv[i][j];
            if (mat_coeff == 0) continue;

            int mat_sign = (mat_coeff < 0);
            uint64_t abs_c = mat_sign ? (uint64_t)(-(uint64_t)mat_coeff) : (uint64_t)mat_coeff;
            int pol_sign = (pol[i][n_limbs - 1] >> (GMP_NUMB_BITS - 1)) & 1;
            
            mp_limb_t tmp_pol[n_limbs];
            if (pol_sign){
                mpn_com(tmp_pol, pol[i], n_limbs);
                mpn_add_1(tmp_pol, tmp_pol, n_limbs, 1);
            } else {
                mpn_copyi(tmp_pol, pol[i], n_limbs);
            }

            prod[n_limbs] = mpn_mul_1(prod, tmp_pol, n_limbs, abs_c);
            (mat_sign ^ pol_sign ? mpn_sub_n : mpn_add_n)(acc, acc, prod, SIZE);
        }
        Q[j] = acc[0];
    }


    for (int j = 0; j < DEGREE; j++) {
        mpn_zero(T[j], SIZE);

        for (int i = 0; i < DEGREE; i++) {
            int64_t mat_coeff = sublattice[i][j];
            int64_t Q_coeff   = Q[i];
            if (mat_coeff == 0 || Q_coeff == 0) continue;
            
            int mat_sign = (mat_coeff < 0);
            int Q_sign   = (Q_coeff < 0);
            uint64_t mat_abs = mat_sign ? -(uint64_t)mat_coeff : (uint64_t)mat_coeff;
            uint64_t Q_abs   = Q_sign   ? (uint64_t)(-Q_coeff) : (uint64_t)Q_coeff;

            mp_limb_t mat_limb = (mp_limb_t)mat_abs;
            mp_limb_t Q_limb   = (mp_limb_t)Q_abs;

            mp_limb_t prod_limbs[SIZE];
            mpn_zero(prod_limbs, SIZE);

            prod_limbs[1] = mpn_mul_1(prod_limbs, &Q_limb, 1, mat_limb);

            int sign_prod = mat_sign ^ Q_sign;
            (sign_prod ? mpn_sub_n : mpn_add_n)(T[j], T[j], prod_limbs, SIZE);
        }
    }

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


void reduction_montgomery_mpz(mpz_t out[DEGREE], mpz_t polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]){
    int64_t Q[DEGREE];
    mpz_t T[DEGREE], tmp, acc, phi;
    mpz_inits(tmp, acc, NULL);
    mpz_init_set_ui(phi, 1);
    mpz_mul_2exp(phi, phi, PHI_POW);  

    for (int j = 0; j < DEGREE; j++) {
        mpz_set_ui(acc, 0);

        for (int i = 0; i < DEGREE; i++) {
            mpz_mul_ui(tmp, polynomial[i], sublattice_inv[i][j]);
            mpz_add(acc, acc, tmp);
        }
        mpz_mod(acc, acc, phi);
        Q[j] = mpz_get_ui(acc);
    }


    for (int j = 0; j < DEGREE; j++) {
        mpz_init_set_ui(T[j], 0);

        for (int i = 0; i < DEGREE; i++) {
            int64_t mat_coeff = sublattice[i][j];
            int64_t Q_coeff   = Q[i];
            if (mat_coeff == 0 || Q_coeff == 0) continue;
            
            __int128 prod = (__int128)mat_coeff * (__int128)Q_coeff;
            unsigned __int128 abs_prod = prod < 0 ? -prod : prod;

            mpz_import(tmp, 2, -1, sizeof(uint64_t), 0, 0, &abs_prod);
            if (prod < 0) mpz_neg(tmp, tmp);
            mpz_add(T[j], T[j], tmp);
        }
    }

    for (int deg = 0; deg < DEGREE; deg++) {
        mpz_init_set_ui(tmp, 0);
        mpz_add(tmp, T[deg], polynomial[deg]);
        mpz_fdiv_q_2exp(out[deg], tmp, PHI_POW);
    }

    mpz_clears(tmp, phi, NULL);
    for (int j=0; j<DEGREE; j++) mpz_clear(T[j]);
}




























void small_toeplitz_vector_matrix(int n,  const __int128 *vector, const int64_t *toeplitz_matrix, __int128 *out, int with_mod_64bits){
    for (int i = 0; i < n; i++) {
        __int128 acc = 0;

        for (int j = 0; j < n; j++) {
            acc += vector[j] * (__int128)toeplitz_matrix[n - 1 + (i - j)];
        }

        out[i] = with_mod_64bits ? (int64_t)acc : acc;
    }
}


void toeplitz_recursive_vector_matrix(int n, __int128 *out, const __int128 *vector, const int64_t *toeplitz_matrix, int with_mod_64bits){
    if ((n % 2 != 0 && n % 3 != 0) || n <=3) {
        small_toeplitz_vector_matrix(n, vector, toeplitz_matrix, out, with_mod_64bits);
        return;
    }

    else if (n % 2 == 0) {
        int m = n / 2;

        // Repartition of the input vector into 3 parts
        const __int128 *v0 = vector;
        const __int128 *v1 = vector + m;
        __int128 v_sum_01[m];

        for (int i = 0; i < m; i++)
            v_sum_01[i] = v0[i] + v1[i];

        const int64_t *M0 = toeplitz_matrix;
        const int64_t *M1 = M0 + m;
        const int64_t *M2 = M1 + m;

        // Create new submatrices for the recursive calls
        int64_t T_sub_M12[2 * m - 1];
        int64_t T_sub_M10[2 * m - 1];

        for (int i = 0; i < 2 * m - 1; i++) {
            T_sub_M12[i] = M1[i] - M2[i];
            T_sub_M10[i] = M1[i] - M0[i];
        }

        // Recursive calls to compute the products with the submatrices
        __int128 p0[m], p1[m], p2[m];
        memset(p0, 0, sizeof(p0));
        memset(p1, 0, sizeof(p1));
        memset(p2, 0, sizeof(p2));

        toeplitz_recursive_vector_matrix(m, p0, v_sum_01, M1, with_mod_64bits);
        toeplitz_recursive_vector_matrix(m, p1, v0, T_sub_M12, with_mod_64bits);
        toeplitz_recursive_vector_matrix(m, p2, v1, T_sub_M10, with_mod_64bits);

        // Combine the results of the recursive calls to form the final output
        for (int i = 0; i < m; i++) {
            out[i] = (with_mod_64bits) ? (__int128)(int64_t)(p0[i] - p2[i]) : p0[i] - p2[i];
            out[i + m] = (with_mod_64bits) ? (__int128)(int64_t)(p0[i] - p1[i]) : p0[i] - p1[i];
        }
    }

    else if (n % 3 == 0) {
        int m = n / 3;

        // Repartition of the input vector into 3 parts
        const __int128 *v0 = vector;
        const __int128 *v1 = vector + m;
        const __int128 *v2 = vector + 2 * m;
        __int128 v_sub_v02[m];
        __int128 v_sub_v12[m];
        __int128 v_sub_v01[m];

        for (int i = 0; i < m; i++) {
            v_sub_v02[i] = v0[i] - v2[i];
            v_sub_v12[i] = v1[i] - v2[i];
            v_sub_v01[i] = v0[i] - v1[i];
        }

        const int64_t *M0 = toeplitz_matrix;
        const int64_t *M1 = M0 + m;
        const int64_t *M2 = M1 + m;
        const int64_t *M3 = M2 + m;
        const int64_t *M4 = M3 + m;

        int64_t T_sum_M012[2 * m - 1];
        int64_t T_sum_M123[2 * m - 1];
        int64_t T_sum_M234[2 * m - 1];

        for (int i = 0; i < 2 * m - 1; i++) {
            T_sum_M012[i] = M0[i] + M1[i] + M2[i];
            T_sum_M123[i] = M1[i] + M2[i] + M3[i];
            T_sum_M234[i] = M2[i] + M3[i] + M4[i];
        }

        __int128 p0[m], p1[m], p2[m], p3[m], p4[m], p5[m];
        memset(p0, 0, sizeof(p0));
        memset(p1, 0, sizeof(p1));
        memset(p2, 0, sizeof(p2));
        memset(p3, 0, sizeof(p3));
        memset(p4, 0, sizeof(p4));
        memset(p5, 0, sizeof(p5));
        

        toeplitz_recursive_vector_matrix(m, p0, v2, T_sum_M012, with_mod_64bits);
        toeplitz_recursive_vector_matrix(m, p1, v1, T_sum_M123, with_mod_64bits);
        toeplitz_recursive_vector_matrix(m, p2, v0, T_sum_M234, with_mod_64bits);

        toeplitz_recursive_vector_matrix(m, p3, v_sub_v02, M2, with_mod_64bits);
        toeplitz_recursive_vector_matrix(m, p4, v_sub_v12, M1, with_mod_64bits);
        toeplitz_recursive_vector_matrix(m, p5, v_sub_v01, M3, with_mod_64bits);

        for (int i = 0; i < m; i++) {
            out[i] = with_mod_64bits ? (__int128)(int64_t)(p0[i] + p3[i] + p4[i]) : p0[i] + p3[i] + p4[i];
            out[i + m] = with_mod_64bits ? (__int128)(int64_t)(p1[i] + p5[i] - p4[i]) : p1[i] + p5[i] - p4[i];
            out[i + 2 * m] = with_mod_64bits ? (__int128)(int64_t)(p2[i] - p3[i] - p5[i]) : p2[i] - p3[i] - p5[i];
        }
    }
}


void reduction_montgomery_toeplitz_recursive(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[2 * DEGREE - 1],  const uint64_t sublattice_inv[2 * DEGREE - 1]){
    __int128 Q[DEGREE] = {0};
    __int128 T[DEGREE] = {0};

    toeplitz_recursive_vector_matrix(DEGREE, Q, polynomial, (const int64_t *)sublattice_inv, 1);
    toeplitz_recursive_vector_matrix(DEGREE, T, Q, sublattice, 0);

    for (int i = 0; i < DEGREE; i++)
        out[i] = (int64_t)((T[i] + polynomial[i]) >> 64);
}