#include "./interfaces/operations_interface.h"
#include "./interfaces/reductions_interface.h"
#include "../params/pmns_params.h"
#include <string.h>

#define TOEPLITZ_THRESHOLD 50
#define KARATSUBA_THRESHOLD 30

// compute the product of an int64_t polynomial with an int64_t coefficient and add the result to an int128 polynomial
void addmul_pol64_int64(__int128 out[DEGREE], const int64_t polynomial[DEGREE], int64_t coeff){
    if (coeff == 0) return;

    for (int deg=0; deg<DEGREE; deg++)
        out[deg] += (__int128)polynomial[deg] * (__int128)coeff;
}

void addmul_pol64_int64_mpn(int n_limbs, mp_limb_t out[DEGREE][n_limbs], const int64_t polynomial[DEGREE], int64_t coeff){
    if (coeff == 0) return;

    for (int deg=0; deg<DEGREE; deg++){
        __int128 product = (__int128)polynomial[deg] * (__int128)coeff;

        int prod_sign = (product < 0);
        unsigned __int128 abs_product = prod_sign ? (-product) : product;

        mp_limb_t mpn_product[n_limbs];
        mpn_zero(mpn_product, n_limbs);
        mpn_product[0] = (mp_limb_t)abs_product;
        mpn_product[1] = (mp_limb_t)(abs_product >> GMP_NUMB_BITS);
        
        (prod_sign ? mpn_sub_n : mpn_add_n)(out[deg], out[deg], mpn_product, n_limbs);
    }
}

// compute the product of an int128 polynomial with an int64_t coefficient and add the result to an mp_limb_t polynomial
// use complement to 2 representation for negative values
void addmul_polmpn_pol64(int n_limbs, mp_limb_t out[DEGREE][n_limbs], mp_limb_t pol_mpn[DEGREE][n_limbs], const int64_t pol64[DEGREE]) {
    int PROD_SIZE = 2 * DEGREE - 1;
    mp_limb_t tmp_pol[PROD_SIZE][n_limbs];

    for (int i=0; i<PROD_SIZE; i++) mpn_zero(tmp_pol[i], n_limbs);

    // polynomial product
    for (int i = 0; i < DEGREE; i++) {
        for (int j = 0; j < DEGREE; j++) {
            if (pol64[j] == 0) continue;

            int sign = (pol64[j] < 0);
            uint64_t abs_c = sign ? (uint64_t)(-(uint64_t)pol64[j]) : (uint64_t)pol64[j];
            (sign ? mpn_submul_1 : mpn_addmul_1)(tmp_pol[i+j], pol_mpn[i], n_limbs, (mp_limb_t)abs_c);
        }
    }

    for (int deg=DEGREE; deg<PROD_SIZE; deg++){
        int affect_deg = deg - DEGREE;

        for (int i=0; i<DEGREE; i++) {
            int64_t mat_coeff = EXT_MAT[affect_deg][i];
            if (mat_coeff == 0) continue;


            int sign = (mat_coeff < 0);
            uint64_t abs_coeff = sign ? (uint64_t)(-(uint64_t)mat_coeff) : (uint64_t)mat_coeff;

            (sign ? mpn_submul_1 : mpn_addmul_1)(out[i], tmp_pol[deg], n_limbs, (mp_limb_t)abs_coeff);
        }
    }

    for (int i = 0; i < DEGREE; i++)
        mpn_add_n(out[i], out[i], tmp_pol[i], n_limbs);
}

// /!\ this function only works if E = X^n - lambda
// multiply polynomial by X^pow and apply external reduction if needed, then add the result to out polynomial    
void addmul_polmpn_Xpow_modE(int n_limbs, mp_limb_t out[DEGREE][n_limbs], mp_limb_t polmpn[DEGREE][n_limbs], unsigned int pow) {
    // extract lambda from the external reduction matrix for the given power
    const int64_t lambda = EXT_MAT[0][0];
    const int lambda_sign = (lambda < 0);
    const uint64_t abs_lambda = lambda_sign ? (uint64_t)(-lambda) : (uint64_t)lambda;

    // multiply by X^pow consist on a register shift, then multiply coefficients by lambda to apply external reduction if needed
    // the result is added to the output polynomial. We use complement to 2 representation for negative values and directly apply
    // the external reduction during the shift if needed to avoid extra computations.
    for (int deg=0; deg<DEGREE; deg++){
        // compute affectation degree
        int curr_deg = deg + pow;
        int prod_sign = 0;

        // correct affectation degree and sign of the product if external reduction is needed
        if (curr_deg >= DEGREE){
            mpn_mul_1(polmpn[deg], polmpn[deg], n_limbs, abs_lambda);

            curr_deg -= DEGREE;
            prod_sign ^= lambda_sign;
        }
        // add the product to the output polynomial with the correct sign
        (prod_sign ? mpn_sub_n : mpn_add_n)(out[curr_deg], out[curr_deg], polmpn[deg], n_limbs);
    }
}


void prod_pol_mat_i64(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t matrix[DEGREE][DEGREE]) {
    for (int i = 0; i < DEGREE; i++) {
        int64_t acc = 0;
        for (int j = 0; j < DEGREE; j++) {
            acc += polynomial[j] * matrix[j][i];
        }
        out[i] = (uint64_t)acc;
    }
}

void prod_pol_mat_i128(__int128 out[DEGREE], int64_t polynomial[DEGREE], const int64_t matrix[DEGREE][DEGREE]) {
    for (int i = 0; i < DEGREE; i++) {
        __int128 acc = 0;
        for (int j = 0; j < DEGREE; j++) {
            acc += (__int128)polynomial[j] * (__int128)matrix[j][i];
        }
        out[i] = acc;
    }
}

void recursive_karatsuba_product(int n, __int128 *out, int64_t *pol_A, int64_t *pol_B) {
    // Compute the product of two polynomials using schoolbook multiplication
    if (n < KARATSUBA_THRESHOLD || ((n & 1) == 1)) {
        for (int i = 0; i < 2 * n - 1; i++) out[i] = 0;

        for (int i = 0; i < n; i++) {
            __int128 ai = (__int128)pol_A[i];
            for (int j = 0; j < n; j++) {
                out[i + j] += ai * pol_B[j];
            }
        }
        return;
    }

    int m = n / 2;

    // Divide the polynomials into two halves
    int64_t *A0 = pol_A;
    int64_t *A1 = pol_A + m;
    int64_t *B0 = pol_B;
    int64_t *B1 = pol_B + m;

    // Compute the sums of the halves to use in the Karatsuba formula
    int64_t A_sum_01[m], B_sum_01[m];
    for (int i = 0; i < m; i++) {
        A_sum_01[i] = A0[i] + A1[i];
        B_sum_01[i] = B0[i] + B1[i];
    }

    // Compute the three products needed for Karatsuba's method
    __int128 Y0[2 * m - 1], Y2[2 * m - 1], Y1[2 * m - 1];
    recursive_karatsuba_product(m, Y0, A0, B0);
    recursive_karatsuba_product(m, Y2, A1, B1);
    recursive_karatsuba_product(m, Y1, A_sum_01, B_sum_01);

    // Initialize the output polynomial and recombine the three products to get the final result
    for (int i = 0; i < 2 * n - 1; i++) out[i] = 0;

    for (int i = 0; i < 2 * m - 1; i++) {
        out[i] += Y0[i];
        out[i + m] += Y1[i] - Y0[i] - Y2[i];
        out[i + 2 * m] += Y2[i];
    }
}



void polynomials_product(__int128 out[DEGREE], int64_t PolA[DEGREE], int64_t PolB[DEGREE]){
    int PROD_SIZE = 2*DEGREE - 1;
    __int128 P[PROD_SIZE];

    recursive_karatsuba_product(DEGREE, P, PolA, PolB);

    // Directly copy the first DEGREE coefficients to the output polynomial
    for (int i = 0; i < DEGREE; i++) out[i] = P[i];

    #if IS_TOEPLITZ_USABLE
    // With this condition, we know our external reduction matrix as a form X^n - lambda, 
    // so we can optimize the external reduction by using the first column of the matrix and 
    // avoid unnecessary multiplications by directly extracting diagonal elements for external reduction.
    int64_t lambda = EXT_MAT[0][0];
    for (int j = DEGREE; j < PROD_SIZE; j++) {
        __int128 p_val = P[j];
        if (p_val == 0) continue;

        out[j-DEGREE] += P[j] * (__int128)lambda;
    }
    #else
    // Apply generic external reduction using the external reduction matrix for the remaining coefficients
    for (int j = 0; j < DEGREE - 1; j++) {
        __int128 p_val = P[DEGREE + j];
        if (p_val == 0) continue;

        for (int i = 0; i < DEGREE; i++) {
            out[i] += p_val * (__int128)EXT_MAT[j][i];
        }
    }
    #endif
}





# if IS_TOEPLITZ_USABLE
void prod_pol_mat_toeplitz_i64(int64_t out[DEGREE], __int128 polynomial[DEGREE], const uint64_t matrix_toeplitz[2*DEGREE - 1]) {
    for (int i = 0; i < DEGREE; i++) {
        __int128 acc = 0;
        for (int j = 0; j < DEGREE; j++) {
            acc += polynomial[j] * (__int128)matrix_toeplitz[DEGREE - 1 - j + i];
        }
        out[i] = (int64_t)acc;
    }
}

void prod_pol_mat_toeplitz_i128(__int128 out[DEGREE], int64_t polynomial[DEGREE], const int64_t matrix_toeplitz[2*DEGREE - 1]) {
    for (int i = 0; i < DEGREE; i++) {
        __int128 acc = 0;
        for (int j = 0; j < DEGREE; j++) {
            acc += (__int128)polynomial[j] * (__int128)matrix_toeplitz[DEGREE - 1 - j + i];
        }
        out[i] = acc;
    }
}


void small_toeplitz_vector_matrix_i64(int n, const __int128 *vector, const uint64_t *toeplitz_matrix, __int128 *out){
    for (int i = 0; i < n; i++) {
        int64_t acc = 0;
        for (int j = 0; j < n; j++) {
            acc += vector[j] * (__int128)toeplitz_matrix[n - 1 + i - j];
        }
        out[i] = acc;
    }
}

void small_toeplitz_vector_matrix_i128(int n, const __int128 *vector, const int64_t *toeplitz_matrix, __int128 *out){
    for (int i = 0; i < n; i++) {
        __int128 acc = 0;

        for (int j = 0; j < n; j++)
            acc += (__int128)vector[j] * (__int128)toeplitz_matrix[n - 1 + i - j];

        out[i] = acc;
    }
}


void prod_pol_mat_toeplitz_recursive_i64(int n, __int128 *out, const __int128 *vector, const uint64_t *toeplitz_matrix){
    if (n < TOEPLITZ_THRESHOLD || ((n%3 !=0) && ((n & 1) == 1))) {
        small_toeplitz_vector_matrix_i64(n, vector, toeplitz_matrix, out);
        return;
    }

    else if (n % 3 == 0) {
        int m = n / 3;

        // Decompose the input vector into three parts
        const __int128 *v0 = vector;
        const __int128 *v1 = vector + m;
        const __int128 *v2 = vector + 2 * m;

        // Additional vectors for the differences of the three parts of the input vector
        __int128 v_sub_02[m], v_sub_12[m], v_sub_01[m];
        for (int i = 0; i < m; i++) {
            v_sub_02[i] = v0[i] - v2[i];
            v_sub_12[i] = v1[i] - v2[i];
            v_sub_01[i] = v0[i] - v1[i];
        }

        // Decompose the Toeplitz matrix into five parts corresponding to the three parts of the input vector
        const uint64_t *M0 = toeplitz_matrix;
        const uint64_t *M1 = toeplitz_matrix + m;
        const uint64_t *M2 = toeplitz_matrix + 2 * m;
        const uint64_t *M3 = toeplitz_matrix + 3 * m;
        const uint64_t *M4 = toeplitz_matrix + 4 * m;

        // Preparation of matrix sums (size of a sub-matrix : 2m - 1)
        uint64_t M_sum_012[2 * m - 1], M_sum_123[2 * m - 1], M_sum_234[2 * m - 1];
        for (int i = 0; i < 2 * m - 1; i++) {
            M_sum_012[i] = M0[i] + M1[i] + M2[i];
            M_sum_123[i] = M1[i] + M2[i] + M3[i];
            M_sum_234[i] = M2[i] + M3[i] + M4[i];
        }

        // Compute the products of the sub-vectors with the corresponding sub-matrices using recursive calls
        __int128 p0[m], p1[m], p2[m], p3[m], p4[m], p5[m];
        prod_pol_mat_toeplitz_recursive_i64(m, p0, v2, M_sum_012);
        prod_pol_mat_toeplitz_recursive_i64(m, p1, v1, M_sum_123);
        prod_pol_mat_toeplitz_recursive_i64(m, p2, v0, M_sum_234);
        prod_pol_mat_toeplitz_recursive_i64(m, p3, v_sub_02, M2);
        prod_pol_mat_toeplitz_recursive_i64(m, p4, v_sub_12, M1);
        prod_pol_mat_toeplitz_recursive_i64(m, p5, v_sub_01, M3);

        // Recombine the results of the sub-products to obtain the final result
        for (int i = 0; i < m; i++) {
            out[i] = (int64_t)(p0[i] + p3[i] + p4[i]);
            out[i + m] = (int64_t)(p1[i] + p5[i] - p4[i]);
            out[i + 2 * m] = (int64_t)(p2[i] - p3[i] - p5[i]);
        }
    }
    
    else if ((n & 1) == 0){
        int m = n / 2;
        const __int128 *v0 = vector;
        const __int128 *v1 = vector + m;
            
        __int128 v_sub_01[m];
        for (int i = 0; i < m; i++)
            v_sub_01[i] = v0[i] - v1[i];

        const uint64_t *M0 = toeplitz_matrix;
        const uint64_t *M1 = toeplitz_matrix + m;
        const uint64_t *M2 = toeplitz_matrix + 2 * m;

        uint64_t T_sum_M01[2 * m - 1];
        uint64_t T_sum_M12[2 * m - 1];

        for (int i = 0; i < 2 * m - 1; i++) {
            T_sum_M01[i] = M0[i] + M1[i];
            T_sum_M12[i] = M1[i] + M2[i];
        }

        __int128 p0[m], p1[m], p2[m];

        prod_pol_mat_toeplitz_recursive_i64(m, p0, v1, T_sum_M01);
        prod_pol_mat_toeplitz_recursive_i64(m, p1, v0, T_sum_M12);
        prod_pol_mat_toeplitz_recursive_i64(m, p2, v_sub_01, M1);

        for (int i = 0; i < m; i++) {
            out[i] = (int64_t)(p0[i] + p2[i]);
            out[i + m] = (int64_t)(p1[i] - p2[i]);
        }
    }
}





void prod_pol_mat_toeplitz_recursive_i128(int n, __int128 *out, const __int128 *vector, const int64_t *toeplitz_matrix){
    if (n < TOEPLITZ_THRESHOLD || ((n%3 !=0) && ((n & 1) == 1)))  {
        small_toeplitz_vector_matrix_i128(n, vector, toeplitz_matrix, out);
        return;
    }

    else if (n % 3 == 0) {
        int m = n / 3;

        // Decompose the input vector into three parts
        const __int128 *v0 = vector;
        const __int128 *v1 = vector + m;
        const __int128 *v2 = vector + 2 * m;

        // Additional vectors for the differences of the three parts of the input vector
        __int128 v_sub_02[m], v_sub_12[m], v_sub_01[m];
        for (int i = 0; i < m; i++) {
            v_sub_02[i] = v0[i] - v2[i];
            v_sub_12[i] = v1[i] - v2[i];
            v_sub_01[i] = v0[i] - v1[i];
        }

        // Decompose the Toeplitz matrix into five parts corresponding to the three parts of the input vector
        const int64_t *M0 = toeplitz_matrix;
        const int64_t *M1 = toeplitz_matrix + m;
        const int64_t *M2 = toeplitz_matrix + 2 * m;
        const int64_t *M3 = toeplitz_matrix + 3 * m;
        const int64_t *M4 = toeplitz_matrix + 4 * m;

        // Preparation of matrix sums (size of a sub-matrix : 2m - 1)
        int64_t M_sum_012[2 * m - 1], M_sum_123[2 * m - 1], M_sum_234[2 * m - 1];
        for (int i = 0; i < 2 * m - 1; i++) {
            M_sum_012[i] = M0[i] + M1[i] + M2[i];
            M_sum_123[i] = M1[i] + M2[i] + M3[i];
            M_sum_234[i] = M2[i] + M3[i] + M4[i];
        }

        // Compute the products of the sub-vectors with the corresponding sub-matrices using recursive calls
        __int128 p0[m], p1[m], p2[m], p3[m], p4[m], p5[m];
        prod_pol_mat_toeplitz_recursive_i128(m, p0, v2, M_sum_012);
        prod_pol_mat_toeplitz_recursive_i128(m, p1, v1, M_sum_123);
        prod_pol_mat_toeplitz_recursive_i128(m, p2, v0, M_sum_234);
        prod_pol_mat_toeplitz_recursive_i128(m, p3, v_sub_02, M2);
        prod_pol_mat_toeplitz_recursive_i128(m, p4, v_sub_12, M1);
        prod_pol_mat_toeplitz_recursive_i128(m, p5, v_sub_01, M3);

        // Recombine the results of the sub-products to obtain the final result
        for (int i = 0; i < m; i++) {
            out[i] = p0[i] + p3[i] + p4[i];
            out[i + m] = p1[i] + p5[i] - p4[i];
            out[i + 2 * m] = p2[i] - p3[i] - p5[i];
        }
    }

    else if ((n & 1) == 0) {
        int m = n / 2;
        const __int128 *v0 = vector;
        const __int128 *v1 = vector + m;
            
        __int128 v_sub_01[m];
        for (int i = 0; i < m; i++)
            v_sub_01[i] = v0[i] - v1[i];

        const int64_t *M0 = toeplitz_matrix;
        const int64_t *M1 = toeplitz_matrix + m;
        const int64_t *M2 = toeplitz_matrix + 2 * m;

        int64_t T_sum_M01[2 * m - 1];
        int64_t T_sum_M12[2 * m - 1];

        for (int i = 0; i < 2 * m - 1; i++) {
            T_sum_M01[i] = M0[i] + M1[i];
            T_sum_M12[i] = M1[i] + M2[i];
        }

        __int128 p0[m], p1[m], p2[m];

        prod_pol_mat_toeplitz_recursive_i128(m, p0, v1, T_sum_M01);
        prod_pol_mat_toeplitz_recursive_i128(m, p1, v0, T_sum_M12);
        prod_pol_mat_toeplitz_recursive_i128(m, p2, v_sub_01, M1);

        for (int i = 0; i < m; i++) {
            out[i] = p0[i] + p2[i];
            out[i + m] = p1[i] - p2[i];
        }
    }
}
# endif




#if IS_BABAI_USABLE
void coeff_shift_i64(int64_t out[DEGREE], __int128 polynomial[DEGREE], int n_shift){
    for (int i=0; i<DEGREE; i++)
        out[i] = (int64_t)(polynomial[i] >> n_shift);
}

void coeff_shift_i128(__int128 out[DEGREE], __int128 polynomial[DEGREE], int n_shift){
    for (int i=0; i<DEGREE; i++)
        out[i] = polynomial[i] >> n_shift;
}
#endif

#if IS_DOUBLE_SPARSE
void linear_prod_pol_lattice_inv_i64(int64_t out[DEGREE], __int128 polynomial[DEGREE]){
    for (int i=0; i<EXTENSION_DEGREE; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + EXTENSION_DEGREE];

        out[i] = (int64_t)(coeff_1 * LAMBDA_INV_MOD + GAMMA_POW_N_LAMBDA_MOD * coeff_2);
    }

    for (int i=EXTENSION_DEGREE; i<DEGREE - EXTENSION_DEGREE; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + EXTENSION_DEGREE];

        out[i] = (int64_t)((-GAMMA_POW_K_MOD * coeff_2) - coeff_1);
    }

    for (int i=DEGREE - EXTENSION_DEGREE; i<DEGREE; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + EXTENSION_DEGREE - DEGREE];

        out[i] = (int64_t)(-coeff_1 - GAMMA_POW_N_LAMBDA_MOD * coeff_2);
    }
}

void linear_prod_pol_lattice_i128(__int128 out[DEGREE], int64_t polynomial[DEGREE]){
    for (int i=0; i<EXTENSION_DEGREE; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + EXTENSION_DEGREE];

        out[i] = -LAMBDA * coeff_1 - GAMMA_POW_K_MOD * coeff_2;
    }

    for (int i=EXTENSION_DEGREE; i<DEGREE - EXTENSION_DEGREE; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + EXTENSION_DEGREE];

        out[i] = (__int128)(coeff_1 - GAMMA_POW_K_MOD * coeff_2);
    }

    for (int i=DEGREE - EXTENSION_DEGREE; i<DEGREE; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + EXTENSION_DEGREE - DEGREE];

        out[i] = (__int128)(coeff_1 + (__int128)(GAMMA_POW_K_MOD) * coeff_2);
    }
}
#endif



