#include "./interfaces/operations_interface.h"
#include "./interfaces/reductions_interface.h"
#include "../params/pmns_params.h"

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

void prod_pol_mat_toeplitz_i64(int64_t out[DEGREE], __int128 polynomial[DEGREE], const uint64_t matrix_toeplitz[2*DEGREE - 1]) {
    for (int i = 0; i < DEGREE; i++) {
        __int128 acc = 0;
        for (int j = 0; j < DEGREE; j++) {
            acc += (__int128)polynomial[j] * (__int128)matrix_toeplitz[DEGREE - 1 - j + i];
        }
        out[i] = (uint64_t)acc;
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

void polynomials_product(__int128 out[DEGREE], int64_t PolA[DEGREE], int64_t PolB[DEGREE]){
    const int PROD_SIZE = 2*DEGREE - 1;
    __int128 P[PROD_SIZE];

    for (int deg=0; deg<PROD_SIZE; deg++) P[deg]=0;

    for (int degA=0; degA<DEGREE; degA++){
        __int128 ai = PolA[degA];
        for (int degB=0; degB<DEGREE; degB++)
            P[degA + degB] += ai * PolB[degB];
    }

    for (int i=0; i<DEGREE; i++) out[i] = P[i];

    for (int i=0; i<DEGREE; i++){
        __int128 acc = 0;
        for (int j = 0; j < DEGREE-1; j++)
            acc += (__int128)P[DEGREE + j] * (__int128)EXT_MAT[j][i];
        out[i] += acc;
    }
}

#ifndef IS_SPARSE
void coeff_shift_i64(int64_t out[DEGREE], __int128 polynomial[DEGREE], int n_shift){
    for (int i=0; i<DEGREE; i++)
        out[i] = (int64_t)(polynomial[i] >> n_shift);
}

void coeff_shift_i128(__int128 out[DEGREE], __int128 polynomial[DEGREE], int n_shift){
    for (int i=0; i<DEGREE; i++)
        out[i] = polynomial[i] >> n_shift;
}
#else

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