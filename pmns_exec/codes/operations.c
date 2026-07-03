#include "./interfaces/operations_interface.h"
#include "./interfaces/reductions_interface.h"
#include "../params/pmns_params.h"
#include <string.h>

#define TOEPLITZ_THRESHOLD 50
#define KARATSUBA_THRESHOLD 30

/*=================================================================
                FUNCTIONS OVER POLYNOMIALS COEFFICIENTS
=================================================================*/
# if IS_BABAI_USABLE
void coeff_shift_i64(int degree, __int128 out[degree], __int128 polynomial[degree], int n_shift){
    /*define a function to shift coefficients of a polynomial by n_shift bits to the right, using int128 type for the coefficients
    and cast the result to int64_t*/
    for (int i=0; i<degree; i++)
        out[i] = (int64_t)(polynomial[i] >> n_shift);
}

void coeff_shift_i128(int degree, __int128 out[degree], __int128 polynomial[degree], int n_shift){
    /*define a function to shift coefficients of a polynomial by n_shift bits to the right, using int128 type for the coefficients*/
    for (int i=0; i<degree; i++)
        out[i] = polynomial[i] >> n_shift;
}
#endif

/*=================================================================
        FUNCTIONS FOR PRODUCTS WITH POLYNOMIALS AND INTEGERS 
=================================================================*/
void addmul_pol64_int64(int degree, __int128 out[degree], const int64_t polynomial[degree], int64_t coeff){
    // compute the product of an int64_t polynomial with an int64_t coefficient and add the result to an int128 polynomial
    if (coeff == 0) return;

    for (int deg=0; deg<degree; deg++)
        out[deg] += (__int128)polynomial[deg] * (__int128)coeff;
}


void addmul_pol64_int64_mpn(int degree, int n_limbs, mp_limb_t out[degree][n_limbs], const int64_t polynomial[degree], int64_t coeff){
    /*compute the product of an int64_t polynomial with an int64_t coefficient and add the result to an mp_limb_t polynomial*/
    if (coeff == 0) return;

    for (int deg=0; deg<degree; deg++){
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


/*=================================================================
        FUNCTIONS FOR PRODUCTS WITH POLYNOMIALS AND POLYNOMIALS 
=================================================================*/

void addmul_polmpn_pol64(int degree, int n_limbs, mp_limb_t out[degree][n_limbs], mp_limb_t pol_mpn[degree][n_limbs], const int64_t pol64[degree]) {
    // compute the product of an int128 polynomial with an int64_t coefficient and add the result to an mp_limb_t polynomial
    // use complement to 2 representation for negative values
    int PROD_SIZE = 2 * degree - 1;
    mp_limb_t tmp_pol[PROD_SIZE][n_limbs];

    for (int i=0; i<PROD_SIZE; i++) mpn_zero(tmp_pol[i], n_limbs);

    // polynomial product
    for (int i = 0; i < degree; i++) {
        for (int j = 0; j < degree; j++) {
            if (pol64[j] == 0) continue;

            int sign = (pol64[j] < 0);
            uint64_t abs_c = sign ? (uint64_t)(-(uint64_t)pol64[j]) : (uint64_t)pol64[j];
            (sign ? mpn_submul_1 : mpn_addmul_1)(tmp_pol[i+j], pol_mpn[i], n_limbs, (mp_limb_t)abs_c);
        }
    }

    for (int deg=degree; deg<PROD_SIZE; deg++){
        int affect_deg = deg - degree;

        for (int i=0; i<degree; i++) {
            int64_t mat_coeff = EXT_MAT[affect_deg][i];
            if (mat_coeff == 0) continue;


            int sign = (mat_coeff < 0);
            uint64_t abs_coeff = sign ? (uint64_t)(-(uint64_t)mat_coeff) : (uint64_t)mat_coeff;

            (sign ? mpn_submul_1 : mpn_addmul_1)(out[i], tmp_pol[deg], n_limbs, (mp_limb_t)abs_coeff);
        }
    }

    for (int i = 0; i < degree; i++)
        mpn_add_n(out[i], out[i], tmp_pol[i], n_limbs);
}


static inline void recursive_karatsuba_product(int degree, __int128 out[2 * degree - 1], int64_t pol_A[degree], int64_t pol_B[degree]) {
    // Compute the product of two polynomials using schoolbook multiplication
    if (degree < KARATSUBA_THRESHOLD || ((degree & 1) == 1)) {
        for (int i = 0; i < 2 * degree - 1; i++) out[i] = 0;

        for (int i = 0; i < degree; i++) {
            __int128 ai = (__int128)pol_A[i];
            for (int j = 0; j < degree; j++) {
                out[i + j] += ai * pol_B[j];
            }
        }
        return;
    }

    int m = degree / 2;

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
    for (int i = 0; i < 2 * degree - 1; i++) out[i] = 0;

    for (int i = 0; i < 2 * m - 1; i++) {
        out[i] += Y0[i];
        out[i + m] += Y1[i] - Y0[i] - Y2[i];
        out[i + 2 * m] += Y2[i];
    }
}


void polynomials_product(int degree, __int128 out[degree], int64_t PolA[degree], int64_t PolB[degree]){
    int PROD_SIZE = 2 * degree - 1;
    __int128 tmp_poly[PROD_SIZE];

    #if IS_E_BINOMIAL
    __int128 sum;            
    for(int i = 0; i < DEGREE-1; i++){
        tmp_poly[i + DEGREE - 1] = PolB[i];
        tmp_poly[i] = PolB[1 + i] * LAMBDA;
    }

    tmp_poly[2*DEGREE - 2] = PolB[DEGREE - 1];

    for(int i = 0; i < DEGREE; i++){
        sum = 0;
        for(int j = 0; j < DEGREE; j++)
            sum += (__int128) PolA[j] * tmp_poly[DEGREE - 1 - j + i];
        out[i] = sum;
    }
    #else
    recursive_karatsuba_product(degree, tmp_poly, PolA, PolB);

    // Directly copy the first degree coefficients to the output polynomial
    for (int i = 0; i < degree; i++) out[i] = tmp_poly[i];
    
    // Apply generic external reduction using the external reduction matrix for the remaining coefficients
    for (int j = 0; j < degree - 1; j++) {
        __int128 p_val = tmp_poly[degree + j];

        for (int i = 0; i < degree; i++) {
            out[i] += p_val * (__int128)EXT_MAT[j][i];
        }
    }
    #endif
}

# if IS_ELEMENTS_IN_GAMMA_BASIS
// /!\ this function only works if E = X^n - lambda
// multiply polynomial by X^pow and apply external reduction if needed, then add the result to out polynomial    
void addmul_polmpn_Xpow_modE(int degree, int n_limbs, mp_limb_t out[degree][n_limbs], mp_limb_t polmpn[degree][n_limbs], unsigned int pow) {
    // extract lambda from the external reduction matrix for the given power
    const int lambda_sign = (LAMBDA < 0);
    const uint64_t abs_lambda = lambda_sign ? (uint64_t)(-LAMBDA) : (uint64_t)LAMBDA;

    // multiply by X^pow consist on a register shift, then multiply coefficients by lambda to apply external reduction if needed
    // the result is added to the output polynomial. We use complement to 2 representation for negative values and directly apply
    // the external reduction during the shift if needed to avoid extra computations.
    for (int deg=0; deg<degree; deg++){
        // compute affectation degree
        int curr_deg = deg + pow;
        int prod_sign = 0;

        // correct affectation degree and sign of the product if external reduction is needed
        if (curr_deg >= degree){
            mpn_mul_1(polmpn[deg], polmpn[deg], n_limbs, abs_lambda);

            curr_deg -= degree;
            prod_sign ^= lambda_sign;
        }
        // add the product to the output polynomial with the correct sign
        (prod_sign ? mpn_sub_n : mpn_add_n)(out[curr_deg], out[curr_deg], polmpn[deg], n_limbs);
    }
}
#endif


/*=================================================================
        FUNCTIONS FOR PRODUCTS WITH POLYNOMIALS AND MATRICES 
=================================================================*/
void prod_polmpn_mat_i64(int degree, int n_limbs, int64_t out[degree], mp_limb_t pol[degree][n_limbs], const int64_t sublattice[degree][degree]) {
    // compute the product of an mp_limb_t polynomial with an int64_t matrix and store the result in an int64_t polynomial
    int SIZE = n_limbs + 1;
    mp_limb_t acc[SIZE];
    mp_limb_t prod[SIZE];

    for (int j = 0; j < DEGREE; j++) {
        mpn_zero(acc, SIZE);

        for (int i = 0; i < DEGREE; i++) {
            int64_t mat_coeff = sublattice[i][j];
            if (mat_coeff == 0) continue;

            int mat_sign = (mat_coeff < 0);
            uint64_t abs_c = mat_sign ? (uint64_t)(-mat_coeff) : (uint64_t)mat_coeff;
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
        out[j] = (int64_t)acc[0];
    }

}


void prod_pol_mat_mpn(int degree, int n_limbs, mp_limb_t out[degree][n_limbs+1], int64_t pol[degree], const int64_t sublattice[degree][degree]){
    // compute the product of an int64_t polynomial with an int64_t matrix and store the result in an mp_limb_t polynomial
    int SIZE = n_limbs + 1;
    for (int j = 0; j < degree; j++) {
        mpn_zero(out[j], SIZE);

        for (int i = 0; i < degree; i++) {
            int64_t mat_coeff = sublattice[i][j];
            int64_t pol_coeff = pol[i];
            if (mat_coeff == 0 || pol_coeff == 0) continue;
            
            int mat_sign = (mat_coeff < 0);
            int pol_sign   = (pol_coeff < 0);
            uint64_t mat_abs = mat_sign ? -(uint64_t)mat_coeff : (uint64_t)mat_coeff;
            uint64_t pol_abs   = pol_sign   ? (uint64_t)(-pol_coeff) : (uint64_t)pol_coeff;

            mp_limb_t mat_limb = (mp_limb_t)mat_abs;
            mp_limb_t pol_limb   = (mp_limb_t)pol_abs;

            mp_limb_t prod_limbs[SIZE];
            mpn_zero(prod_limbs, SIZE);

            prod_limbs[1] = mpn_mul_1(prod_limbs, &pol_limb, 1, mat_limb);

            int sign_prod = mat_sign ^ pol_sign;
            (sign_prod ? mpn_sub_n : mpn_add_n)(out[j], out[j], prod_limbs, SIZE);
        }
    }
}


void prod_pol_mat_mpz(int degree, mpz_t out[degree], int64_t pol[degree], const int64_t sublattice[degree][degree]) {
    // compute the product of an int64_t polynomial with an int64_t matrix and store the result in an mpz_t polynomial
    mpz_t tmp;
    mpz_init(tmp);

    for (int j = 0; j < degree; j++) {
        mpz_init_set_ui(out[j], 0);

        for (int i = 0; i < degree; i++) {
            int64_t mat_coeff = sublattice[i][j];
            int64_t pol_coeff = pol[i];
            if (mat_coeff == 0 || pol_coeff == 0) continue;
            
            __int128 prod = (__int128)mat_coeff * (__int128)pol_coeff;
            unsigned __int128 abs_prod = prod < 0 ? -prod : prod;

            mpz_import(tmp, 2, -1, sizeof(uint64_t), 0, 0, &abs_prod);
            if (prod < 0) mpz_neg(tmp, tmp);
            mpz_add(out[j], out[j], tmp);
        }
    }

    mpz_clear(tmp);
}


void prod_polmpz_mat_i64(int degree, int64_t out[degree], mpz_t pol[degree], const int64_t sublattice[degree][degree], mpz_t phi) {
    // compute the product of an mpz_t polynomial with an int64_t matrix and store the result in an int64_t polynomial
    mpz_t tmp, acc;
    mpz_inits(tmp, acc, NULL);

     for (int j = 0; j < DEGREE; j++) {
        mpz_set_ui(acc, 0);

        for (int i = 0; i < DEGREE; i++) {
            mpz_mul_ui(tmp, pol[i], sublattice[i][j]);
            mpz_add(acc, acc, tmp);
        }
        mpz_mod(acc, acc, phi);
        out[j] = mpz_get_ui(acc);
    }
    mpz_clears(tmp, acc, NULL);
}


void prod_pol_mat_i64(int degree, __int128 out[degree], __int128 polynomial[degree], const int64_t matrix[degree][degree]) {
    // compute the product of an int128 polynomial with an int64_t matrix mod 2^PHI_POW and store the result in an int128 polynomial
    for (int i=0; i<degree; i++){
        int64_t acc = 0;
        for (int j=0; j<degree; j++)
            acc += (int64_t)polynomial[j] * (int64_t)matrix[j][i];
        out[i] = acc;
    }
}


void prod_pol_mat_i128(int degree, __int128 out[degree], __int128 polynomial[degree], const int64_t matrix[degree][degree]) {
    // compute the product of an int128 polynomial with an int64_t matrix and store the result in an int128 polynomial
    for (int i=0; i<degree; i++){
        __int128 acc = 0;
        for (int j=0; j<degree; j++)
            acc += (__int128)polynomial[j] * (__int128)matrix[j][i];
        out[i] = acc;
    }
}


# if IS_TOEPLITZ_USABLE
void prod_pol_mat_toeplitz_i64(int degree, __int128 out[degree], const __int128 polynomial[degree], const uint64_t matrix_toeplitz[2*degree - 1]) {
    // compute the product of an int128 polynomial with a uint64_t Toeplitz matrix mod 2^PHI_POW and store the result in an int128 polynomial
    for (int i = 0; i < degree; i++) {
        int64_t acc = 0;
        for (int j = 0; j < degree; j++)
            acc += (int64_t)polynomial[j] * (int64_t)matrix_toeplitz[degree - 1 - j + i];
        out[i] = acc;
    }
}


void prod_pol_mat_toeplitz_i128(int degree, __int128 out[degree], const __int128 polynomial[degree], const int64_t matrix_toeplitz[2*degree - 1]) {
    // compute the product of an int128 polynomial with an int64_t Toeplitz matrix and store the result in an int128 polynomial
    for (int i = 0; i < degree; i++) {
        __int128 acc = 0;
        for (int j = 0; j < degree; j++)
            acc += (__int128)polynomial[j] * (__int128)matrix_toeplitz[degree - 1 - j + i];
        out[i] = acc;
    }
}


static inline void small_prod_pol_mat_toeplitz_i64(int degree, __int128 out[degree], const __int128 polynomial[degree], const uint64_t matrix_toeplitz[2*degree - 1]){
    for (int i = 0; i < degree; i++) {
        int64_t acc = 0;
        for (int j = 0; j < degree; j++)
            acc += (int64_t)polynomial[j] * (int64_t)matrix_toeplitz[degree - 1 - j + i];
        out[i] = acc;
    }
}


static inline void small_prod_pol_mat_toeplitz_i128(int degree, __int128 out[degree], const __int128 polynomial[degree], const int64_t matrix_toeplitz[2*degree - 1]){
    for (int i = 0; i < degree; i++) {
        __int128 acc = 0;
        for (int j = 0; j < degree; j++)
            acc += (__int128)polynomial[j] * (__int128)matrix_toeplitz[degree - 1 - j + i];
        out[i] = acc;
    }
}


static inline void internal_prod_pol_mat_toeplitz_recursive_i64(int degree, __int128 out[degree], const __int128 vector[degree], const uint64_t toeplitz_matrix[2*degree - 1]){
    if ((degree < TOEPLITZ_THRESHOLD) || ((degree % 3 != 0) && (!(degree & 1)))){
        small_prod_pol_mat_toeplitz_i64(degree, out, vector, toeplitz_matrix);
        return;                                                           
    }                                                                     
                                                                          
    else if (degree % 3 == 0) {                                         
        int m = degree / 3;                                             
                                                                          
        const __int128 *v0 = vector;                                 
        const __int128 *v1 = vector + m;                              
        const __int128 *v2 = vector + 2 * m;                          
                                                                          
        __int128 v_sub_02[m], v_sub_12[m], v_sub_01[m];                   
        for (int i = 0; i < m; i++) {                                     
            v_sub_02[i] = v0[i] - v2[i];                                  
            v_sub_12[i] = v1[i] - v2[i];                                  
            v_sub_01[i] = v0[i] - v1[i];                                  
        }                                                                
                                                                          
        const uint64_t *M0 = toeplitz_matrix;                             
        const uint64_t *M1 = toeplitz_matrix + m;                         
        const uint64_t *M2 = toeplitz_matrix + 2 * m;                     
        const uint64_t *M3 = toeplitz_matrix + 3 * m;                     
        const uint64_t *M4 = toeplitz_matrix + 4 * m;                     
                                                                          
        uint64_t M_sum_012[2 * m - 1];                                    
        uint64_t M_sum_123[2 * m - 1];                                    
        uint64_t M_sum_234[2 * m - 1];                                    
                                                                          
        for (int i = 0; i < 2 * m - 1; i++) {                             
            M_sum_012[i] = M0[i] + M1[i] + M2[i];                         
            M_sum_123[i] = M1[i] + M2[i] + M3[i];                         
            M_sum_234[i] = M2[i] + M3[i] + M4[i];                         
        }                                                                
                                                                          
        __int128 p0[m], p1[m], p2[m], p3[m], p4[m], p5[m];               
                                                                          
        internal_prod_pol_mat_toeplitz_recursive_i64(m, p0, v2, M_sum_012);                             
        internal_prod_pol_mat_toeplitz_recursive_i64(m, p1, v1, M_sum_123);                             
        internal_prod_pol_mat_toeplitz_recursive_i64(m, p2, v0, M_sum_234);                             
        internal_prod_pol_mat_toeplitz_recursive_i64(m, p3, v_sub_02, M2);                              
        internal_prod_pol_mat_toeplitz_recursive_i64(m, p4, v_sub_12, M1);                              
        internal_prod_pol_mat_toeplitz_recursive_i64(m, p5, v_sub_01, M3);                              
                                                                          
        for (int i = 0; i < m; i++) {                                    
            out[i] = (int64_t)(p0[i] + p3[i] + p4[i]);        
            out[i + m] = (int64_t)(p1[i] + p5[i] - p4[i]);        
            out[i + 2 * m] = (int64_t)(p2[i] - p3[i] - p5[i]);        
        }                                                                
        return;                                                          
    }                                                                     
                                                                          
    else if (!(degree & 1)) {                                                                
        int m = (degree) / 2;                                             
                                                                          
        const __int128 *v0 = vector;                                 
        const __int128 *v1 = vector + m;                              
                                                                          
        __int128 v_sub_01[m];                                             
        for (int i = 0; i < m; i++)                                       
            v_sub_01[i] = v0[i] - v1[i];                                  
                                                                          
        const uint64_t *M0 = toeplitz_matrix;                             
        const uint64_t *M1 = toeplitz_matrix + m;                         
        const uint64_t *M2 = toeplitz_matrix + 2 * m;                     
                                                                          
        uint64_t M_sum_M01[2 * m - 1];                                    
        uint64_t M_sum_M12[2 * m - 1];                                    
                                                                          
        for (int i = 0; i < 2 * m - 1; i++) {                             
            M_sum_M01[i] = M0[i] + M1[i];                                 
            M_sum_M12[i] = M1[i] + M2[i];                                 
        }                                                                
                                                                        
        __int128 p0[m], p1[m], p2[m];                                     
                                                                        
        internal_prod_pol_mat_toeplitz_recursive_i64(m, p0, v1, M_sum_M01);                             
        internal_prod_pol_mat_toeplitz_recursive_i64(m, p1, v0, M_sum_M12);                             
        internal_prod_pol_mat_toeplitz_recursive_i64(m, p2, v_sub_01, M1);                              
                                                                        
        for (int i = 0; i < m; i++) {                                     
            out[i] = (int64_t)(p0[i] + p2[i]);                    
            out[i + m] = (int64_t)(p1[i] - p2[i]);                    
        }                                                                
        return;                                                           
    }
}

static inline void internal_prod_pol_mat_toeplitz_recursive_i128(int degree, __int128 out[degree], const __int128 vector[degree], const int64_t toeplitz_matrix[2*degree - 1]){
    if ((degree < TOEPLITZ_THRESHOLD) || ((degree % 3 != 0) && (!(degree & 1)))){
        small_prod_pol_mat_toeplitz_i128(degree, out, vector, toeplitz_matrix);
        return;                                                           
    }                                                                     
                                                                          
    else if (degree % 3 == 0) {                                         
        int m = degree / 3;                                             
                                                                          
        const __int128 *v0 = vector;                                 
        const __int128 *v1 = vector + m;                              
        const __int128 *v2 = vector + 2 * m;                          
                                                                          
        __int128 v_sub_02[m], v_sub_12[m], v_sub_01[m];                   
        for (int i = 0; i < m; i++) {                                     
            v_sub_02[i] = v0[i] - v2[i];                                  
            v_sub_12[i] = v1[i] - v2[i];                                  
            v_sub_01[i] = v0[i] - v1[i];                                  
        }                                                                
                                                                          
        const int64_t *M0 = toeplitz_matrix;                             
        const int64_t *M1 = toeplitz_matrix + m;                         
        const int64_t *M2 = toeplitz_matrix + 2 * m;                     
        const int64_t *M3 = toeplitz_matrix + 3 * m;                     
        const int64_t *M4 = toeplitz_matrix + 4 * m;                     
                                                                          
        int64_t M_sum_012[2 * m - 1];                                    
        int64_t M_sum_123[2 * m - 1];                                    
        int64_t M_sum_234[2 * m - 1];                                    
                                                                          
        for (int i = 0; i < 2 * m - 1; i++) {                             
            M_sum_012[i] = M0[i] + M1[i] + M2[i];                         
            M_sum_123[i] = M1[i] + M2[i] + M3[i];                         
            M_sum_234[i] = M2[i] + M3[i] + M4[i];                         
        }                                                                
                                                                          
        __int128 p0[m], p1[m], p2[m], p3[m], p4[m], p5[m];               
                                                                          
        internal_prod_pol_mat_toeplitz_recursive_i128(m, p0, v2, M_sum_012);                             
        internal_prod_pol_mat_toeplitz_recursive_i128(m, p1, v1, M_sum_123);                             
        internal_prod_pol_mat_toeplitz_recursive_i128(m, p2, v0, M_sum_234);                             
        internal_prod_pol_mat_toeplitz_recursive_i128(m, p3, v_sub_02, M2);                              
        internal_prod_pol_mat_toeplitz_recursive_i128(m, p4, v_sub_12, M1);                              
        internal_prod_pol_mat_toeplitz_recursive_i128(m, p5, v_sub_01, M3);                              
                                                                          
        for (int i = 0; i < m; i++) {                                    
            out[i] = (__int128)(p0[i] + p3[i] + p4[i]);        
            out[i + m] = (__int128)(p1[i] + p5[i] - p4[i]);        
            out[i + 2 * m] = (__int128)(p2[i] - p3[i] - p5[i]);        
        }                                                                
        return;                                                          
    }                                                                     
                                                                          
    else if (!(degree & 1)) {                                                                
        int m = (degree) / 2;                                             
                                                                          
        const __int128 *v0 = vector;                                 
        const __int128 *v1 = vector + m;                              
                                                                          
        __int128 v_sub_01[m];                                             
        for (int i = 0; i < m; i++)                                       
            v_sub_01[i] = v0[i] - v1[i];                                  
                                                                          
        const int64_t *M0 = toeplitz_matrix;                             
        const int64_t *M1 = toeplitz_matrix + m;                         
        const int64_t *M2 = toeplitz_matrix + 2 * m;                     
                                                                          
        int64_t M_sum_M01[2 * m - 1];                                    
        int64_t M_sum_M12[2 * m - 1];                                    
                                                                          
        for (int i = 0; i < 2 * m - 1; i++) {                             
            M_sum_M01[i] = M0[i] + M1[i];                                 
            M_sum_M12[i] = M1[i] + M2[i];                                 
        }                                                                
                                                                        
        __int128 p0[m], p1[m], p2[m];                                     
                                                                        
        internal_prod_pol_mat_toeplitz_recursive_i128(m, p0, v1, M_sum_M01);                             
        internal_prod_pol_mat_toeplitz_recursive_i128(m, p1, v0, M_sum_M12);                             
        internal_prod_pol_mat_toeplitz_recursive_i128(m, p2, v_sub_01, M1);                              
                                                                        
        for (int i = 0; i < m; i++) {                                     
            out[i] = (int64_t)(p0[i] + p2[i]);                    
            out[i + m] = (int64_t)(p1[i] - p2[i]);                    
        }                                                                
        return;                                                           
    }
}
void prod_pol_mat_toeplitz_recursive_i64(int degree, __int128 out[degree], const __int128 vector[degree], const uint64_t toeplitz_matrix[2*degree - 1]){
    // compute the product of an int128 polynomial with a uint64_t Toeplitz matrix and store the result in an int128 polynomial using a recursive approach
    internal_prod_pol_mat_toeplitz_recursive_i64(degree, out, vector, toeplitz_matrix);
}


void prod_pol_mat_toeplitz_recursive_i128(int degree, __int128 out[degree], const __int128 vector[degree], const int64_t toeplitz_matrix[2*degree - 1]){
    // compute the product of an int128 polynomial with an int64_t Toeplitz matrix and store the result in an int128 polynomial using a recursive approach
    internal_prod_pol_mat_toeplitz_recursive_i128(degree, out, vector, toeplitz_matrix);
}
# endif


# if IS_DOUBLE_SPARSE
void prod_pol_mat_linear_i64(int extension_degree, int degree, int64_t out[degree], __int128 polynomial[degree]) {
    /* compute the product of an int128 polynomial with a linear matrix and store the result in an int64_t polynomial
    This function is specialized for the double sparse structure */
    for (int i=0; i<extension_degree; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + extension_degree];

        out[i] = (int64_t)(coeff_1 * LAMBDA_INV_MOD + GAMMA_POW_N_LAMBDA_MOD * coeff_2);
    }

    for (int i=extension_degree; i<degree - extension_degree; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + extension_degree];

        out[i] = (int64_t)((-GAMMA_POW_K_MOD * coeff_2) - coeff_1);
    }

    for (int i=degree - extension_degree; i<degree; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + extension_degree - degree];

        out[i] = (int64_t)(-coeff_1 - GAMMA_POW_N_LAMBDA_MOD * coeff_2);
    }
}

void prod_pol_mat_linear_i128(int extension_degree, int degree, __int128 out[degree], int64_t polynomial[degree]){
    /* compute the product of an int64_t polynomial with a linear matrix and store the result in an int128 polynomial.
     This function is specialized for the double sparse structure */

    for (int i=0; i<extension_degree; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + extension_degree];

        out[i] = -LAMBDA * coeff_1 - GAMMA_POW_K_MOD * coeff_2;
    }

    for (int i=extension_degree; i<degree - extension_degree; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + extension_degree];

        out[i] = (__int128)(coeff_1 - GAMMA_POW_K_MOD * coeff_2);
    }

    for (int i=degree - extension_degree; i<degree; i++){
        __int128 coeff_1 = polynomial[i];
        __int128 coeff_2 = polynomial[i + extension_degree - degree];

        out[i] = (__int128)(coeff_1 + (__int128)(GAMMA_POW_K_MOD) * coeff_2);
    }
}
#endif
