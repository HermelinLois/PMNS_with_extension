# include "../codes/interfaces/conversions_interface.h"
# include "../codes/interfaces/reductions_interface.h"
# include "../codes/interfaces/operations_interface.h"
# include "../params/conversions_params.h"
# include "../params/reductions_params.h"
# include <stdio.h>
# include <gmp.h>

/*=================================================================
                FUNCTIONS USED FOR CONVERSIONS 
=================================================================*/
static inline uint64_t apply_mask(const mp_limb_t element_data[N_LIMBS], uint64_t mask_size, size_t *mask_start_pos) {
    // apply a mask of size mask_size on the element data starting from the position given by mask_start_pos and 
    // return the result as a uint64_t, then increase the mask_start_pos by mask_size
    uint64_t start = *mask_start_pos;
    uint64_t limb_idx = start / GMP_NUMB_BITS;

    uint64_t mask = (1ULL << mask_size) - 1ULL;
    uint64_t shift = start % GMP_NUMB_BITS;
    uint64_t chunk = 0;

    // apply mask which is lower than phi_po
    if (limb_idx < N_LIMBS) {
        chunk = (uint64_t)(element_data[limb_idx] >> shift);

        // apply a second mask on higher register in case the mask apply on 2 register at the same time
        if (mask_size > (GMP_NUMB_BITS - shift) && (limb_idx + 1 < N_LIMBS)) {
            chunk |= (uint64_t)(element_data[limb_idx + 1] << (GMP_NUMB_BITS - shift));
        }
    }
    // increase starting point using given mask_size
    *mask_start_pos += mask_size;
    return chunk & mask;
}


static inline void init_element_times_phipow(int extension_degree, int degree, mpz_t out[degree], const mp_limb_t data[extension_degree][N_LIMBS], unsigned int power) {
    // init mpz_t with element data and multiply the values by phi^power mod p
    mpz_t prime, coeff, acc;
    mpz_inits(prime, coeff, acc, NULL);
    mpz_import(prime, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, P);
    
    for (int deg = 0; deg < extension_degree; deg++) {
        mpz_init(out[deg]);
        mpz_import(coeff, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, data[deg]);
        mpz_mul_2exp(coeff, coeff, PHI_POW * power);
        mpz_mod(out[deg], coeff, prime);
    }
    
    for (int j=extension_degree; j<degree; j++) mpz_init_set_ui(out[j], 0);
    mpz_clears(prime, coeff, acc, NULL);
}


# if !IS_ELEMENTS_IN_GAMMA_BASIS
static inline void element_change_basis(int extension_degree, int degree, mpz_t out[degree], mpz_t extension_element[degree]){
    // compute the product of the given element with the transition matrix and apply modulus p to the result
    // with our implementation, element are express in gamma basis, no need to convert element
    mpz_t acc, coeff, prime, tmp;
    mpz_inits(acc, coeff, prime, tmp, NULL);
    mpz_import(prime, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, P);

    for (int j=0; j<extension_degree; j++){
        mpz_init_set_ui(tmp, 0);

        for (int i=0; i<extension_degree; i++){
            mpz_import(coeff, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, TRANSITION_MATRIX[i][j]);
            mpz_mul(coeff, extension_element[i], coeff);
            mpz_add(tmp, tmp, coeff);
        }
        mpz_mod(out[j], tmp, prime);
    }
    mpz_clears(acc, coeff, prime, tmp, NULL);
}
# endif


/*=================================================================
                        CONVERSIONS FUNCTIONS
=================================================================*/
void convert_element_to_pmns_exact(int extension_degree, int degree, int64_t out[degree], const mp_limb_t element_data[extension_degree][N_LIMBS]){
    // convert an element of the extension field given to its representation in PMNS using the exact method
    // This method is slower than the fast method but does not require precomputation of the PMNS representation of the basis elements
    mpz_t vector[degree];
    
    init_element_times_phipow(extension_degree, degree, vector, element_data, N_INT_RED_CLASSICAL);

    # if !IS_ELEMENTS_IN_GAMMA_BASIS
    mpz_t intermediate_vector[degree];
    for (int deg=0; deg<degree; deg++) mpz_init(intermediate_vector[deg]);

    element_change_basis(extension_degree, degree, intermediate_vector, vector);

    for (int deg=0; deg<degree; deg++) {
        mpz_set(vector[deg], intermediate_vector[deg]);
        mpz_clear(intermediate_vector[deg]);
    }
    # endif
    for (int i=0; i<N_INT_RED_CLASSICAL; i++)
        reduction_montgomery_mpz(degree, vector, vector, L, L_INV);

    for (int deg=0; deg<degree; deg++)
        out[deg] = mpz_get_si(vector[deg]);
}


void convert_element_to_pmns_fast(int extension_degree, int degree, int64_t out[degree], const mp_limb_t element_data[extension_degree][N_LIMBS]){
    // convert an element of the extension field given to its representation in PMNS using a fast method
    // This is the fastest conversion method but requires an important number of precomputations
    __int128 polynomial[degree];    
    memset(polynomial, 0, sizeof(polynomial));
    for (int deg=0; deg<extension_degree; deg++){
        size_t mask_pos = 0;

        for (int i=0; i<N_POL; i++){
            uint64_t part = apply_mask(element_data[deg], THETA_POW, &mask_pos);
            addmul_pol64_int64(degree, polynomial, PMNS_THETA_FAST[deg][i], part);
        }
    }
    reduction_montgomery_int128(degree, out, polynomial, L, L_INV);
}


void convert_element_to_pmns_pseudo_fast(int extension_degree, int degree, int64_t out[degree], const mp_limb_t element_data[extension_degree][N_LIMBS]){
    // convert an element of the extension field given to its representation in PMNS using a pseudo-fast 
    // This method require less precomputation than the fast method and is slower but for large element, it is faster than the exact converison method
    int POL_LIMBS = N_INT_RED_PSEUDO_FAST + 1;
    mp_limb_t vector[degree][POL_LIMBS];
    mp_limb_t partial_polynomial[degree][POL_LIMBS];

    for (int i=0; i<degree; i++) mpn_zero(vector[i], POL_LIMBS);

    // convert each integer in PMNS representation and multiply theme by a representation
    for (int deg=0; deg<extension_degree; deg++){
        size_t mask_pos = 0;

        for (int i=0; i<degree; i++) mpn_zero(partial_polynomial[i], POL_LIMBS);

        for (int i=0; i<N_POL; i++){
            uint64_t part = apply_mask(element_data[deg], THETA_POW, &mask_pos);
            addmul_pol64_int64_mpn(degree, POL_LIMBS, partial_polynomial, PMNS_THETA_PSEUDO_FAST[i], part);
        }

    # if IS_ELEMENTS_IN_GAMMA_BASIS
        addmul_polmpn_Xpow_modE(degree, POL_LIMBS, vector, partial_polynomial, deg);
    # else
        addmul_polmpn_pol64(degree, POL_LIMBS, vector, partial_polynomial, PMNS_FIELD_ROOTS[deg]);
    # endif
    }

    for (int it = 0; it < N_INT_RED_PSEUDO_FAST; it++)
        reduction_montgomery_mpn(POL_LIMBS, degree, vector, vector, L, L_INV);

    for (int i=0; i<degree; i++)
        out[i] = (int64_t)vector[i][0];
}


void convert_pmns_to_element(int extension_degree, int degree, mp_limb_t out[extension_degree][N_LIMBS], int64_t polynomial[degree]) {
    // Convert a polynomial in PMNS representation to its representation in the extension field using 
    // mpz_t to handle the large integers and apply the external reduction if needed easily. 
    // The result is stored in the out array.
    mpz_t prime, tmp_gamma_pow_part;
    mpz_init(prime);
    mpz_import(prime, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, P);

    mpz_t tmp_result[extension_degree];
    for (int j = 0; j < extension_degree; j++) mpz_init(tmp_result[j]);

    mpz_init(tmp_gamma_pow_part);
    for (int i = 0; i < degree; i++) {
        int64_t coeff = polynomial[i];
        if (coeff == 0) continue;

        for (int j = 0; j < extension_degree; j++) {
            mpz_import(tmp_gamma_pow_part, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, GAMMA_POW[i][j]);
            mpz_mul_si(tmp_gamma_pow_part, tmp_gamma_pow_part, coeff);
            mpz_add(tmp_result[j], tmp_result[j], tmp_gamma_pow_part);
        }
    }

    for (int j = 0; j < extension_degree; j++) {
        mpz_mod(tmp_result[j], tmp_result[j], prime);
        mpn_zero(out[j], N_LIMBS);

        mpn_copyi(out[j], mpz_limbs_read(tmp_result[j]), N_LIMBS);
        mpz_clear(tmp_result[j]);
    }

    mpz_clear(tmp_gamma_pow_part);
    mpz_clear(prime);
}