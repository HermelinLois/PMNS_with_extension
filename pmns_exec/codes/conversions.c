# include "../codes/interfaces/conversions_interface.h"
# include "../codes/interfaces/reductions_interface.h"
# include "../codes/interfaces/operations_interface.h"
# include "../params/conversions_params.h"
# include "../params/reductions_params.h"
# include <stdio.h>
# include <gmp.h>


// apply a mask of size mask_size on the element data starting from the position given by mask_start_pos and 
// return the result as a uint64_t, then increase the mask_start_pos by mask_size
static inline uint64_t apply_mask(const mp_limb_t element_data[N_LIMBS], uint64_t mask_size, size_t *mask_start_pos) {
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


// convert an element of the extension field given to its representation in PMNS using a fast method
void convert_element_to_pmns_fast(int64_t out[DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]){
    __int128 polynomial[DEGREE] = {0};    

    for (int deg=0; deg<EXTENSION_DEGREE; deg++){
        size_t mask_pos = 0;

        for (int i=0; i<N_POL; i++){
            uint64_t part = apply_mask(element_data[deg], THETA_POW, &mask_pos);
            addmul_pol64_int64(polynomial, PMNS_THETA_FAST[deg][i], part);
        }
    }
    
    reduction_montgomery_int128(out, polynomial, L, L_INV);
}



// convert an element of the extension field given to its representation in PMNS using a pseudo-fast method
void convert_element_to_pmns_pseudo_fast(int64_t out[DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]){
    int POL_LIMBS = N_INT_RED_PSEUDO_FAST + 1;
    mp_limb_t vector[DEGREE][POL_LIMBS];
    mp_limb_t partial_polynomial[DEGREE][POL_LIMBS];

    for (int i=0; i<DEGREE; i++) mpn_zero(vector[i], POL_LIMBS);

    // convert each integer in PMNS representation and multiply theme by a representation
    for (int deg=0; deg<EXTENSION_DEGREE; deg++){
        size_t mask_pos = 0;

        for (int i=0; i<DEGREE; i++) mpn_zero(partial_polynomial[i], POL_LIMBS);

        for (int i=0; i<N_POL; i++){
            uint64_t part = apply_mask(element_data[deg], THETA_POW, &mask_pos);
            addmul_pol64_int64_mpn(POL_LIMBS, partial_polynomial, PMNS_THETA_PSEUDO_FAST[i], part);
        }

    # if ELEMENTS_EXPRESS_IN_GAMMA_BASIS
        addmul_polmpn_Xpow_modE(POL_LIMBS, vector, partial_polynomial, deg);
    # else
        addmul_polmpn_pol64(POL_LIMBS, vector, partial_polynomial, PMNS_FIELD_ROOTS[deg]);
    # endif
    }

    for (int it = 0; it < N_INT_RED_PSEUDO_FAST; it++)
        reduction_montgomery_mpn(POL_LIMBS, vector, vector, L, L_INV);

    for (int i=0; i<DEGREE; i++)
        out[i] = (int64_t)vector[i][0];
}

// construct a PMNS vector
void convert_element_to_pmns_vector(int64_t out[EXTENSION_DEGREE][DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]){
    for (int deg=0; deg<EXTENSION_DEGREE; deg++){
        __int128 polynomial[DEGREE] = {0};
        size_t mask_pos = 0;

        for (int i=0; i<N_POL; i++){
            uint64_t part = apply_mask(element_data[deg], THETA_POW, &mask_pos);
            addmul_pol64_int64(polynomial, PMNS_THETA_FAST[deg][i], part);
        }
        reduction_montgomery_int128(out[deg], polynomial, L, L_INV);
    }
}

// convert a PMNS représentation to the represented element in the extension field
void convert_pmns_to_element(mp_limb_t out[EXTENSION_DEGREE][N_LIMBS], int64_t polynomial[DEGREE]) {
    mpz_t prime, tmp_gamma_pow_part;
    mpz_init(prime);
    mpz_import(prime, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, P);

    mpz_t tmp_result[EXTENSION_DEGREE];
    for (int j = 0; j < EXTENSION_DEGREE; j++) mpz_init(tmp_result[j]);

    mpz_init(tmp_gamma_pow_part);
    for (int i = 0; i < DEGREE; i++) {
        int64_t coeff = polynomial[i];
        if (coeff == 0) continue;

        for (int j = 0; j < EXTENSION_DEGREE; j++) {
            mpz_import(tmp_gamma_pow_part, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, GAMMA_POW[i][j]);
            mpz_mul_si(tmp_gamma_pow_part, tmp_gamma_pow_part, coeff);
            mpz_add(tmp_result[j], tmp_result[j], tmp_gamma_pow_part);
        }
    }

    for (int j = 0; j < EXTENSION_DEGREE; j++) {
        mpz_mod(tmp_result[j], tmp_result[j], prime);
        mpn_zero(out[j], N_LIMBS);

        mpn_copyi(out[j], mpz_limbs_read(tmp_result[j]), N_LIMBS);
        mpz_clear(tmp_result[j]);
    }

    mpz_clear(tmp_gamma_pow_part);
    mpz_clear(prime);
}

// init mpz_t with element data and multiply the values by phi^power mod p
static inline void init_element_times_phipow(mpz_t out[DEGREE], const mp_limb_t data[EXTENSION_DEGREE][N_LIMBS], unsigned int power) {
    mpz_t prime, coeff, acc;
    mpz_inits(prime, coeff, acc, NULL);
    mpz_import(prime, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, P);
    
    for (int deg = 0; deg < EXTENSION_DEGREE; deg++) {
        mpz_init(out[deg]);
        mpz_import(coeff, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, data[deg]);
        mpz_mul_2exp(coeff, coeff, PHI_POW * power);
        mpz_mod(out[deg], coeff, prime);
    }
    
    for (int j=EXTENSION_DEGREE; j<DEGREE; j++) mpz_init_set_ui(out[j], 0);
    mpz_clears(prime, coeff, acc, NULL);
}


# if !ELEMENTS_EXPRESS_IN_GAMMA_BASIS
// compute the product of the given element with the transition matrix and apply modulus p to the result
// with our implementation, element are express in gamma basis, no need to convert element

static inline void element_change_basis(mpz_t out[DEGREE], mpz_t extension_element[DEGREE]){
    mpz_t acc, coeff, prime, tmp;
    mpz_inits(acc, coeff, prime, tmp, NULL);
    mpz_import(prime, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, P);

    for (int j=0; j<EXTENSION_DEGREE; j++){
        mpz_init_set_ui(tmp, 0);

        for (int i=0; i<EXTENSION_DEGREE; i++){
            mpz_import(coeff, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, TRANSITION_MATRIX[i][j]);
            mpz_mul(coeff, extension_element[i], coeff);
            mpz_add(tmp, tmp, coeff);
        }
        mpz_mod(out[j], tmp, prime);
    }
    mpz_clears(acc, coeff, prime, tmp, NULL);
}
# endif

// convert an element of the extension field given to its representation in PMNS using the classical method
void convert_element_to_pmns_exact(int64_t out[DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]){
    mpz_t vector[DEGREE];
    
    init_element_times_phipow(vector, element_data, N_INT_RED_CLASSICAL);

# if !ELEMENTS_EXPRESS_IN_GAMMA_BASIS
    element_change_basis(vector, vector);
# endif

    for (int i=0; i<N_INT_RED_CLASSICAL; i++)
        reduction_montgomery_mpz(vector, vector, L, L_INV);

    for (int deg=0; deg<DEGREE; deg++)
        out[deg] = mpz_get_si(vector[deg]);
}