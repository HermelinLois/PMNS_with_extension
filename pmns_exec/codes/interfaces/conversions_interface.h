#ifndef CONVERSION_API
#define CONVERSION_API

#include "../params/pmns_params.h"

void convert_element_to_pmns_exact(int extension_degree, int degree, int64_t out[degree], const mp_limb_t element_data[extension_degree][N_LIMBS]);

void convert_element_to_pmns_pseudo_fast(int extension_degree, int degree, int64_t out[degree], const mp_limb_t element_data[extension_degree][N_LIMBS]);

void convert_element_to_pmns_fast(int extension_degree, int degree, int64_t out[degree], const mp_limb_t element_data[extension_degree][N_LIMBS]);

void convert_pmns_to_element(int extension_degree, int degree, mp_limb_t out[extension_degree][N_LIMBS], int64_t polynomial[degree]);
#endif