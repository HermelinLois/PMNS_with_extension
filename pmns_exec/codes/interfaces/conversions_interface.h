#ifndef CONVERSION_API
#define CONVERSION_API

#include "../params/pmns_params.h"

void convert_element_to_pmns_exact(int64_t out[DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]);

void convert_element_to_pmns_pseudo_fast(int64_t out[DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]);

void convert_element_to_pmns_fast(int64_t out[DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]);

void convert_pmns_to_element(mp_limb_t out[EXTENSION_DEGREE][N_LIMBS], int64_t polynomial[DEGREE]);
#endif