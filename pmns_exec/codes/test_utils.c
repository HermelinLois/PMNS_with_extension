# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <gmp.h>
# include "../params/pmns_params.h"
# include "../codes/interfaces/test_utils_interface.h"
# include "../codes/interfaces/conversions_interface.h"

void rand_field_element(int extension_degree, mp_limb_t out[extension_degree][N_LIMBS], gmp_randstate_t state){
    /* Generate a random field element as extension_degree coefficients.
    Each coefficient is sampled uniformly in [2^(P_SIZE-1), p-1], to ensure
    the product bewteen two elements need a reduction modulus the prime p. */
	mpz_t p_mpz, rand_val, pow2, tmp;
    mpz_init(p_mpz);
    mpz_init(rand_val);
    mpz_init(pow2);
    mpz_init(tmp);

    mpz_setbit(pow2, P_SIZE-1);
	mpz_import(p_mpz, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, P);
    mpz_sub_ui(tmp, p_mpz, 1);
    mpz_sub(tmp, tmp, pow2);

    for (int i=0; i<extension_degree; i++){
        mpz_urandomm(rand_val, state, tmp);
        mpz_add(rand_val, rand_val, pow2);
        mpn_copyi(out[i], mpz_limbs_read(rand_val), N_LIMBS);
    }

    mpz_clear(p_mpz);
    mpz_clear(rand_val);
    mpz_clear(pow2);
    mpz_clear(tmp);
}


void print_pol(int degree, int64_t P[degree]){
    /*Display degree coefficients of the polynomial*/
    for (int idx=0; idx<degree; idx++)
        printf("%lld  ", (long long)P[idx]);
    printf("\n");
}


void check_equality(int degree, int64_t P[degree], int64_t C[degree], const char* method){
    /*Check if the two polynomials P and C are equal by comparing each coefficient. If they are not equal, 
    print an error message and the two polynomials, then exit the program.*/
    for (int i=0; i<degree; i++){
        if(C[i] != P[i]) {
            printf("/!\\ Error has been found /!\\ \n");
            printf("Generate with Python :\n");
            print_pol(degree, C);
            printf("Generate with C (with %s) :\n", method);
            print_pol(degree, P);
            exit(1);
        }
    }
}

void reset_polynomial(int degree, int64_t polynomial[degree]){
    /*Reset the polynomial to zero by setting each coefficient to 0*/
    for (int i=0; i<degree; i++) 
        polynomial[i]=0;
}