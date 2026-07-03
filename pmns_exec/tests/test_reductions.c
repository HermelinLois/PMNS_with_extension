# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "../codes/interfaces/test_utils_interface.h"
# include "../codes/interfaces/operations_interface.h"
# include "../codes/interfaces/reductions_interface.h"
# include "reductions_values.h"


void test_equality(){
    /* Test the different reduction methods and compare the result
    to the expected values given by python*/
    if (N_TESTS == 0) {
        printf("No tests to run. To add tests, run 'make new-tests NTESTS=<number_of_tests>'\n");
        exit(1);
    }

     __int128 polynomial[DEGREE];
    int64_t out[DEGREE];

    for (int idx=0; idx<N_TESTS; idx++){
        for (int i=0; i<DEGREE; i++) out[i] = 0;
        // compute the polynomial product without internal reduction
        polynomials_product(polynomial, POL_A[idx], POL_B[idx]);

        // compute different reductions and compare the result to the expected values
        reduction_montgomery_lattice(DEGREE, out, polynomial, L, L_INV);
        check_equality(DEGREE, out, MONTGOMERY_PROD_RED[idx], "Montgomery");

        # if IS_DOUBLE_SPARSE
        reduction_montgomery_linear(EXTENSION_DEGREE, DEGREE, out, polynomial);
        check_equality(DEGREE, out, MONTGOMERY_PROD_RED[idx], "Montgomery linear");
        #endif

        # if IS_BABAI_USABLE
        reduction_babai(DEGREE, out, polynomial, L, L_INV_BABAI);
        check_equality(DEGREE, out, BABAI_PROD_RED[idx], "Babai");
        #endif

        # if IS_TOEPLITZ_USABLE
        reduction_montgomery_toeplitz(DEGREE, out, polynomial, TOEPLITZ_MAT_M, TOEPLITZ_MAT_N);
        check_equality(DEGREE, out, MONTGOMERY_PROD_RED_TOEPLITZ[idx], "Montgomery Toeplitz");

        reduction_montgomery_toeplitz_recursive(DEGREE, out, polynomial, TOEPLITZ_MAT_M, TOEPLITZ_MAT_N);
        check_equality(DEGREE, out, MONTGOMERY_PROD_RED_TOEPLITZ[idx], "Montgomery Toeplitz recursive");
        #endif
    }

# if IS_BABAI_USABLE
    printf("Montgomery and Babai reductions seems to work with given parameters\n");
# else 
    printf("Montgomery reductions seems to work with given parameters\n");
# endif
}


int main(){
    test_equality();
    return 0;
}