#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <gmp.h>

#include "../params/comparaison_params.h"
#include "../codes/interfaces/measurement_utils_interface.h"
#include "../params/reductions_params.h"

#include <flint/fq.h>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>


extern "C" {
void rand_field_element(int extension_degree, mp_limb_t (*a)[N_LIMBS], gmp_randstate_t state);

void convert_element_to_pmns_fast(int extension_degree, int degree, int64_t *out, const mp_limb_t (*element_data)[N_LIMBS]);

void polynomials_product(int degree, __int128 *out, const int64_t *PolA, const int64_t *PolB);

void reduction_montgomery_int128(int degree, int64_t *out, __int128 *polynomial, const int64_t (*sublattice)[], const int64_t (*sublattice_inv)[]);

# if IS_TOEPLITZ_USABLE
void reduction_montgomery_toeplitz(int degree, int64_t *out, __int128 *polynomial, const int64_t *sublattice, const uint64_t *sublattice_inv);

void reduction_montgomery_toeplitz_recursive(int degree, int64_t *out, __int128 *polynomial, const int64_t *sublattice, const uint64_t *sublattice_inv);
#endif

# if IS_DOUBLE_SPARSE
void reduction_montgomery_linear(int extension_degree, int degree, int64_t *out, __int128 *polynomial);
#endif
}

using namespace NTL;
static fq_ctx_t field;
static fmpz_mod_ctx_t ctx;



// =====================================================================
//                        FORMATS & OPERATIONS
// =====================================================================

// ============================ FORMAT =================================
static inline void pmns_format(int64_t poly_res[DEGREE],  mp_limb_t a[EXTENSION_DEGREE][N_LIMBS]){
    // Convert the field element a to its PMNS representation and store it in poly_res
    convert_element_to_pmns_fast(EXTENSION_DEGREE, DEGREE, poly_res, a);
}

static inline void flint_format(fq_t out, mp_limb_t a[EXTENSION_DEGREE][N_LIMBS]){
    // Convert the field element a to its FLINT representation and store it in out
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    fmpz_t coeff;
    mpz_t tmp;
    fmpz_init(coeff);
    mpz_init(tmp);

    for (size_t i = 0; i < EXTENSION_DEGREE; i++) {
        mpz_import(tmp, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, a[i]);
        fmpz_set_mpz(coeff, tmp);

        fmpz_poly_set_coeff_fmpz(poly, i, coeff);
    }

    fq_set_fmpz_poly(out, poly, field);

    fmpz_clear(coeff);
    mpz_clear(tmp);
    fmpz_poly_clear(poly);
}

static inline void ntl_format(ZZ_pE& out, mp_limb_t a[EXTENSION_DEGREE][N_LIMBS]){
    // Convert the field element a to its NTL representation and store it in out
    mpz_t coeff_mpz;
    mpz_init(coeff_mpz);

    ZZ_pX poly;
    for (long i = 0; i < EXTENSION_DEGREE; i++) {
        mpz_import(coeff_mpz, N_LIMBS, -1, sizeof(mp_limb_t), 0, 0, a[i]);
        char* coeff = mpz_get_str(NULL, 10, coeff_mpz);
        ZZ z = conv<ZZ>(coeff);
        free(coeff);

        SetCoeff(poly, i, conv<ZZ_p>(z));
    }

    out = conv<ZZ_pE>(poly);
    mpz_clear(coeff_mpz);
}

// ============================ OPERATIONS ===============================
static inline void pmns_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    /* Define the operation that will be measured for benchmarking. In this case, it's the polynomial multiplication 
     followed by reduction in the PMNS representation.*/
    __int128_t tmp[DEGREE];
    polynomials_product(DEGREE, tmp, poly_a, poly_b);
    reduction_montgomery_int128(DEGREE, poly_res, tmp, L, L_INV);
}

# if IS_TOEPLITZ_USABLE
static inline void pmns_toeplitz_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    /* Define the operation that will be measured for benchmarking. In this case, it's the polynomial multiplication 
     followed by reduction in the PMNS representation.*/
    __int128_t tmp[DEGREE];
    polynomials_product(DEGREE, tmp, poly_a, poly_b);
    reduction_montgomery_toeplitz(DEGREE, poly_res, tmp, TOEPLITZ_MAT_M, TOEPLITZ_MAT_N);
}

static inline void pmns_toeplitz_recursive_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    /* Define the operation that will be measured for benchmarking. In this case, it's the polynomial multiplication 
     followed by reduction in the PMNS representation.*/
    __int128_t tmp[2 * DEGREE - 1];
    polynomials_product(DEGREE, tmp, poly_a, poly_b);
    reduction_montgomery_toeplitz_recursive(DEGREE, poly_res, tmp, TOEPLITZ_MAT_M, TOEPLITZ_MAT_N);
}
#endif

# if IS_DOUBLE_SPARSE
static inline void pmns_linear_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    /* Define the operation that will be measured for benchmarking. In this case, it's the polynomial multiplication 
     followed by reduction in the PMNS representation.*/
    __int128_t tmp[DEGREE];
    polynomials_product(DEGREE, tmp, poly_a, poly_b);
    reduction_montgomery_linear(EXTENSION_DEGREE, DEGREE, poly_res, tmp);
}
#endif

static inline void flint_operation(fq_t out, fq_t a, fq_t b){
    /* Define the operation that will be measured for benchmarking. In this case, it's the multiplication 
    of two field elements */
    fq_mul(out, a, b, field);
}

static inline void ntl_operation(ZZ_pE& res, const ZZ_pE& a, const ZZ_pE& b){
    /* Define the operation that will be measured for benchmarking. In this case, it's the multiplication 
    of two field elements */
    res = a * b;
}


// =====================================================================
//                             BENCHMARKS
// =====================================================================

// ======================= TEST VALUES GENERATION ======================
static inline void gen_tests_pool(mp_limb_t pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], gmp_randstate_t state){
    // Generate random field elements to be used for benchmarking
	for (int i = 0; i < N_BENCH_SAMPLES; i++) {
		rand_field_element(EXTENSION_DEGREE, pool[i], state);
		rand_field_element(EXTENSION_DEGREE, pool[i + N_BENCH_SAMPLES], state);
	}
}

// ======================== TESTS GENERATION ===========================
/*We define a generic benchmarking macro to measure the performance of different operations with the different libraries*/
#define BENCHMARKS_CORE(TEST_POOL, INIT_ELEMENTS, CLEAR_ELEMENTS, METHOD_NAME, FORMAT_FUNC, OPERATION_FUNC){		\
	uint64_t *cycles = (uint64_t *)calloc(N_BENCH_TESTS,sizeof(uint64_t)), *statTimer;	\
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,meanTimermax = 0;	\
    uint64_t t1,t2, diff_t;																\
																						\
	INIT_ELEMENTS																		\
																						\
	for(int i=0;i<N_BENCH_TESTS; i++){													\
        FORMAT_FUNC(obj_a, TEST_POOL[i]);												\
        FORMAT_FUNC(obj_b, TEST_POOL[i + N_BENCH_SAMPLES]);							\
																						\
        OPERATION_FUNC(obj_res, obj_a, obj_b);										\
	}																					\
																						\
	for(int i=0; i<N_BENCH_SAMPLES; i++){												\
		FORMAT_FUNC(obj_a, TEST_POOL[i]);												\
        FORMAT_FUNC(obj_b, TEST_POOL[i + N_BENCH_SAMPLES]);							\
																						\
		timermin = (uint64_t)0x1<<63;													\
		timermax = 0;																	\
																						\
		memset(cycles,0,N_BENCH_TESTS*sizeof(uint64_t));								\
																						\
		for(int j=0;j<N_BENCH_TESTS;j++){												\
			t1 = cpucyclesStart();														\
			OPERATION_FUNC(obj_res, obj_a, obj_b);										\
			t2 = cpucyclesStop();														\
																						\
			diff_t = cycles_diff(t1, t2);												\
			if(timermin > diff_t) timermin = diff_t;									\
			if(timermax < diff_t) timermax = diff_t;									\
			cycles[j]=diff_t;															\
		}																				\
		meanTimermin += timermin;														\
		meanTimermax += timermax;														\
																						\
		statTimer = quartiles(cycles,N_BENCH_TESTS);									\
																						\
		medianTimer += statTimer[1];													\
		free(statTimer);																\
	}																					\
	CLEAR_ELEMENTS																		\
	free(cycles);																		\
																						\
    printf("====================================================\n");					\
    printf("|%*s%*s|\n", 25 + (int)(strlen(METHOD_NAME)/2), METHOD_NAME, 25 - (int)(strlen(METHOD_NAME)/2), "");\
    printf("====================================================\n");					\
    printf("| %-20s : %-26llu|\n", "Mean Minimum cycles", (unsigned long long)(meanTimermin/N_BENCH_SAMPLES));\
    printf("| %-20s : %-26llu|\n", "Median cycles", (unsigned long long)(medianTimer/N_BENCH_SAMPLES));\
    printf("| %-20s : %-26llu|\n", "Mean Maximum cycles", (unsigned long long)(meanTimermax/N_BENCH_SAMPLES));\
    printf("====================================================\n");					\
}

// =================== PARAMETERS INITAILISATION AND CLEAR ===========================
/*We define macro to initialize and clear elements of each approach*/
#define INIT_PMNS_ELEMENTS int64_t obj_a[DEGREE], obj_b[DEGREE], obj_res[DEGREE];
#define INIT_NTL_ELEMENTS ZZ_pE obj_a, obj_b, obj_res;
#define NO_ELEMENT_TO_CLEAR

#define INIT_FLINT_ELEMENTS 	\
    fq_t obj_a, obj_b, obj_res; \
    fq_init(obj_a, field); 		\
    fq_init(obj_b, field); 		\
    fq_init(obj_res, field);

#define CLEAR_FLINT_ELEMENTS 	\
    fq_clear(obj_a, field); 	\
    fq_clear(obj_b, field); 	\
    fq_clear(obj_res, field);


// ============================ BENCHMARKS DEFINITION =================================
void do_pmns_generic_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
    /* Do the bench march for the product between the vector and a generic matrix*/
	BENCHMARKS_CORE(tests_pool, INIT_PMNS_ELEMENTS, NO_ELEMENT_TO_CLEAR, method_name, pmns_format, pmns_operation);
}

# if IS_TOEPLITZ_USABLE
void do_pmns_toeplitz_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
    /* Do the bench march for the product using Toeplitz matrix structure */
	BENCHMARKS_CORE(tests_pool, INIT_PMNS_ELEMENTS, NO_ELEMENT_TO_CLEAR, method_name, pmns_format, pmns_toeplitz_operation);
}

void do_pmns_toeplitz_recursive_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
    /* Do the bench march for the product using Toeplitz matrix structure and recursive algorithm */
	BENCHMARKS_CORE(tests_pool, INIT_PMNS_ELEMENTS, NO_ELEMENT_TO_CLEAR, method_name, pmns_format, pmns_toeplitz_recursive_operation);
}
#endif

# if IS_DOUBLE_SPARSE
void do_pmns_linear_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
    /* Do the bench march for the product using linear matrix structure */
	BENCHMARKS_CORE(tests_pool, INIT_PMNS_ELEMENTS, NO_ELEMENT_TO_CLEAR, method_name, pmns_format, pmns_linear_operation);
}
#endif

void do_ntl_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
    /* Do the bench march for the product using ntl library */
    BENCHMARKS_CORE(tests_pool, INIT_NTL_ELEMENTS, NO_ELEMENT_TO_CLEAR, method_name, ntl_format, ntl_operation);
}

void do_flint_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
    /* Do the bench march for the product using flint library */
	BENCHMARKS_CORE(tests_pool, INIT_FLINT_ELEMENTS, CLEAR_FLINT_ELEMENTS, method_name, flint_format, flint_operation);
}

// =====================================================================
//                             MAIN
// =====================================================================

void test_libraries_vs_pmns(){
    /*Compare the performance of different libraries for a specific operation*/
    // generate a random seed and initialize testing pool
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    mp_limb_t pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS];
    gen_tests_pool(pool, state);

    // Run benchmarks for each method
    do_pmns_generic_bench(pool, "PMNS (Generic matrix)");

	# if IS_TOEPLITZ_USABLE
    do_pmns_toeplitz_bench(pool, "PMNS (Toeplitz matrix)");
	do_pmns_toeplitz_recursive_bench(pool, "PMNS (Toeplitz recursive matrix)");
	#endif
	
	# if IS_DOUBLE_SPARSE
    do_pmns_linear_bench(pool, "PMNS (Linear representation)");
	#endif

    do_flint_bench(pool, "FLINT");
    do_ntl_bench(pool, "NTL");

    gmp_randclear(state);
}

int main() {
    // construct the fields for NTL library and FLINT library
    construct_ntl_field();
    construct_flint_field(field, ctx);
    test_libraries_vs_pmns();
    free_flint_field(field, ctx);
    return 0;

}