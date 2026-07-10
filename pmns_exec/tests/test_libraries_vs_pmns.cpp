#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <gmp.h>

#include "../params/comparaison_params.h"
#include "../codes/interfaces/measurement_utils_interface.h"
#include "../params/reductions_params.h"

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fq.h>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>


extern "C" {
void rand_field_element(mp_limb_t (*a)[N_LIMBS], gmp_randstate_t state);

void convert_element_to_pmns_fast(int64_t *out, const mp_limb_t (*element_data)[N_LIMBS]);

void polynomials_product(__int128 *out, const int64_t *PolA, const int64_t *PolB);

void reduction_montgomery_lattice(int64_t *out, __int128 *polynomial, const int64_t (*sublattice)[], const int64_t (*sublattice_inv)[]);

# if TOEPLITZ_IS_USABLE
void reduction_montgomery_toeplitz(int64_t *out, __int128 *polynomial, const int64_t *sublattice, const uint64_t *sublattice_inv);

void reduction_montgomery_toeplitz_recursive(int64_t *out, __int128 *polynomial, const int64_t *sublattice, const uint64_t *sublattice_inv);
#endif

# if LATTICE_IS_DOUBLE_SPARSE
void reduction_montgomery_linear(int64_t *out, __int128 *polynomial);
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
    convert_element_to_pmns_fast(poly_res, a);
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
static inline void pmns_lattice_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    /* Define the operation that will be measured for benchmarking. In this case, it's the polynomial multiplication 
     followed by reduction in the PMNS representation.*/
    __int128 tmp[DEGREE];
    polynomials_product(tmp, poly_a, poly_b);
    reduction_montgomery_lattice(poly_res, tmp, L, L_INV);
}

# if TOEPLITZ_IS_USABLE
static inline void pmns_toeplitz_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    /* Define the operation that will be measured for benchmarking. In this case, it's the polynomial multiplication 
     followed by reduction in the PMNS representation.*/
    __int128 tmp[DEGREE];
    polynomials_product(tmp, poly_a, poly_b);
    reduction_montgomery_toeplitz(poly_res, tmp, TOEPLITZ_MAT_M, TOEPLITZ_MAT_N);
}

static inline void pmns_toeplitz_recursive_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    /* Define the operation that will be measured for benchmarking. In this case, it's the polynomial multiplication 
     followed by reduction in the PMNS representation.*/
    __int128 tmp[DEGREE];
    polynomials_product(tmp, poly_a, poly_b);
    reduction_montgomery_toeplitz_recursive(poly_res, tmp, TOEPLITZ_MAT_M, TOEPLITZ_MAT_N);
}
#endif

# if LATTICE_IS_DOUBLE_SPARSE
static inline void pmns_linear_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    /* Define the operation that will be measured for benchmarking. In this case, it's the polynomial multiplication 
     followed by reduction in the PMNS representation.*/
    __int128 tmp[DEGREE];
    polynomials_product(tmp, poly_a, poly_b);
    reduction_montgomery_linear(poly_res, tmp);
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
		rand_field_element(pool[i], state);
		rand_field_element(pool[i + N_BENCH_SAMPLES], state);
	}
}

// ======================== TESTS GENERATION ===========================
void do_pmns_bench(void (*operation)(int64_t out[DEGREE], int64_t pol_a[DEGREE], int64_t pol_b[DEGREE]), mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
    /* Do the bench march for pmns representation*/
	uint64_t *cycles = (uint64_t *)calloc(N_BENCH_TESTS,sizeof(uint64_t)), *statTimer;	
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,meanTimermax = 0;	
    uint64_t t1,t2, diff_t;																
																						
	int64_t pol_a[DEGREE], pol_b[DEGREE];
    int64_t pol_res[DEGREE];								
																						
	for(int i=0;i<N_BENCH_TESTS; i++){													
        pmns_format(pol_a, tests_pool[i]);												
        pmns_format(pol_b, tests_pool[i + N_BENCH_SAMPLES]);							
																						
        operation(pol_res, pol_a, pol_b);										
	}																					
																						
	for(int i=0; i<N_BENCH_SAMPLES; i++){												
		pmns_format(pol_a, tests_pool[i]);												
        pmns_format(pol_b, tests_pool[i + N_BENCH_SAMPLES]);							
																						
		timermin = (uint64_t)0x1<<63;													
		timermax = 0;																	
																						
		memset(cycles,0,N_BENCH_TESTS*sizeof(uint64_t));								
																						
		for(int j=0;j<N_BENCH_TESTS;j++){												
			t1 = cpucyclesStart();														
			operation(pol_res, pol_a, pol_b);										
			t2 = cpucyclesStop();													
																						
			diff_t = cycles_diff(t1, t2);												
			if(timermin > diff_t) timermin = diff_t;									
			if(timermax < diff_t) timermax = diff_t;									
			cycles[j]=diff_t;															
		}																				
		meanTimermin += timermin;														
		meanTimermax += timermax;														
																						
		statTimer = quartiles(cycles,N_BENCH_TESTS);									
																						
		medianTimer += statTimer[1];													
		free(statTimer);																
	}																					
	free(cycles);																		
																						
    printf("====================================================\n");					
    printf("|%*s%*s|\n", 25 + (int)(strlen(method_name)/2), method_name, 25 - (int)(strlen(method_name)/2), "");
    printf("====================================================\n");					
    printf("| %-20s : %-26llu|\n", "Mean Minimum cycles", (unsigned long long)(meanTimermin/N_BENCH_SAMPLES));
    printf("| %-20s : %-26llu|\n", "Median cycles", (unsigned long long)(medianTimer/N_BENCH_SAMPLES));
    printf("| %-20s : %-26llu|\n", "Mean Maximum cycles", (unsigned long long)(meanTimermax/N_BENCH_SAMPLES));
    printf("====================================================\n");					
}



void do_flint_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
    /* Do the bench march for flint representation*/
	uint64_t *cycles = (uint64_t *)calloc(N_BENCH_TESTS,sizeof(uint64_t)), *statTimer;	
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,meanTimermax = 0;	
    uint64_t t1,t2, diff_t;																
								
    
	fq_t elmnt_a, elmnt_b, elmnt_res;
    fq_init(elmnt_a, field), fq_init(elmnt_b, field), fq_init(elmnt_res, field);							
																						
	for(int i=0;i<N_BENCH_TESTS; i++){													
        flint_format(elmnt_a, tests_pool[i]);												
        flint_format(elmnt_b, tests_pool[i + N_BENCH_SAMPLES]);							
																						
        flint_operation(elmnt_res, elmnt_a, elmnt_b);										
	}																					
																						
	for(int i=0; i<N_BENCH_SAMPLES; i++){												
		flint_format(elmnt_a, tests_pool[i]);												
        flint_format(elmnt_b, tests_pool[i + N_BENCH_SAMPLES]);							
																						
		timermin = (uint64_t)0x1<<63;													
		timermax = 0;																	
																						
		memset(cycles,0,N_BENCH_TESTS*sizeof(uint64_t));								
																						
		for(int j=0;j<N_BENCH_TESTS;j++){												
			t1 = cpucyclesStart();														
			flint_operation(elmnt_res, elmnt_a, elmnt_b);										
			t2 = cpucyclesStop();													
																						
			diff_t = cycles_diff(t1, t2);												
			if(timermin > diff_t) timermin = diff_t;									
			if(timermax < diff_t) timermax = diff_t;									
			cycles[j]=diff_t;															
		}																				
		meanTimermin += timermin;														
		meanTimermax += timermax;														
																						
		statTimer = quartiles(cycles,N_BENCH_TESTS);									
																						
		medianTimer += statTimer[1];													
		free(statTimer);																
	}																					
	free(cycles);	
    fq_clear(elmnt_a, field), fq_clear(elmnt_b, field), fq_clear(elmnt_res, field);																	
																						
    printf("====================================================\n");					
    printf("|%*s%*s|\n", 25 + (int)(strlen(method_name)/2), method_name, 25 - (int)(strlen(method_name)/2), "");
    printf("====================================================\n");					
    printf("| %-20s : %-26llu|\n", "Mean Minimum cycles", (unsigned long long)(meanTimermin/N_BENCH_SAMPLES));
    printf("| %-20s : %-26llu|\n", "Median cycles", (unsigned long long)(medianTimer/N_BENCH_SAMPLES));
    printf("| %-20s : %-26llu|\n", "Mean Maximum cycles", (unsigned long long)(meanTimermax/N_BENCH_SAMPLES));
    printf("====================================================\n");					
}



void do_ntl_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
    /* Do the bench march for ntl representation*/
	uint64_t *cycles = (uint64_t *)calloc(N_BENCH_TESTS,sizeof(uint64_t)), *statTimer;	
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,meanTimermax = 0;	
    uint64_t t1,t2, diff_t;																
								
    
	ZZ_pE elmnt_a, elmnt_b, elmnt_res;						
																						
	for(int i=0;i<N_BENCH_TESTS; i++){													
        ntl_format(elmnt_a, tests_pool[i]);												
        ntl_format(elmnt_b, tests_pool[i + N_BENCH_SAMPLES]);							
																						
        ntl_operation(elmnt_res, elmnt_a, elmnt_b);										
	}																					
																						
	for(int i=0; i<N_BENCH_SAMPLES; i++){												
		ntl_format(elmnt_a, tests_pool[i]);												
        ntl_format(elmnt_b, tests_pool[i + N_BENCH_SAMPLES]);							
																						
		timermin = (uint64_t)0x1<<63;													
		timermax = 0;																	
																						
		memset(cycles,0,N_BENCH_TESTS*sizeof(uint64_t));								
																						
		for(int j=0;j<N_BENCH_TESTS;j++){												
			t1 = cpucyclesStart();														
			ntl_operation(elmnt_res, elmnt_a, elmnt_b);										
			t2 = cpucyclesStop();													
																						
			diff_t = cycles_diff(t1, t2);												
			if(timermin > diff_t) timermin = diff_t;									
			if(timermax < diff_t) timermax = diff_t;									
			cycles[j]=diff_t;															
		}																				
		meanTimermin += timermin;														
		meanTimermax += timermax;														
																						
		statTimer = quartiles(cycles,N_BENCH_TESTS);									
																						
		medianTimer += statTimer[1];													
		free(statTimer);																
	}																					
	free(cycles);																
																						
    printf("====================================================\n");					
    printf("|%*s%*s|\n", 25 + (int)(strlen(method_name)/2), method_name, 25 - (int)(strlen(method_name)/2), "");
    printf("====================================================\n");					
    printf("| %-20s : %-26llu|\n", "Mean Minimum cycles", (unsigned long long)(meanTimermin/N_BENCH_SAMPLES));
    printf("| %-20s : %-26llu|\n", "Median cycles", (unsigned long long)(medianTimer/N_BENCH_SAMPLES));
    printf("| %-20s : %-26llu|\n", "Mean Maximum cycles", (unsigned long long)(meanTimermax/N_BENCH_SAMPLES));
    printf("====================================================\n");					
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
    do_pmns_bench(pmns_lattice_operation, pool, "PMNS (Generic matrix)");

	# if TOEPLITZ_IS_USABLE
    do_pmns_bench(pmns_toeplitz_operation, pool, "PMNS (Toeplitz matrix)");

    do_pmns_bench(pmns_toeplitz_recursive_operation, pool, "PMNS (Toeplitz recursive matrix)");
	#endif
	
	# if LATTICE_IS_DOUBLE_SPARSE
    do_pmns_bench(pmns_linear_operation, pool, "PMNS (Linear representation)");
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