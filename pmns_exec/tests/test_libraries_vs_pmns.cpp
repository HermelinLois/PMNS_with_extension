#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <gmp.h>

#include "../params/comparaison_params.h"
#include "../codes/interfaces/measurement_utils_interface.h"
# include "../params/reductions_params.h"

extern "C" {
void rand_field_element(mp_limb_t a[EXTENSION_DEGREE][N_LIMBS], gmp_randstate_t state);

void convert_element_to_pmns_fast(int64_t out[DEGREE], const mp_limb_t element_data[EXTENSION_DEGREE][N_LIMBS]);

void polynomials_product(__int128 out[DEGREE], int64_t PolA[DEGREE], int64_t PolB[DEGREE]);

void reduction_montgomery_int128(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[DEGREE][DEGREE], const int64_t sublattice_inv[DEGREE][DEGREE]);

void reduction_montgomery_toeplitz(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[2*DEGREE - 1], const uint64_t sublattice_inv[2*DEGREE - 1]);

void reduction_montgomery_toeplitz_recursive(int64_t out[DEGREE], __int128 polynomial[DEGREE], const int64_t sublattice[2*DEGREE - 1], const uint64_t sublattice_inv[2*DEGREE - 1]);

#ifdef IS_SPARSE
void reduction_montgomery_linear(int64_t out[DEGREE], __int128 polynomial[DEGREE]);
#endif
}

#include <flint/fq.h>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>

using namespace NTL;
static fq_ctx_t field;
static fmpz_mod_ctx_t ctx;



// =====================================================================
//                        FORMATS & OPERATIONS
// =====================================================================
static inline void pmns_format(int64_t poly_res[DEGREE], mp_limb_t a[EXTENSION_DEGREE][N_LIMBS]){
    convert_element_to_pmns_fast(poly_res, a);
}

static inline void pmns_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    __int128_t tmp[DEGREE];
    polynomials_product(tmp, poly_a, poly_b);
    reduction_montgomery_int128(poly_res, tmp, L, L_INV);
}

static inline void pmns_toeplitz_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    __int128_t tmp[DEGREE];
    polynomials_product(tmp, poly_a, poly_b);
    reduction_montgomery_toeplitz(poly_res, tmp, TOEPLITZ_MAT_M, TOEPLITZ_MAT_N);
}

static inline void pmns_toeplitz_recursive_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    __int128_t tmp[DEGREE];
    polynomials_product(tmp, poly_a, poly_b);
    reduction_montgomery_toeplitz_recursive(poly_res, tmp, TOEPLITZ_MAT_M, TOEPLITZ_MAT_N);
}

#ifdef IS_SPARSE
static inline void pmns_linear_operation(int64_t poly_res[DEGREE], int64_t poly_a[DEGREE], int64_t poly_b[DEGREE]){
    __int128_t tmp[DEGREE];
    polynomials_product(tmp, poly_a, poly_b);
    reduction_montgomery_linear(poly_res, tmp);
}
#endif

static inline void flint_format(fq_t out, mp_limb_t a[EXTENSION_DEGREE][N_LIMBS]){
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

static inline void flint_operation(fq_t out, fq_t a, fq_t b){
    fq_mul(out, a, b, field);
}


static inline void ntl_format(ZZ_pE& out, mp_limb_t a[EXTENSION_DEGREE][N_LIMBS]){
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

static inline void ntl_operation(ZZ_pE& res, const ZZ_pE& a, const ZZ_pE& b){
    res = a * b;
}




// =====================================================================
//                             BENCHMARKS
// =====================================================================

static inline void gen_tests_pool(mp_limb_t pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], gmp_randstate_t state) {
	for (int i = 0; i < N_BENCH_SAMPLES; i++) {
		rand_field_element(pool[i], state);
		rand_field_element(pool[i + N_BENCH_SAMPLES], state);
	}
}

void do_ntl_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
    uint64_t *cycles = (uint64_t*)calloc(N_BENCH_TESTS, sizeof(uint64_t)), *statTimer;
	uint64_t timermin, timermax, meanTimermin = 0, medianTimer = 0, meanTimermax = 0;
	uint64_t t1, t2, diff_t;

    ZZ_pE fa, fb, fres;

	for (int i = 0; i < N_BENCH_TESTS; i++) {
        mp_limb_t (*a)[N_LIMBS] = tests_pool[i];
		mp_limb_t (*b)[N_LIMBS] = tests_pool[i + N_BENCH_SAMPLES];

        ntl_format(fa, a);
        ntl_format(fb, b);

        ntl_operation(fres, fa, fb);
    }

	for (int i = 0; i < N_BENCH_SAMPLES; i++) {
		mp_limb_t (*a)[N_LIMBS] = tests_pool[i];
		mp_limb_t (*b)[N_LIMBS] = tests_pool[i + N_BENCH_SAMPLES];

		ntl_format(fa, a);
		ntl_format(fb, b);

		timermin = (uint64_t)1ULL << 63;
		timermax = 0;

		memset(cycles, 0, N_BENCH_TESTS * sizeof(uint64_t));

		for (int j = 0; j < N_BENCH_TESTS; j++) {
			t1 = cpucyclesStart();
			ntl_operation(fres, fa, fb);
			t2 = cpucyclesStop();

			diff_t = cycles_diff(t1, t2);
			if (timermin > diff_t) timermin = diff_t;
			if (timermax < diff_t) timermax = diff_t;
			cycles[j] = diff_t;
		}

		meanTimermin += timermin;
		meanTimermax += timermax;
		statTimer = quartiles(cycles, N_BENCH_TESTS);

		medianTimer += statTimer[1];
		free(statTimer);
	}

	free(cycles);

	printf("====================================================\n");
	printf("|%*s%*s|\n", 25 + (int)(strlen(method_name) / 2), method_name, 25 - (int)(strlen(method_name) / 2), "");
	printf("====================================================\n");
	printf("| %-20s : %-26llu|\n", "Mean Minimum cycles", (unsigned long long)(meanTimermin / N_BENCH_SAMPLES));
	printf("| %-20s : %-26llu|\n", "Median cycles", (unsigned long long)(medianTimer / N_BENCH_SAMPLES));
	printf("| %-20s : %-26llu|\n", "Mean Maximum cycles", (unsigned long long)(meanTimermax / N_BENCH_SAMPLES));
	printf("====================================================\n");
}


void do_flint_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
	uint64_t *cycles = (uint64_t *)calloc(N_BENCH_TESTS,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,meanTimermax = 0;
    uint64_t t1,t2, diff_t;

    fq_t fa, fb, fres;

    fq_init(fa, field);
    fq_init(fb, field);
    fq_init(fres, field);

    for(int i=0;i<N_BENCH_TESTS;i++){
		mp_limb_t (*a)[N_LIMBS] = tests_pool[i];
        mp_limb_t (*b)[N_LIMBS] = tests_pool[i + N_BENCH_SAMPLES];

        flint_format(fa, a);
        flint_format(fb, b);
        
        flint_operation(fres, fa, fb);
	}

    fq_clear(fa, field);
    fq_clear(fb, field);
    fq_clear(fres, field);


	for(int i=0;i<N_BENCH_SAMPLES;i++){
		mp_limb_t (*a)[N_LIMBS] = tests_pool[i];
        mp_limb_t (*b)[N_LIMBS] = tests_pool[i + N_BENCH_SAMPLES];

        fq_init(fa, field);
        fq_init(fb, field);
        fq_init(fres, field);

        flint_format(fa, a);
        flint_format(fb, b);

		timermin = (uint64_t)0x1<<63;
		timermax = 0;
        
		memset(cycles, 0, N_BENCH_TESTS*sizeof(uint64_t));

		for(int j=0;j<N_BENCH_TESTS;j++){
			t1 = cpucyclesStart();
			flint_operation(fres, fa, fb);
			t2 = cpucyclesStop();

			diff_t = cycles_diff(t1, t2);
			if(timermin > diff_t) timermin = diff_t;
			if(timermax < diff_t) timermax = diff_t;
			cycles[j]=diff_t;
		}
        fq_clear(fa, field);
        fq_clear(fb, field);
        fq_clear(fres, field);
        
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


void do_pmns_generic_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
	uint64_t *cycles = (uint64_t *)calloc(N_BENCH_TESTS,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,meanTimermax = 0;
    uint64_t t1,t2, diff_t;

    int64_t poly_a[DEGREE], poly_b[DEGREE], poly_res[DEGREE];
	
	for(int i=0;i<N_BENCH_TESTS;i++){
		mp_limb_t (*a)[N_LIMBS] = tests_pool[i];
        mp_limb_t (*b)[N_LIMBS] = tests_pool[i + N_BENCH_SAMPLES];

        pmns_format(poly_a, a);
        pmns_format(poly_b, b);
		
        pmns_operation(poly_res, poly_a, poly_b);
	}
	
	for(int i=0;i<N_BENCH_SAMPLES;i++){
		pmns_format(poly_a, tests_pool[i]);
        pmns_format(poly_b, tests_pool[i + N_BENCH_SAMPLES]);

		timermin = (uint64_t)0x1<<63;
		timermax = 0;
        
		memset(cycles,0,N_BENCH_TESTS*sizeof(uint64_t));

		for(int j=0;j<N_BENCH_TESTS;j++){
			t1 = cpucyclesStart();
			pmns_operation(poly_res, poly_a, poly_b);
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


void do_pmns_toeplitz_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
	uint64_t *cycles = (uint64_t *)calloc(N_BENCH_TESTS,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,meanTimermax = 0;
    uint64_t t1,t2, diff_t;

    int64_t poly_a[DEGREE], poly_b[DEGREE], poly_res[DEGREE];
	
	for(int i=0;i<N_BENCH_TESTS;i++){
		mp_limb_t (*a)[N_LIMBS] = tests_pool[i];
        mp_limb_t (*b)[N_LIMBS] = tests_pool[i + N_BENCH_SAMPLES];

        pmns_format(poly_a, a);
        pmns_format(poly_b, b);
		
        pmns_toeplitz_operation(poly_res, poly_a, poly_b);
	}
	
	for(int i=0;i<N_BENCH_SAMPLES;i++){
		pmns_format(poly_a, tests_pool[i]);
        pmns_format(poly_b, tests_pool[i + N_BENCH_SAMPLES]);

		timermin = (uint64_t)0x1<<63;
		timermax = 0;
        
		memset(cycles,0,N_BENCH_TESTS*sizeof(uint64_t));

		for(int j=0;j<N_BENCH_TESTS;j++){
			t1 = cpucyclesStart();
			pmns_toeplitz_operation(poly_res, poly_a, poly_b);
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



void do_pmns_toeplitz_recursive_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
	uint64_t *cycles = (uint64_t *)calloc(N_BENCH_TESTS,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,meanTimermax = 0;
    uint64_t t1,t2, diff_t;

    int64_t poly_a[DEGREE], poly_b[DEGREE], poly_res[DEGREE];
	
	for(int i=0;i<N_BENCH_TESTS;i++){
		mp_limb_t (*a)[N_LIMBS] = tests_pool[i];
        mp_limb_t (*b)[N_LIMBS] = tests_pool[i + N_BENCH_SAMPLES];

        pmns_format(poly_a, a);
        pmns_format(poly_b, b);
		
        pmns_toeplitz_recursive_operation(poly_res, poly_a, poly_b);
	}
	
	for(int i=0;i<N_BENCH_SAMPLES;i++){
		pmns_format(poly_a, tests_pool[i]);
        pmns_format(poly_b, tests_pool[i + N_BENCH_SAMPLES]);

		timermin = (uint64_t)0x1<<63;
		timermax = 0;
        
		memset(cycles,0,N_BENCH_TESTS*sizeof(uint64_t));

		for(int j=0;j<N_BENCH_TESTS;j++){
			t1 = cpucyclesStart();
			pmns_toeplitz_recursive_operation(poly_res, poly_a, poly_b);
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

#ifdef IS_SPARSE
void do_pmns_linear_bench(mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], const char* method_name){
	uint64_t *cycles = (uint64_t *)calloc(N_BENCH_TESTS,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,meanTimermax = 0;
    uint64_t t1,t2, diff_t;

    int64_t poly_a[DEGREE], poly_b[DEGREE], poly_res[DEGREE];
	
	for(int i=0;i<N_BENCH_TESTS;i++){
		mp_limb_t (*a)[N_LIMBS] = tests_pool[i];
        mp_limb_t (*b)[N_LIMBS] = tests_pool[i + N_BENCH_SAMPLES];

        pmns_format(poly_a, a);
        pmns_format(poly_b, b);
		
        pmns_linear_operation(poly_res, poly_a, poly_b);
	}
	
	for(int i=0;i<N_BENCH_SAMPLES;i++){
		pmns_format(poly_a, tests_pool[i]);
        pmns_format(poly_b, tests_pool[i + N_BENCH_SAMPLES]);

		timermin = (uint64_t)0x1<<63;
		timermax = 0;
        
		memset(cycles,0,N_BENCH_TESTS*sizeof(uint64_t));

		for(int j=0;j<N_BENCH_TESTS;j++){
			t1 = cpucyclesStart();
			pmns_linear_operation(poly_res, poly_a, poly_b);
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
#endif

// =====================================================================
//                             MAIN
// =====================================================================

void test_libraries_vs_pmns(){

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    mp_limb_t pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS];
    gen_tests_pool(pool, state);


    do_pmns_generic_bench(pool, "PMNS (Generic matrix)");
    do_pmns_toeplitz_bench(pool, "PMNS (Toeplitz matrix)");
	do_pmns_toeplitz_recursive_bench(pool, "PMNS (Toeplitz recursive matrix)");
	
#ifdef IS_SPARSE
    do_pmns_linear_bench(pool, "PMNS (Linear representation)");
#endif

    do_flint_bench(pool, "FLINT");
    do_ntl_bench(pool, "NTL");

    gmp_randclear(state);
}

int main() {
    construct_ntl_field();
    construct_flint_field(field, ctx);
    test_libraries_vs_pmns();
    free_flint_field(field, ctx);
    return 0;

}