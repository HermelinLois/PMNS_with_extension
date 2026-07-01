# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <string.h>

# include "conversions_values.h"
# include "../codes/interfaces/test_utils_interface.h"
# include "../codes/interfaces/conversions_interface.h"
# include "../codes/interfaces/test_utils_interface.h"
# include "../codes/interfaces/measurement_utils_interface.h"



void test_equality(){
    if (N_TESTS == 0) {
        printf("No tests to run. To add tests, run 'make new-tests NTESTS=<number_of_tests>'\n");
        exit(1);
    }

    int64_t polynomial[DEGREE];

    for (int idx=0; idx<N_TESTS; idx++){
        reset_polynomial(DEGREE, polynomial);
        convert_element_to_pmns_exact(EXTENSION_DEGREE, DEGREE, polynomial, EXTENSION_FIELD_ELEMENTS[idx]);
        check_equality(DEGREE, CONVERTED_ELEMENTS_EXACT[idx], polynomial, "exact");

        reset_polynomial(DEGREE, polynomial);
        convert_element_to_pmns_pseudo_fast(EXTENSION_DEGREE, DEGREE, polynomial, EXTENSION_FIELD_ELEMENTS[idx]);
        check_equality(DEGREE, CONVERTED_ELEMENTS_PSEUDO_FAST[idx], polynomial, "pseudo-fast");

        reset_polynomial(DEGREE, polynomial);
        convert_element_to_pmns_fast(EXTENSION_DEGREE, DEGREE, polynomial, EXTENSION_FIELD_ELEMENTS[idx]);
        check_equality(DEGREE, CONVERTED_ELEMENTS_FAST[idx], polynomial, "fast");
    }
    printf("Conversions seems to work with given parameters\n");
}


void do_bench(void (*to_pmns)(int extension_degree, int degree, int64_t pmns[degree], const mp_limb_t element_data[extension_degree][N_LIMBS]), mp_limb_t tests_pool[N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], char* method_name, gmp_randstate_t state){
	uint64_t *cycles = (uint64_t *)calloc(N_BENCH_TESTS,sizeof(uint64_t)), *statTimer;
	uint64_t timermin , timermax, meanTimermin =0,	medianTimer = 0,meanTimermax = 0;
    uint64_t t1,t2, diff_t;

	mp_limb_t a[EXTENSION_DEGREE][N_LIMBS];
    int64_t polynomial[DEGREE];
	
	for(int i=0;i<N_BENCH_TESTS;i++){
		rand_field_element(EXTENSION_DEGREE, a, state);
		to_pmns(EXTENSION_DEGREE, DEGREE, polynomial, a);
	}
	
	for(int i=0;i<N_BENCH_SAMPLES;i++){
		mp_limb_t (*element_data)[N_LIMBS] = tests_pool[i];

		timermin = (uint64_t)0x1<<63;
		timermax = 0;
        
		memset(cycles,0,N_BENCH_TESTS*sizeof(uint64_t));

		for(int j=0;j<N_BENCH_TESTS;j++){
            reset_polynomial(DEGREE, polynomial);

			t1 = cpucyclesStart();
			to_pmns(EXTENSION_DEGREE, DEGREE, polynomial, element_data);
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
    printf("| %-20s : %-26llu|\n", "Minimum cycles", (unsigned long long)(meanTimermin/N_BENCH_SAMPLES));
    printf("| %-20s : %-26llu|\n", "Median cycles", (unsigned long long)(medianTimer/N_BENCH_SAMPLES));
    printf("| %-20s : %-26llu|\n", "Maximum cycles", (unsigned long long)(meanTimermax/N_BENCH_SAMPLES));
    printf("====================================================\n");
}



static inline void gen_tests_pool(mp_limb_t pool[N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS], gmp_randstate_t state){
    for (int i=0; i<N_BENCH_SAMPLES; i++)
        rand_field_element(EXTENSION_DEGREE, pool[i], state);
}

void test_speed(){
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    mp_limb_t tests_pool[2 * N_BENCH_SAMPLES][EXTENSION_DEGREE][N_LIMBS];
    gen_tests_pool(tests_pool, state);

    do_bench(convert_element_to_pmns_exact, tests_pool, "Exact", state);
    do_bench(convert_element_to_pmns_pseudo_fast, tests_pool, "Pseudo-Fast", state);
    do_bench(convert_element_to_pmns_fast, tests_pool, "Fast", state);
    
    gmp_randclear(state);
}



int main(int argc, char *argv[]){
    if (argc > 1 && strcmp(argv[1], "equality") == 0) {
        test_equality();
        return 0;
    }

    if (argc > 1 && strcmp(argv[1], "speed") == 0) {
        test_speed();
        return 0;
    }

    printf("Usable parameters: [equality|speed]\n");
    return 0;
}