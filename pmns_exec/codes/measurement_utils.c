#include "./interfaces/measurement_utils_interface.h"
# include <stdio.h>
# include <gmp.h>
#include <string.h>

uint64_t cpucyclesStart(void){
	unsigned hi, lo;
	__asm__ __volatile__ (	"CPUID\n    "
			"RDTSC\n    "
			"mov %%edx, %0\n    "
			"mov %%eax, %1\n    "
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

uint64_t cpucyclesStop(void){	
	unsigned hi, lo;
	__asm__ __volatile__(	"RDTSCP\n    "
			"mov %%edx, %0\n    "
			"mov %%eax, %1\n    "
			"CPUID\n    "
			: "=r" (hi), "=r" (lo)
			:
			: "%rax", "%rbx", "%rcx", "%rdx");
	
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);
}

uint64_t cycles_diff(uint64_t t1, uint64_t t2){
	uint64_t diff_t;
	 if (t2 < t1){
        diff_t = 18446744073709551615ULL - t1;
        return t2 + diff_t + 1;
    }
	return t2 - t1;
}

void quicksort(uint64_t* t, int n){
	if (n > 0){
        int i, temp;
        int j=0;

        for (i=0; i<n-1; i++)
            if (t[i] < t[n-1]){
                temp = t[j];
                t[j] = t[i];
                t[i] = temp;
                j++;
            }
            
        temp = t[j];
        t[j] = t[n-1];
        t[n-1] = temp;
        
        quicksort(t, j);
        quicksort(&t[j+1], n-j-1);
	}
}

uint64_t *quartiles(uint64_t *tab, uint64_t size){
	uint64_t *result;
	uint64_t aux;
	
	result = malloc(3*sizeof(uint64_t));
	quicksort(tab, size);
	aux = size >> 2;
	if (size % 4)
		aux++;
	// Q1
	result[0] = tab[aux-1];

	// Mediane
	result[1]	= tab[(size+1)/2 - 1];

	// Q3
	aux = (3*size) >> 2;

	if ((3*size) % 4) 
		aux++;

	result[2]	= tab[aux - 1];
	return result;
}