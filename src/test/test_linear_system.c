#include <obliv_common.h>
#include <errno.h>
#include <obliv.h>
#include <stdio.h>
#include "linear.h"
#include "util.h"
#include "check_error.h"
#include "input.h"
#include "tridiag.h"
int precision = 54;

int read_ls_from_file(int party, const char *filepath, linear_system_t *ls) {
	FILE *file = NULL;
	int res;
	matrix_t A, A_mask;
	vector_t b, b_mask;
	A.value = A_mask.value = b.value = b_mask.value = NULL;

	check(ls && filepath, "Arguments may not be null.");
	file = fopen(filepath, "r");
	//check(file, "Could not open file: %s.", strerror(errno));

	// Read linear system from file
	res = read_matrix(file, &A, precision, false, 0);
	check(!res, "Could not read A.");
	res = read_vector(file, &b, precision, false, 0);
	check(!res, "Could not read b.");

	// generate random masks
	A_mask.value = calloc(A.d[0] * A.d[1] , sizeof(fixed_t));
	b_mask.value = calloc(b.len , sizeof(fixed_t));
	A_mask.d[0] = A.d[0];
	A_mask.d[1] = A.d[1];
	b_mask.len = b.len;

	for(size_t i = 0; i < A.d[0]; i++) {
		for(size_t j = 0; j < A.d[1]; j++) {
			size_t index = i*A.d[1]+j;
			//if(party ==2) printf("%g ", fixed_to_double(A.value[index], precision));
			A_mask.value[index] = 123456; // Of course in a real application random masks would be used
			A.value[index] = (uint64_t) A.value[index] - (uint64_t) A_mask.value[index];
		}
		//if(party==2)printf("\n");
		b_mask.value[i] = 0xDEADBEEF;
		b.value[i] = (uint64_t) b.value[i] - (uint64_t) b_mask.value[i];
	}

	fclose(file);
	file = NULL;

	// Construct instance ls
	ls->precision = precision;
	ls->self = NULL;
	if(party != 1) {
		ls->a = A;
		ls->b = b;
		ls->beta.value = NULL;
		ls->beta.len = -1;
		free(A_mask.value);
		free(b_mask.value);
	} else {
		ls->a = A_mask;
		ls->b = b_mask;
		//ls->beta.value = malloc( A_mask.d[1]*(A_mask.d[1] + 1)/2*sizeof(fixed_t));
		ls->beta.value = malloc(A_mask.d[1]*sizeof(fixed_t));
		//ls->beta.len = A_mask.d[1]*(A_mask.d[1]+1)/2;
		ls->beta.len = A_mask.d[1];
		free(A.value);
		free(b.value);
	}
	return 0;

error:	// for some reason, oblivc removes this label if the stuff
	// below isn't commented out. TODO: fix this and do proper cleanup
	/* fclose(file);
	free(A_mask.value);
	free(b_mask.value);
	free(A.value);
	free(b.value);*/
	return 1;
}


int main(int argc, char **argv) {
	check(argc != 6, "Usage: %s [Port] [Party] [Input file] [Algorithm] [Num. iterations CGD] [Precision]", argv[0]);
	char *algorithm = argv[4];
	check(!strcmp(algorithm, "cholesky") || !strcmp(algorithm, "ldlt")  || !strcmp(algorithm, "cgd"),
	      "Algorithm must be cholesky, ldlt, or cgd.");
	//check(strcmp(algorithm, "cgd") || argc == 6, "Number of iterations for CGD must be provided");
	char *end;
	precision = (int) strtol(argv[6], &end, 10);
	check(!errno, "strtol: %s", strerror(errno));
	check(!*end, "Precision must be a number");
	int party = 0;
	if(!strcmp(argv[2], "1")) {
		party = 1;
	} else if(!strcmp(argv[2], "2")) {
		party = 2;
	}
	check(party > 0, "Party must be either 1 or 2.");

	linear_system_t ls;
	read_ls_from_file(party, argv[3], &ls);
	if(!strcmp(algorithm, "cgd")){
	       ls.num_iterations = atoi(argv[5]);
	} else {
	       ls.num_iterations = 0;
	}

	ProtocolDesc pd;
	ocTestUtilTcpOrDie(&pd, party==1, argv[1]);
	setCurrentParty(&pd, party);

	double time = wallClock();
	if(party == 2) {
	      printf("\n");
	      printf("Algorithm: %s\n", algorithm);
	}
	void (*algorithms[])(void *) = {cholesky, ldlt, cgd};
	int alg_index;
	if(!strcmp(algorithm, "cholesky")) {
              alg_index = 0;
	} else if (!strcmp(algorithm, "ldlt")) {
              alg_index = 1;
	} else {
	      alg_index = 2;
	}


        execYaoProtocol(&pd, tridiag, &ls);
	//execYaoProtocol(&pd, algorithms[alg_index], &ls);
	//execDebugProtocol(&pd, algorithms[alg_index], &ls);

	if(party == 2) {
	  //check(ls.beta.len == d, "Computation error.");
	  printf("Time elapsed: %f\n", wallClock() - time);
	  printf("Number of gates: %lld\n", ls.gates);
	  printf("Vector u: ");
	  for(size_t i = 0; i < ls.beta.len; i++) {
	    printf("%20.15f ", fixed_to_double(ls.beta.value[i], precision));
	  }
	  printf("\n");
          printf("Matrix A: ");
	  for(size_t i = 0; i < ls.beta.len; i++) {
		for (size_t j = 0; j <= i; j++){
	  		printf("%20.15f ", fixed_to_double(ls.a.value[idx(i,j)], precision)); 
		}  
	  }
	  printf("\n");
	}

	/* This code is to run all three algorithms
	void (*algorithms[])(void *) = {cholesky, ldlt, cgd};
	char *algorithm_names[] = {"Cholesky", "LDL^T", "Conjugate Gradient Descent"};
	for(int i = 0; i < 3; i++) {
		double time = wallClock();
		if(party == 2) {
			printf("\n");
			printf("Algorithm: %s\n", algorithm_names[i]);
		}
		execYaoProtocol(&pd, algorithms[i], &ls);

		if(party == 2) {
		  //check(ls.beta.len == d, "Computation error.");
			printf("Time elapsed: %f\n", wallClock() - time);
			printf("Number of gates: %d\n", ls.gates);
			printf("Result: ");
			for(size_t i = 0; i < ls.beta.len; i++) {
				printf("%20.15f ", fixed_to_double(ls.beta.value[i], precision));
			}
			printf("\n");
		}
	}
	*/
	cleanupProtocol(&pd);
	//free(ls.a.value);
	//free(ls.b.value);
	//if (ls.beta.value) free(ls.beta.value);

	return 0;
error:
	return 1;
}
