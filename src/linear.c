#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fixed.h"
#include "linear.h"
#include "check_error.h"


size_t idx(size_t i, size_t j) {
	if(j > i) {
		i^=j; j^=i; i^=j;
	}
	return (i * (i + 1)) / 2 + j;
}

fixed_t inner_product(vector_t *vector_1, vector_t *vector_2){
	fixed_t res = 0;
	for(size_t i = 0; i < vector_1->len; i++) {
		res += vector_1->value[i]*vector_2->value[i];
	}
	return res;
}




int read_matrix(FILE *file, matrix_t *matrix, int precision,
		bool normalize, double normalizer) {
	int n, m, res;
	check(matrix && file, "Arguments may not be null.");
	matrix->value = NULL;

	res = fscanf(file, "%d", &n);
	check(res == 1, "fscanf: %s.", strerror(errno));
	res = fscanf(file, "%d", &m);
	check(res == 1, "fscanf: %s.", strerror(errno));

	matrix->d[0] = n;
	matrix->d[1] = m;
	matrix->value = malloc(n*m*sizeof(fixed_t));

	// printf("A = \n");
	for(size_t i = 0; i < n; i++) {
		for(size_t j = 0; j < m; j++) {
			double val;
			res = fscanf(file, "%lf", &val);
			if(normalize){
				val /= normalizer;
			}
			check(res == 1, "fscanf: %s.", strerror(errno));
			matrix->value[i*m+j] = double_to_fixed(val, precision);
			//printf("%3.8f ", val);
		}
		//printf("\n");
	}

	return 0;
	
error:
	if(matrix) {
		matrix->d[0] = matrix->d[1] = 0;
		free(matrix->value);
		matrix->value = NULL;
	}
	return 1;	
}


int read_sym_mat(FILE *file, symmetric_matrix_t *matrix, int precision) {
	int n, res;
	check(matrix && file, "Arguments may not be null.");
	matrix->value = NULL;

	res = fscanf(file, "%d", &n);
	check(res == 1, "fscanf: %s.", strerror(errno));

	matrix->d = n;
	matrix->value = malloc(n*(n+1)/2*sizeof(fixed_t));

	// printf("A = \n");
	for(size_t i = 0; i < n; i++) {
		for(size_t j = 0; j <= i; j++) {
			double val;
			res = fscanf(file, "%lf", &val);
			check(res == 1, "fscanf: %s.", strerror(errno));
			matrix->value[idx(i,j)] = double_to_fixed(val, precision);
			//printf("%3.8f ", val);
		}
		//printf("\n");
	}
	return 0;
	
error:
	if(matrix) {
		matrix->d = 0;
		free(matrix->value);
		matrix->value = NULL;
	}
	return 1;	
}


int read_vector(FILE *file, vector_t *vector,
		int precision, bool normalize, double normalizer) {
	int l, res;
	check(vector && file, "Arguments may not be null.");

	res = fscanf(file, "%d", &l);
	check(res == 1, "fscanf: %s.", strerror(errno));

	vector->len = l;
	vector->value = malloc(l * sizeof(fixed_t));

	//printf("l = %d, b = \n", l);
	for(size_t i = 0; i < l; i++) {
		double val;
		res = fscanf(file, "%lf", &val);
		if (normalize){
			val /= normalizer;
		}
		check(res == 1, "fscanf: %s.", strerror(errno));
		vector->value[i] = double_to_fixed(val, precision);
		// printf("%3.8f ", val);
	}
	//printf("\n");

	return 0;

error:
	if(vector) {
		vector->len = 0;
		free(vector->value);
		vector->value = NULL;
	}
	return 1;

}
