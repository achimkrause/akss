#ifndef ABELIAN_H
#define ABELIAN_H

#include "matrix.h"
#include "smith.h"
#include <gmp.h>
#include <stdio.h>

////////////////////////structs///////////////////////////

struct abelian {
	int *orders; //this is an array of length tor_rank containing the nonzero orders.
	int tor_rank;
	int free_rank;
};



///////////////////////homological algebra////////////////////////////////////////////////////////////
void kernel(int p, struct matrix *f, struct abelian x, struct abelian y, struct matrix_arr to_X, struct matrix_arr from_X, struct abelian *k, struct matrix_arr *to_K, struct matrix_arr *from_K);
void cokernel(int p, struct matrix *f, struct abelian y, struct matrix_arr to_Y, struct matrix_arr from_Y, struct abelian *c, struct matrix_arr *to_C, struct matrix_arr *from_C);
void epi_mono(int p, struct matrix *f, struct abelian x, struct abelian y, struct abelian *img, struct matrix **proj, struct matrix **inc);

/////////////////////diagonal helper functions///////////////////////////////////////////
void compose_diag_p_power(int p, struct matrix *f, int *exponents, struct matrix **res);
void lift_diag_p_power(int p, struct matrix *f, int *exponents, struct matrix **res);

///////////////////////matrix helper functions, maybe factor into matrix.h////////////////////////
void compose(struct matrix *g, struct matrix *f, struct matrix **gf);


///////////////////////lifting, not used as of now/////////////////////////////////////////////
int lift_diag(struct matrix *f, struct matrix *d, struct matrix **res);



void abelian_init(struct abelian *x, int tor_rank, int free_rank);
void abelian_clear(struct abelian x);

//int lift(int p, struct matrix *f, struct matrix *g, struct matrix **res);



#endif