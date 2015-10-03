#ifndef ABELIAN_H
#define ABELIAN_H

#include "matrix.h"
#include "smith.h"
#include <gmp.h>
#include <stdio.h>

////////////////////////structs///////////////////////////
 typedef struct {
  int *orders;  //this is an array of length tor_rank containing the nonzero orders.
  int tor_rank;
  int free_rank;
} AbelianGroup;

///////////////////////homological algebra////////////////////////////////////////////////////////////
void kernel(int p, Matrix *f, AbelianGroup x, AbelianGroup y,
    MatrixArray to_X, MatrixArray from_X, AbelianGroup *k,
    MatrixArray *to_K, MatrixArray *from_K);
void cokernel(int p, Matrix *f, AbelianGroup y, MatrixArray to_Y,
    MatrixArray from_Y, AbelianGroup *c, MatrixArray *to_C,
    MatrixArray *from_C);
void epi_mono(int p, Matrix *f, AbelianGroup x, AbelianGroup y,
    AbelianGroup *img, Matrix **proj, Matrix **inc);

/////////////////////diagonal helper functions///////////////////////////////////////////
void compose_diag_p_power(int p, Matrix *f, int *exponents,
    Matrix **res);
void lift_diag_p_power(int p, Matrix *f, int *exponents,
    Matrix **res);

///////////////////////matrix helper functions, maybe factor into matrix.h////////////////////////
void compose(Matrix *g, Matrix *f, Matrix **gf);

///////////////////////lifting, not used as of now/////////////////////////////////////////////
int lift_diag(Matrix *f, Matrix *d, Matrix **res);

void abelian_init(AbelianGroup *x, int tor_rank, int free_rank);
void abelian_clear(AbelianGroup x);

//int lift(int p, Matrix *f, Matrix *g, Matrix **res);

#endif
