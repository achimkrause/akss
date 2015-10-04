#ifndef SMITH_H
#define SMITH_H

#include <gmp.h>
#include "matrix.h"

//////////////////////////////////valuation helper functions//////////////////////////////////////////////////////
int valuation(const int p, const mpq_t lam);
int valuation_int(const int p, const mpz_t z);
/////////////////////////////////////p-local smith////////////////////////////////////////////////////////////////////////
void smith(int p, Matrix *mat, MatrixArray to_X,
    MatrixArray from_X_0, MatrixArray to_Y_0,
    MatrixArray from_Y);
#endif
