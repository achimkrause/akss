#ifndef SMITH_H
#define SMITH_H

#include <gmp.h>
#include "matrix.h"


//////////////////////////////////valuation helper functions//////////////////////////////////////////////////////
int valuation(const int p, const mpq_t lam);
int valuation_int(const int p, const mpz_t z);
/////////////////////////////////////p-local smith////////////////////////////////////////////////////////////////////////
void smith(int p, struct matrix *mat, struct matrix_arr to_X, struct matrix_arr from_X_0,struct matrix_arr to_Y_0, struct matrix_arr from_Y);
#endif
