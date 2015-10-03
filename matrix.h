#ifndef MATRIX_H
#define MATRIX_H

//#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

////////////////////////structs///////////////////////////////////////////////

typedef struct {
  mpq_t *entries;
  int height;
  int width;
} Matrix;

typedef struct {
  Matrix **entries;
  int length;
} MatrixArray;

///////////matrix: init and clear /////////////////////////////////////////////////////////////

Matrix* matrix_init(int height, int width);
void matrix_clear(Matrix *mat);

////////matrix_arr//////////////////////////////////////////////////////

void matrix_arr_clear(MatrixArray arr);
void matrix_arr_init(MatrixArray *arr, int length);

////////elementary manipulations//////////////////////////////////////////////////////////

void add_row(Matrix *mat, const int i1, const int i2, const mpq_t lam);
void add_col(Matrix *mat, const int j1, const int j2, const mpq_t lam);
void sub_row(Matrix *mat, const int i1, const int i2, const mpq_t lam);
void sub_col(Matrix *mat, const int j1, const int j2, const mpq_t lam);
void swap_row(Matrix *mat, const int i1, const int i2);
void swap_col(Matrix *mat, const int j1, const int j2);
void mul_row(Matrix *mat, const int i, const mpq_t lam);
void mul_col(Matrix *mat, const int j, const mpq_t lam);
void div_row(Matrix *mat, const int i, const mpq_t lam);
void div_col(Matrix *mat, const int j, const mpq_t lam);

//////////basis manipulation functions/////////////////////////

void add_base(MatrixArray to_X, MatrixArray from_X, int i1, int i2,
    mpq_t lam);
void swap_base(MatrixArray to_X, MatrixArray from_X, int i1, int i2);
void mul_base(MatrixArray to_X, MatrixArray from_X, int i,
    mpq_t lam);

///////////matrix construction functions////////////////////////////////////////////////////////////////////

void set_unit(int i0, int j0, int i_range, int j_range, Matrix *mat);
void set_submatrix(int i0_source, int j0_source, int i0_target, int j0_target,
    int i_range, int j_range, Matrix *source, Matrix *target);
void set_diag_p_powers(int p, int i0_target, int j0_target, int range,
    int *source, Matrix *target);

////////////////////////////////////////////////////////////////////////////////////

#endif
