#ifndef MATRIX_H
#define MATRIX_H

//#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

////////////////////////structs///////////////////////////////////////////////

struct matrix {
	mpq_t *entries;
	int height;
	int width;
};

struct matrix_arr {
	struct matrix **entries;
	int length;
};

///////////matrix: init and clear /////////////////////////////////////////////////////////////

struct matrix* matrix_init(int height, int width);
void matrix_clear(struct matrix *mat);

////////matrix_arr//////////////////////////////////////////////////////

void matrix_arr_clear(struct matrix_arr arr);
void matrix_arr_init(struct matrix_arr *arr, int length);

////////elementary manipulations//////////////////////////////////////////////////////////

void add_row(struct matrix *mat, const int i1, const int i2, const mpq_t lam);
void add_col(struct matrix *mat, const int j1, const int j2, const mpq_t lam);
void sub_row(struct matrix *mat, const int i1, const int i2, const mpq_t lam);
void sub_col(struct matrix *mat, const int j1, const int j2, const mpq_t lam);
void swap_row(struct matrix *mat, const int i1, const int i2);
void swap_col(struct matrix *mat, const int j1, const int j2);
void mul_row(struct matrix *mat, const int i, const mpq_t lam);
void mul_col(struct matrix *mat, const int j, const mpq_t lam);
void div_row(struct matrix *mat, const int i, const mpq_t lam);
void div_col(struct matrix *mat, const int j, const mpq_t lam);

//////////basis manipulation functions/////////////////////////

void add_base(struct matrix_arr to_X, struct matrix_arr from_X, int i1, int i2, mpq_t lam);
void swap_base(struct matrix_arr to_X, struct matrix_arr from_X, int i1, int i2);
void mul_base(struct matrix_arr to_X, struct matrix_arr from_X, int i, mpq_t lam);

///////////matrix construction functions////////////////////////////////////////////////////////////////////

void set_unit(int i0, int j0, int i_range, int j_range, struct matrix *mat);
void set_submatrix(int i0_source, int j0_source, int i0_target, int j0_target, int i_range, int j_range, struct matrix *source, struct matrix *target);
void set_diag_p_powers(int p, int i0_target, int j0_target, int range, int *source, struct matrix *target);


////////////////////////////////////////////////////////////////////////////////////

#endif