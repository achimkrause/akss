#include "matrix.h"

Matrix* matrix_alloc()
{
  return((Matrix*) malloc(sizeof(Matrix)));
}

void matrix_init(Matrix *mat, const int height, const int width) {
  mat->height = height;
  mat->width = width;
  mat->entries = (mpq_t*) malloc(sizeof(mpq_t) * height * width);
  for (int i = 0; i < height * width; i++) {
    mpq_init(mat->entries[i]);
  }
}

void matrix_clear_entries(Matrix *mat) {
  for (int i = 0; i < mat->height * mat->width; i++) {
    mpq_clear(mat->entries[i]);
  }
  free(mat->entries);
}

void matrix_clear(Matrix *mat) {
  matrix_clear_entries(mat);
  free(mat);
}

Matrix* matrix_copy(const Matrix *mat) {
  Matrix* new_mat = matrix_alloc();
  matrix_copy_to(mat, new_mat);
  return new_mat;
}

void matrix_copy_to(const Matrix *source, Matrix *target) {
  matrix_init(target, source->height, source->width);
  for(int i=0; i<source->height*source->width; i++) {
    mpq_set(target->entries[i], source->entries[i]);
  }
}

void matrix_arr_clear(MatrixArray *arr) {
  
  for (int i = 0; i < arr->length; i++) {
    matrix_clear_entries(arr->entries[i]);
  free(arr->entries);
  free(arr);
  }
}


//copies a MatrixArray (also copies the entries, not just references.)
MatrixArray* matrix_arr_copy(const MatrixArray *arr){
  if(arr == NULL) {
    return NULL;
  }
  MatrixArray *result;
  result = matrix_arr_alloc();
  matrix_arr_init(result, arr->length);
  for(int i=0; i<arr->length; i++)
  {
    result->entries[i] = matrix_copy(arr->entries[i]);
  }
  return result;
}

MatrixArray* matrix_arr_alloc() {
  return((MatrixArray*) malloc(sizeof(MatrixArray)));
}

void matrix_arr_init(MatrixArray *arr, const int length) {

  arr->length = length;
  arr->entries = malloc(sizeof(Matrix*) * length);
}


void add_row(Matrix *mat, const int i1, const int i2, const mpq_t lam) {
  mpq_t prod;
  mpq_init(prod);
  mpq_t *row1_entry = mat->entries + i1 * (mat->width);
  mpq_t *row2_entry = mat->entries + i2 * (mat->width);

  for (int j = 0; j < mat->width; j++) {
    mpq_mul(prod, *row1_entry, lam);
    mpq_add(*row2_entry, *row2_entry, prod);
    row1_entry++;
    row2_entry++;
  }
  mpq_clear(prod);
}

void add_col(Matrix *mat, const int j1, const int j2, const mpq_t lam) {
  mpq_t prod;
  mpq_init(prod);

  mpq_t *col1_entry = mat->entries + j1;
  mpq_t *col2_entry = mat->entries + j2;

  for (int i = 0; i < mat->height; i++) {
    mpq_mul(prod, *col1_entry, lam);
    mpq_add(*col2_entry, *col2_entry, prod);
    col1_entry += mat->width;
    col2_entry += mat->width;

  }
  mpq_clear(prod);
}

void sub_row(Matrix *mat, const int i1, const int i2, const mpq_t lam) {
  mpq_t prod;
  mpq_init(prod);
  mpq_t *row1_entry = mat->entries + i1 * (mat->width);
  mpq_t *row2_entry = mat->entries + i2 * (mat->width);

  for (int j = 0; j < mat->width; j++) {
    mpq_mul(prod, *row1_entry, lam);
    mpq_sub(*row2_entry, *row2_entry, prod);
    row1_entry++;
    row2_entry++;
  }
  mpq_clear(prod);
}

void sub_col(Matrix *mat, const int j1, const int j2, const mpq_t lam) {
  mpq_t prod;
  mpq_init(prod);

  mpq_t *col1_entry = mat->entries + j1;
  mpq_t *col2_entry = mat->entries + j2;

  for (int i = 0; i < mat->height; i++) {
    mpq_mul(prod, *col1_entry, lam);
    mpq_sub(*col2_entry, *col2_entry, prod);
    col1_entry += mat->width;
    col2_entry += mat->width;

  }
  mpq_clear(prod);
}

void swap_row(Matrix *mat, const int i1, const int i2) {
  mpq_t *row1_entry = mat->entries + i1 * (mat->width);
  mpq_t *row2_entry = mat->entries + i2 * (mat->width);

  for (int j = 0; j < mat->width; j++) {
    mpq_swap(*row1_entry, *row2_entry);
    row1_entry++;
    row2_entry++;
  }
}

void swap_col(Matrix *mat, const int j1, const int j2) {
  mpq_t *col1_entry = mat->entries + j1;
  mpq_t *col2_entry = mat->entries + j2;

  for (int i = 0; i < mat->height; i++) {

    mpq_swap(*col1_entry, *col2_entry);
    col1_entry += mat->width;
    col2_entry += mat->width;

  }
}

void mul_row(Matrix *mat, const int i, const mpq_t lam) {
  mpq_t *row_entry = mat->entries + i * (mat->width);
  for (int j = 0; j < mat->width; j++) {
    mpq_mul(*row_entry, *row_entry, lam);
    row_entry++;
  }
}

void mul_col(Matrix *mat, const int j, const mpq_t lam) {
  mpq_t *col_entry = mat->entries + j;
  for (int i = 0; i < mat->height; i++) {
    mpq_mul(*col_entry, *col_entry, lam);
    col_entry += mat->width;
  }
}

void div_row(Matrix *mat, const int i, const mpq_t lam) {
  mpq_t *row_entry = mat->entries + i * (mat->width);
  for (int j = 0; j < mat->width; j++) {
    mpq_div(*row_entry, *row_entry, lam);
    row_entry++;
  }
}

void div_col(Matrix *mat, const int j, const mpq_t lam) {
  mpq_t *col_entry = mat->entries + j;
  for (int i = 0; i < mat->height; i++) {
    mpq_div(*col_entry, *col_entry, lam);
    col_entry += mat->width;
  }
}


void add_base(MatrixArray *to_X, MatrixArray *from_X, int i1, int i2,
    mpq_t lam) {
  //For a free module X with basis x_i, modifies that basis by setting x_i2' = x_i2 + lam*x_i1, x_i' = x_i else. Base-changes maps from and to X accordingly, passed in lists.
  //on maps to X, this has the following effect:
  // expanding f(x) = a_i x_i = a_i' x_i', we see that a_i1 = a_i1' + lam*a_i2', and a_i = a_i' else, so a_i1' = a_i1-lam*a_i2.
  //     ==>  We add -lam times the i2-th row to the i1-th row.
  //on maps from X:
  // expanding f(x_i2') = f(x_i2) + lam*f(x_i1)
  //     ==> We add lam times the i1-th row to the i2-th row.

  if(to_X!=NULL) {
    for (int i = 0; i < to_X->length; i++) {
      sub_row(to_X->entries[i], i2, i1, lam);
    }
  }

  if(from_X!=NULL) {
    for (int i = 0; i < from_X->length; i++) {
      add_col(from_X->entries[i], i1, i2, lam);
    }
  } 
}

void swap_base(MatrixArray *to_X, MatrixArray *from_X, int i1, int i2) {
  //for X a free module with basis x_i, this modifies the basis by setting x_i2'=x_i1, x_i1'=x_i2. Base-changes maps accordingly:
  // maps TO X have their rows i1 and i2 swapped.
  // maps FROM X have their columns i1 and i2 swapped.
  if(to_X!=NULL) {
    for (int i = 0; i < to_X->length; i++) {
      swap_row(to_X->entries[i], i1, i2);
    }
  }

  if(from_X!=NULL) {
    for (int i = 0; i < from_X->length; i++) {
      swap_col(from_X->entries[i], i1, i2);
    }
  } 
}

void mul_base(MatrixArray *to_X, MatrixArray *from_X, int i1,
    mpq_t lam) {
  //for X a free module with basis x_i, this modifies the basis by setting x_i' = lam*x_i. Base-changes maps accordingly:
  // maps to X:
  // f(x) = a_i x_i = a_i' x_i', so a_i = lam a_i', so a_i' = a_i/lam.
  // maps from X:
  // f(x_i') = lam*f(x_i)
  if(to_X!=NULL) {
    for (int i = 0; i < to_X->length; i++) {
      div_row(to_X->entries[i], i1, lam);
    } 
  }

  if(from_X!=NULL) {
    for (int i = 0; i < from_X->length; i++) {
      mul_col(from_X->entries[i], i1, lam);
    }
  } 
}

///////////matrix construction functions////////////////////////////////////////////////////////////////////

void set_unit(Matrix *mat) {
  set_unit_range(0, 0, mat->height, mat->width, mat);
}

void set_unit_range(int i0, int j0, int i_range, int j_range, Matrix *mat) {
  for (int i = i0; i < i0+i_range; i++) {
    for (int j = j0; j < j0+j_range; j++) {
      if (i == j) {
        mpq_set_ui(mat->entries[i*mat->width+j], 1, 1);
      } else {
        mpq_set_ui(mat->entries[i*mat->width+j], 0, 1);
      }
    }
  }
}

void set_submatrix(int i0_source, int j0_source, int i0_target, int j0_target,
    int i_range, int j_range, const Matrix *source, Matrix *target) {
  for (int i = 0; i < i_range; i++) {
    for (int j = 0; j < j_range; j++) {
      mpq_set(target->entries[(i0_target+i)*target->width+j0_target + j], 
               source->entries[(i0_source + i)*source->width + j0_source + j]);
    }
  }
}
//TODO: maybe add aliases like copy etc here.

void set_diag_p_powers(int p, int i0_target, int j0_target, int range,
    int *source, Matrix *target) {

  for (int i = 0; i < range; i++) {
    for (int j = 0; j < range; j++) {
      if (i == j) {
        mpz_ui_pow_ui(mpq_numref(target->entries[(i0_target + i)*target->width + j0_target + j]), 
                                  p,
                                  source[i]);
        mpz_set_ui(mpq_denref(target->entries[(i0_target + i)*target->width + j0_target + j]), 1);
      } else {
        mpq_set_ui(target->entries[(i0_target + i)*target->width + j0_target + j], 0, 1);
      }
    }
  }
}

void compose(const Matrix *g, const Matrix *f, Matrix *gf) {
  matrix_init(gf, g->height, f->width);
  mpq_t prod;
  mpq_init(prod);

  for (int i = 0; i < gf->height; i++) {
    for (int j = 0; j < gf->width; j++) {

      mpq_set_ui(gf->entries[i*gf->width+j],0,1);
      for (int k = 0; k < f->height; k++) {
        mpq_mul(prod, f->entries[k*f->width+j], g->entries[i*g->width + k]);
        mpq_add(gf->entries[i*gf->width+j], gf->entries[i*gf->width+j], prod);
      }
    }
  }
  mpq_clear(prod);
}
