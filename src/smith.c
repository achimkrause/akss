#include "smith.h"
#include "test.h" //only for debug

////p-valuation helper functions////////////////////////////////////////////////////////////////////////7
int valuation_int(const int p, const mpz_t z) {
  //fprintf(stderr, "valuation_int called on %u, %u\n", p, mpz_get_ui(z));
  mpz_t rem;
  mpz_init(rem);
  int result = 0;
  mpz_set(rem, z);
  while (mpz_divisible_ui_p(rem, p)) {
    result++;
    mpz_divexact_ui(rem, rem, p);
  }
  mpz_clear(rem);
  return result;
}

int valuation(const int p, const mpq_t lam) {
  mpz_t z;
  mpz_init(z);
  mpz_set(z, mpq_numref(lam));
  int result;
  result = valuation_int(p, z);
  if (result > 0) {
    mpz_clear(z);
    return result;
  } else {
    mpz_set(z, mpq_denref(lam));
    result = -valuation_int(p, z);
    mpz_clear(z);
    return result;
  }
}


void smith(const int p, Matrix *mat, MatrixArray *to_X,
    MatrixArray *from_X_0, MatrixArray *to_Y_0,
    MatrixArray *from_Y) {
  //mat represents a map X -> Y. to_X, from_X, to_Y, from_Y are representatives for maps to/from X/Y.
  //smith brings mat into smith normal form by applying base changes to X and Y and transforming to_X and to_Y accordingly.
  mpq_t min;
  mpq_init(min);
  mpq_t lam;
  mpq_init(lam);

  MatrixArray *from_X = matrix_arr_alloc();
  MatrixArray *to_Y = matrix_arr_alloc();

  if(from_X_0!=NULL) {

    matrix_arr_init(from_X, from_X_0->length+1);

    for (int i = 0; i < from_X_0->length; i++) {
      from_X->entries[i] = from_X_0->entries[i];
    }  
  }
  else {
    matrix_arr_init(from_X,1);
  }
  
  if(to_Y_0!=NULL) {
    matrix_arr_init(to_Y, to_Y_0->length+1);
    for (int i = 0; i < to_Y_0->length; i++) {
      to_Y->entries[i] = to_Y_0->entries[i];
    }
  }
  else {
    matrix_arr_init(to_Y,1);
  }

  from_X->entries[from_X->length -1] = mat;
  to_Y->entries[to_Y->length-1] = mat;

  int block = 0;

  mpq_t zero;
  mpq_init(zero);
  mpq_set_ui(zero,0,1);

  while (block < mat->width && block < mat->height) {
    //find smallest in block:
    int i_min = -1;
    int j_min = -1;
    int min_val = 0;
    for (int i = block; i < mat->height; i++) {
      for (int j = block; j < mat->width; j++) {
        if (!mpq_equal(mat->entries[i*mat->width+j],zero)) {
          if (i_min == -1) {
            i_min = i;
            j_min = j;
            mpq_set(min, mat->entries[i*mat->width+j]);
            min_val = valuation(p, min);
          } else {
            if (valuation(p, mat->entries[i*mat->width+j]) < min_val) {
              i_min = i;
              j_min = j;
              mpq_set(min, mat->entries[i*mat->width+j]);
              min_val = valuation(p, min);
            }
          }
        }
      }
    }
    if (i_min == -1) {
      mpq_clear(min);
      mpq_clear(lam);
      //Caller still owns mat and the entries of from_X_0, to_Y_0, so we can't clear them.
      free(from_X->entries);
      free(from_X);
      free(to_Y->entries);
      free(to_Y);
      return;
    }

    for (int i = block; i < mat->height; i++) {
      if (i != i_min) {
        mpq_div(lam, mat->entries[i*mat->width+j_min], min);
        add_base(to_Y, from_Y, i, i_min, lam);
        //has the effect of subtracting lam times the i_min-row from the i-row.
      }
    }

    for (int j = block; j < mat->width; j++) {
      if (j != j_min) {
        mpq_div(lam, mat->entries[i_min*mat->width + j], min);
        mpq_neg(lam, lam);
        add_base(to_X, from_X, j_min, j, lam);
        //has the effect of adding lam times the j_min col to the j-col.
      }
    }
    swap_base(to_Y, from_Y, i_min, block);
    swap_base(to_X, from_X, j_min, block);

    mpz_ui_pow_ui(mpq_numref(lam), p, min_val);
    mpz_set_ui(mpq_denref(lam), 1);
    mpq_div(lam, lam, min);
    mul_base(to_X, from_X, block, lam);

    block++;
  }
  mpq_clear(min);
  mpq_clear(lam);
  //Caller still owns mat and the entries of from_X_0, to_Y_0, so we can't clear them.
  free(from_X->entries); 
  free(from_X);
  free(to_Y->entries);
  free(to_Y);
}

