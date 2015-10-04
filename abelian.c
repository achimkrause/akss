#include "abelian.h"

/*
 given A -> B and resolutions rA -> gA -> A, rB -> gB -> B, and a map given as gA -> gB
 (the resolutions are in Smith form and therefore just given by int-arrays. 
 The convention is that sth like (p^{i_1},...,p^{i_n},0,..,0) corresponds to the nxm matrix with entries the p-powers where n is the length and m is the length minus # of zeros).

 Lift to rmap: rA -> rB (easy: this just multiplies columns with the entries of orders_A and divides rows by entries of orders_B)
 Now form maps rA -> rB + gA -> gB. The first map is f=(rmap, -orders_A)^t, the second g=(orders_B, map).
 Compute homology: smithify f, and base-change g accordingly. The non-zero entries of f now lead to zero-columns in g. These are left alone (except for permutations) when we smithify g.
 If additional zero columns appear, we get free summands in homology, otherwise everything is torsion.
 So if we know to only obtain torsion, we can ignore the second run of smith.
 Similarly, if A is torsion-free, we need to compute just the kernel of A + rB -> gB, which takes only one application of smith as well.
 */

void kernel(int p, Matrix *f, AbelianGroup x, AbelianGroup y, MatrixArray to_X,
    MatrixArray from_X, AbelianGroup *k, MatrixArray *to_K, MatrixArray *from_K) {

  //first, compute for f*rel_x a lift l along rel_y. Then a map R_x -> R_y + X is given by (-l, rel_x)^t,
  //while a map R_y +X -> Y is given by (r_Y,f).
  //The condition on the g_i is that f*g_i lifts over rel_y. Denote these lifts by l_i. Then maps -> R_y + X are given by (-l_i, g_i)^t
  //Now compute the kernel of R_y + X -> Y (denote K), lift the map (-l, rel_x)^t as well as the (-l_i, g_i)^t.
  //the lift of the rel thing determines the orders in K. Smith it, transforming the maps to K as well as the composition K -> R_y+X -> X accordingly.
  //we can drop the generators corresponding to order 1, killing the corresponding columns in the map from K and the corresponding rows in the maps to K.
  //we obtain the orders of a final K, an inclusion map K->X, lifts g_lift_i.

  //multiply cols of f with p^orders_x, divide rows by p^orders_y
  Matrix *rel_free = matrix_init(y.tor_rank + f->width, x.tor_rank);  //build from -l and rel_x. ALLOCATES

  int *orders_y = y.orders;

  int *orders_x = x.orders;

  mpq_t *entry_rel_free = rel_free->entries;
  mpq_t *entry_f = f->entries;
  mpq_t lam_col;
  mpq_t lam_row;
  mpq_init(lam_col);  //ALLOCATES
  mpq_init(lam_row);  //ALLOCATES
  mpz_set_ui(mpq_denref(lam_col), 1);
  mpz_set_ui(mpq_denref(lam_row), 1);

  for (int i = 0; i < y.tor_rank; i++) {

    mpz_ui_pow_ui(mpq_numref(lam_row), p, *orders_y);
    for (int j = 0; j < x.tor_rank; j++) {

      mpz_ui_pow_ui(mpq_numref(lam_col), p, *orders_x);  //could replace this by a second loop, traversing column-wise.(depends on how expensive the powering is).
      mpq_mul(*entry_rel_free, *entry_f, lam_col);
      mpq_div(*entry_rel_free, *entry_rel_free, lam_row);
      mpq_neg(*entry_rel_free, *entry_rel_free);

      entry_rel_free++;
      entry_f++;
      orders_x++;
    }
    orders_x = x.orders;
    entry_f = entry_f - x.tor_rank + f->width;
    orders_y++;
  }

  set_diag_p_powers(p, y.tor_rank, 0, x.tor_rank, x.orders, rel_free);

  Matrix *f_free = matrix_init(f->height, y.tor_rank + f->width);  //build f_free from rel_x, f. ALLOCATES

  set_diag_p_powers(p, 0, 0, y.tor_rank, y.orders, f_free);
  set_submatrix(0, 0, 0, y.tor_rank, f->height, f->width, f, f_free);

  //next, compute the l_i by lifting f*g_i over rel_y.

  Matrix *f_to_X_i;

  MatrixArray to_Ry_X;
  matrix_arr_init(&to_Ry_X, to_X.length + 1);  //reserve additional entry for rel_free!

  Matrix **to_X_i = to_X.entries;
  Matrix **to_Ry_X_i = to_Ry_X.entries;

  for (int i = 0; i < to_X.length; i++) {
    compose(f, (*to_X_i), &f_to_X_i);

    *to_Ry_X_i = matrix_init(y.tor_rank + f->width, (*to_X_i)->width);

    mpq_t *entry_f_to_X = f_to_X_i->entries;
    mpq_t *entry_to_Ry_X = (*to_Ry_X_i)->entries;

    orders_y = y.orders;

    for (int row = 0; row < y.tor_rank; row++) {
      mpz_ui_pow_ui(mpq_numref(lam_row), p, *orders_y);

      for (int col = 0; col < (*to_X_i)->width; col++) {
        mpq_div(*entry_to_Ry_X, *entry_f_to_X, lam_row);
        mpq_neg(*entry_to_Ry_X, *entry_to_Ry_X);
        entry_to_Ry_X++;
        entry_f_to_X++;
      }
      orders_y++;
    }
    set_submatrix(0, 0, y.tor_rank, 0, (*to_X_i)->height, (*to_X_i)->width,
        (*to_X_i), (*to_Ry_X_i));

    matrix_clear(f_to_X_i);
    to_X_i++;
    to_Ry_X_i++;
  }
  *to_Ry_X_i = rel_free;

  MatrixArray from_Ry_X;
  matrix_arr_init(&from_Ry_X, from_X.length);
  Matrix **from_Ry_X_i = from_Ry_X.entries;
  Matrix **from_X_i = from_X.entries;
  for (int i = 0; i < from_X.length; i++) {
    *from_Ry_X_i = matrix_init((*from_X_i)->height, y.tor_rank + f->width);
    set_submatrix(0, 0, 0, y.tor_rank, (*from_X_i)->height, f->width, *from_X_i,
        *from_Ry_X_i);
  }

  MatrixArray null_arr;
  null_arr.length = 0;
  null_arr.entries = NULL;

  //now, smith f_free:
  smith(p, f_free, to_Ry_X, from_Ry_X, null_arr, null_arr);

  to_Ry_X.length = to_Ry_X.length - 1;  //we don't need rel_free in this list anymore since next call to smith will be on rel_free.

  int n = 0;
  mpq_t *entry_f_free = f_free->entries;
  mpq_t zero;
  mpq_init(zero);
  mpq_set_ui(zero, 0, 1);

  while (n < f_free->width && n < f_free->height
      && !mpq_equal(*entry_f_free, zero)) {
    n++;
    entry_f_free = entry_f_free + f_free->width + 1;
  }

  Matrix *rel_K = matrix_init(rel_free->height - n, rel_free->width);
  set_submatrix(n, 0, 0, 0, rel_K->height, rel_K->width, rel_free, rel_K);

  MatrixArray to_free_K;
  matrix_arr_init(&to_free_K, to_Ry_X.length);
  to_Ry_X_i = to_Ry_X.entries;
  Matrix **to_free_K_i = to_free_K.entries;

  for (int i = 0; i < to_Ry_X.length; i++) {
    *to_free_K_i = matrix_init(rel_free->height - n, (*to_Ry_X_i)->width);
    set_submatrix(n, 0, 0, 0, (*to_free_K_i)->height, (*to_free_K_i)->width,
        (*to_Ry_X_i), (*to_free_K_i));
    to_Ry_X_i++;
    to_free_K_i++;
  }

  MatrixArray from_free_K;
  matrix_arr_init(&from_free_K, from_Ry_X.length);

  from_Ry_X_i = from_Ry_X.entries;
  Matrix **from_free_K_i = from_free_K.entries;

  for (int i = 0; i < from_Ry_X.length; i++) {
    *from_free_K_i = matrix_init((*from_Ry_X_i)->height, rel_free->height - n);
    set_submatrix(0, n, 0, 0, (*from_free_K_i)->height, (*from_free_K_i)->width,
        (*from_Ry_X_i), (*from_free_K_i));
    from_Ry_X_i++;
    from_free_K_i++;
  }

  AbelianGroup free_K;
  abelian_init(&free_K, 0, rel_free->height - n);

  cokernel(p, rel_K, free_K, to_free_K, from_free_K, k, to_K, from_K);
  abelian_clear(&free_K);

  matrix_clear(rel_free);
  matrix_clear(rel_K);
  matrix_clear(f_free);
  mpq_clear(lam_col);
  mpq_clear(lam_row);
  mpq_clear(zero);

  matrix_arr_clear(to_free_K);
  matrix_arr_clear(from_free_K);
  matrix_arr_clear(from_Ry_X);
  matrix_arr_clear(to_Ry_X);
}

void cokernel(int p, Matrix *f, AbelianGroup y, MatrixArray to_Y,
    MatrixArray from_Y, AbelianGroup *c, MatrixArray *to_C, MatrixArray *from_C) {

  Matrix *f_rel_y = matrix_init(f->height, f->width + y.tor_rank);

  set_submatrix(0, 0, 0, 0, f->height, f->width, f, f_rel_y);
  set_diag_p_powers(p, 0, f->width, y.tor_rank, y.orders, f_rel_y);

  //now smith
  MatrixArray null_arr;
  null_arr.length = 0;
  null_arr.entries = NULL;

  smith(p, f_rel_y, null_arr, null_arr, to_Y, from_Y);

  mpq_t one;
  mpq_init(one);
  mpq_set_ui(one, 1, 1);

  mpq_t *entry_diag = f_rel_y->entries;
  int n = 0;

  while (n < f_rel_y->height && n < f_rel_y->width
      && mpq_equal(*entry_diag, one)) {
    n++;
    entry_diag = entry_diag + f_rel_y->width + 1;
  }

  //n is now the first non-1 entry.
  //so the g_factor guys come from the columns of the g beginning with n,
  //and the projection is given by the rows beginning with n.
  if (c != NULL) {

    mpq_t zero;
    mpq_init(zero);
    mpq_set_ui(zero, 0, 1);

    int tor_rank_c = 0;
    while (tor_rank_c < f_rel_y->height - n && tor_rank_c < f_rel_y->width - n
        && !mpq_equal(*entry_diag, zero)) {
      tor_rank_c++;
      entry_diag = entry_diag + f_rel_y->width + 1;
    }

    abelian_init(c, tor_rank_c, f_rel_y->height - n - tor_rank_c);

    int *orders_entry = c->orders;

    entry_diag = f_rel_y->entries + n + n * (f_rel_y->width);
    for (int i = 0; i < c->tor_rank; i++) {

      *orders_entry = valuation(p, *entry_diag);

      orders_entry++;
      entry_diag = entry_diag + f_rel_y->width + 1;
    }

    mpq_clear(zero);
  }

  if (to_C != NULL) {
    matrix_arr_init(to_C, to_Y.length);

    Matrix **entry_to_C = to_C->entries;
    Matrix **entry_to_Y = to_Y.entries;

    for (int i = 0; i < to_C->length; i++) {
      *entry_to_C = matrix_init((*entry_to_Y)->height - n,
          (*entry_to_Y)->width);
      set_submatrix(n, 0, 0, 0, (*entry_to_C)->height, (*entry_to_C)->width,
          (*entry_to_Y), (*entry_to_C));
      entry_to_C++;
    }
  }

  if (from_C != NULL) {
    matrix_arr_init(from_C, from_Y.length);

    Matrix **entry_from_C = from_C->entries;
    Matrix **entry_from_Y = from_Y.entries;

    for (int i = 0; i < from_C->length; i++) {
      *entry_from_C = matrix_init((*entry_from_Y)->height,
          (*entry_from_Y)->width - n);
      set_submatrix(0, n, 0, 0, (*entry_from_C)->height, (*entry_from_C)->width,
          (*entry_from_Y), (*entry_from_C));
      entry_from_C++;
    }
  }

  matrix_clear(f_rel_y);
  mpq_clear(one);
}

void epi_mono(int p, Matrix *f, AbelianGroup x, AbelianGroup y,
    AbelianGroup *img, Matrix **proj, Matrix **inc) {
  //first, compute kernel of f, tracking the inclusion K -> X.
  //then, compute the cokernel of K -> X, tracking the projection X -> C and a factorization of f over it.

  MatrixArray from_X;
  matrix_arr_init(&from_X, 1);
  from_X.entries[0] = matrix_init(f->width, f->width);
  set_unit(from_X.entries[0]);

  MatrixArray to_X;
  matrix_arr_init(&to_X, 0);

  MatrixArray from_K;

  kernel(p, f, x, y, to_X, from_X, NULL, NULL, &from_K);  //check that this doesn't change f (I think it doesn't).

  MatrixArray to_X_2;
  matrix_arr_init(&to_X_2, 1);
  to_X_2.entries[0] = matrix_init(f->width, f->width);
  set_unit(to_X_2.entries[0]);

  MatrixArray from_X_2;
  matrix_arr_init(&from_X_2, 1);
  from_X_2.entries[0] = f;

  MatrixArray to_C;
  MatrixArray from_C;
  cokernel(p, from_K.entries[0], x, to_X_2, from_X_2, img, &to_C, &from_C);

  *proj = matrix_copy(to_C.entries[0]);
  *inc = matrix_copy(from_C.entries[0]);

  matrix_arr_clear(to_X);
  matrix_arr_clear(to_X_2);
  matrix_arr_clear(from_X);
  matrix_arr_clear(from_X_2);
  matrix_arr_clear(from_K);
  matrix_arr_clear(to_C);
  matrix_arr_clear(from_C);
}

void compose_diag_p_power(int p, Matrix *f, int *exponents, Matrix **res) {
  int m = 0;
  int *exps = exponents;
  while (*exps != 0 && m < f->width) {
    m++;
    exps++;
  }

  *res = matrix_init(f->height, m);
  mpq_t *entry_res = (*res)->entries;
  mpq_t *entry_f = f->entries;
  exps = exponents;

  mpq_t power;
  mpq_init(power);
  mpz_set_ui(mpq_denref(power), 1);

  for (int j = 0; j < m; j++) {
    mpz_ui_pow_ui(mpq_numref(power), p, *exponents);
    for (int i = 0; i < f->height; i++) {
      mpq_mul(*entry_res, *entry_f, power);
      entry_res += m;
      entry_f += f->width;
    }
    entry_res = entry_res - m * f->height + 1;
    entry_f = entry_f - f->width * f->height + 1;
  }

  mpq_clear(power);

}
void lift_diag_p_power(int p, Matrix *f, int *exponents, Matrix **res) {
  //*exponents is supposed to be the normalized orders of a Z_(p) module, i.e.:
  // -no p^0 appear
  // -0 codifies actual zeroes. Those are supposed to be grouped in the end
  //The lift is computed against the nxm matrix where n is the full rank of f's target (f->height), m is the number of nonzero entries in *exponents
  // and the entries of this nxm matrix are the nonzero entries in *exponents.
  int m = 0;
  int *exps = exponents;
  while (*exps != 0 && m < f->height) {
    m++;
    exps++;
  }

  *res = matrix_init(m, f->width);

  mpq_t *entry_res = (*res)->entries;
  mpq_t *entry_f = f->entries;
  exps = exponents;

  /*mpq_t zero;//DEBUG ONLY
   mpq_init(zero);//DEBUG ONLY
   mpq_set_ui(zero, 0, 1);//DEBUG ONLY
   */

  mpq_t power;
  mpq_init(power);
  mpz_set_ui(mpq_denref(power), 1);

  for (int i = 0; i < m; i++) {
    mpz_ui_pow_ui(mpq_numref(power), p, *exps);

    if (*exps == 0) {
      for (int j = 0; j < f->width; j++) {

        //assert(mpq_equal(*entry_f,zero));
        mpq_set_ui(*entry_res, 0, 1);

        entry_res++;
        entry_f++;
      }
    } else {
      for (int j = 0; j < f->width; j++) {
        //assert(*exps <= valuation(p, *entry_f));
        mpq_div(*entry_res, *entry_f, power);

        entry_res++;
        entry_f++;
      }
    }
    exps++;
  }

  /*for(int i=m; i<f->height; i++)//DEBUG ONLY
   {
   assert(mpq_equal(*entry_f,zero));
   }*/

  mpq_clear(power);
  //mpq_clear(zero);//DEBUG ONLY
}

void compose(Matrix *g, Matrix *f, Matrix **gf) {
  *gf = matrix_init(g->height, f->width);
  Matrix *result = *gf;
  mpq_t *entry_gf = result->entries;
  mpq_t prod;
  mpq_init(prod);

  mpq_t *entry_g = g->entries;
  mpq_t *entry_f = f->entries;

  for (int i = 0; i < result->height; i++) {
    for (int j = 0; j < result->width; j++) {

      for (int k = 0; k < f->height; k++) {

        //fprintf(stderr, "i=%u, j=%u, k=%u\n",i,j,k);

        mpq_mul(prod, *entry_f, *entry_g);
        mpq_add(*entry_gf, *entry_gf, prod);
        entry_g++;
        entry_f += f->width;
      }
      entry_gf++;
      entry_g = entry_g - (g->width);  //set g back to beginning of row.
      entry_f = entry_f - (f->height) * (f->width) + 1;  //set f back to beginning of col, add 1.
    }
    //here, f is 1 past the first row, set back.
    entry_f = entry_f - (f->width);
    //g needs to be set to the next row.
    entry_g = entry_g + (g->width);
  }
  mpq_clear(prod);
}

int lift_diag(Matrix *f, Matrix *d, Matrix **res) {
  Matrix *lift = matrix_init(d->width, f->width);
  mpq_t *entry_d = d->entries;
  mpq_t *entry_f = f->entries;
  mpq_t *entry_l = lift->entries;

  for (int i = 0; i < f->height; i++) {
    if (i < d->width) {
      for (int j = 0; j < f->width; j++) {
        if (mpq_sgn(*entry_d) == 0) {
          if (mpq_sgn(*entry_f) != 0) {
            matrix_clear(lift);
            return 0;
          }
        } else {
          mpq_div(*entry_l, *entry_f, *entry_d);
        }
        entry_f++;
        entry_l++;
      }
    } else {
      for (int j = 0; j < f->width; j++) {
        if (mpq_sgn(*entry_f) != 0) {
          matrix_clear(lift);
          return 0;
        }
        entry_f++;
      }
    }
    entry_d = entry_d + (d->width) + 1;
  }

  (*res) = lift;
  return 1;
}

void abelian_init(AbelianGroup *x, int tor_rank, int free_rank) {
  x->tor_rank = tor_rank;
  x->free_rank = free_rank;
  x->orders = malloc(sizeof(int) * tor_rank);
}

void abelian_clear(AbelianGroup *x) {
  free(x->orders);
}
/*
 int lift(int p, Matrix *f, Matrix *g, Matrix **res)
 {
 //lifts f along g, i.e. computes res with g*res = f, or returns 0.
 //first, smith g, putting f into to_Y, id into from_X    (g'= a*g*b,f' =a*f. Lift f' along g' to l, a*g*b*l = a*f, so g*b*l = f. Compute b*l)



 Matrix *g_copy;
 Matrix *f_copy;
 Matrix *b;
 g_copy = copy(g);
 f_copy = copy(f);
 b      = matrix_init(g->width, g->width);
 set_unit(b);

 matrix_ll *f_ll = matrix_ll_init(f_copy);
 matrix_ll *b_ll = matrix_ll_init(b);

 smith(p,g_copy,NULL,b_ll,f_ll,NULL);

 free(f_ll);
 free(b_ll);

 mpq_t *entry_f = f_copy->entries;

 Matrix *lift;

 if(!lift_diag(f_copy,g_copy,&lift))
 {
 clear(f_copy);
 clear(g_copy);
 clear(b);
 return 0;
 }

 *res = matrix_init(g->width, f->width);

 compose(b,lift,&res);//OPTIMIZE: could do inplace compose!

 clear(g_copy);
 clear(f_copy);
 clear(b);
 clear(lift);
 return 1;
 }*/
