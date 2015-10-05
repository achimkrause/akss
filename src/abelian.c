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

//computes kernel of map f. f can't be NULL. if x or y are NULL, 
// the corresponding groups are taken to be free of the correct rank.
// if to_X, from_X are NULL, the corresponding lists are taken to be empty.
// if k, to_K or from_K are NULL, the corresponding return value is discarded.
void kernel(const int p, const Matrix *f, const AbelianGroup *x, const AbelianGroup *y, const MatrixArray *to_X,
    const MatrixArray *from_X, AbelianGroup *k, MatrixArray *to_K, MatrixArray *from_K) {

  //first, compute for f*rel_x a lift l along rel_y. Then a map R_x -> R_y + X is given by (-l, rel_x)^t,
  //while a map R_y +X -> Y is given by (r_Y,f).
  //The condition on the g_i is that f*g_i lifts over rel_y. Denote these lifts by l_i. Then maps -> R_y + X are given by (-l_i, g_i)^t
  //Now compute the kernel of R_y + X -> Y (denote K), lift the map (-l, rel_x)^t as well as the (-l_i, g_i)^t.
  //the lift of the rel thing determines the orders in K. Smith it, transforming the maps to K as well as the composition K -> R_y+X -> X accordingly.
  //we can drop the generators corresponding to order 1, killing the corresponding columns in the map from K and the corresponding rows in the maps to K.
  //we obtain the orders of a final K, an inclusion map K->X, lifts g_lift_i.

  mpq_t lam_col;
  mpq_t lam_row;
  mpq_init(lam_col);  
  mpq_init(lam_row);  
  mpz_set_ui(mpq_denref(lam_col), 1);
  mpz_set_ui(mpq_denref(lam_row), 1);

  //to construct -l, we multiply the rows of f with rel_y, then divide the columns by rel_x.
  int tor_rank_y;
  if(y != NULL) {
    tor_rank_y = 0;
  }
  else {
    tor_rank_y = y->tor_rank;
  }

  int tor_rank_x;
  if(x != NULL) {
    tor_rank_x = 0;
  }
  else {
    tor_rank_x = x->tor_rank;
  }
  Matrix *rel_free = matrix_alloc();
  matrix_init(rel_free, tor_rank_y + f->width, tor_rank_x); //build from -l and rel_x.


  for (int i = 0; i < tor_rank_y; i++) {

    mpz_ui_pow_ui(mpq_numref(lam_row), p, y->orders[i]);
    for (int j = 0; j < tor_rank_x; j++) {

      mpz_ui_pow_ui(mpq_numref(lam_col), p, x->orders[j]);  //could replace this by a second loop, traversing column-wise.(depends on how expensive the powering is).
      mpq_mul(rel_free->entries[i*rel_free->width + j], 
               f->entries[i*f->width+j], lam_col);
      mpq_div(rel_free->entries[i*rel_free->width+j], 
               rel_free->entries[i*rel_free->width+j], lam_row);
      mpq_neg(rel_free->entries[i*rel_free->width+j], 
               rel_free->entries[i*rel_free->width+j]);

    }
  }
  //set the second block of rel_free to rel_x
  if(x!= NULL) { 
    set_diag_p_powers(p, tor_rank_y, 0, tor_rank_x, x->orders, rel_free);
  }

  Matrix *f_free = matrix_alloc();
  matrix_init(f_free, f->height, tor_rank_y + f->width);  //build f_free from rel_y, f.

  if(y!=NULL) {
    set_diag_p_powers(p, 0, 0, y->tor_rank, y->orders, f_free);
  }
  set_submatrix(0, 0, 0, tor_rank_y, f->height, f->width, f, f_free);

  //next, compute the l_i by lifting f*g_i over rel_y.

  Matrix *f_to_X_i = matrix_alloc();
  MatrixArray *to_Ry_X = matrix_arr_alloc();
  matrix_arr_init(to_Ry_X, to_X->length + 1);  //reserve additional entry for rel_free!

  for (int i = 0; i < to_X->length; i++) {
    compose(f, to_X->entries[i], f_to_X_i);
    to_Ry_X->entries[i]=matrix_alloc();
    matrix_init(to_Ry_X->entries[i],tor_rank_y + f->width, to_X->entries[i]->width);

    //this loop lifts f*to_X_i over rel_y to create l_i, 
    //and installs that as first block in to_Ry_X.
    for (int row = 0; row < tor_rank_y; row++) {
      mpz_ui_pow_ui(mpq_numref(lam_row), p, y->orders[row]);

      for (int col = 0; col < to_X->entries[i]->width; col++) {
        mpq_div(to_Ry_X->entries[i]->entries[row*to_Ry_X->entries[i]->width + col], 
                 f_to_X_i->entries[row*f_to_X_i->width + col], lam_row);
        mpq_neg(to_Ry_X->entries[i]->entries[row*to_Ry_X->entries[i]->width + col],
                 to_Ry_X->entries[i]->entries[row*to_Ry_X->entries[i]->width + col]);
      }
    }

    //set second block of to_Ry_X to to_X.
    set_submatrix(0, 0, tor_rank_y, 0, to_X->entries[i]->height, to_X->entries[i]->width,
        to_X->entries[i], to_Ry_X->entries[i]);

    //we have to clear f_to_X_i by hand since matrix_clear frees the pointer itself, and we want to reuse it.
    for (int d = 0; d < f_to_X_i->height * f_to_X_i->width; d++) {
      mpq_clear(f_to_X_i->entries[d]);
    }
    free(f_to_X_i->entries);
  }
  free(f_to_X_i);

  to_Ry_X->entries[to_Ry_X->length - 1] = rel_free;

  MatrixArray *from_Ry_X = matrix_arr_alloc();
  matrix_arr_init(from_Ry_X, from_X->length);

  for (int i = 0; i < from_X->length; i++) {
    from_Ry_X->entries[i]=matrix_alloc();
    matrix_init(from_Ry_X->entries[i],from_X->entries[i]->height, tor_rank_y + f->width);
    set_submatrix(0, 0, 0, tor_rank_y, from_X->entries[i]->height, f->width, from_X->entries[i],
        from_Ry_X->entries[i]);
  }

  //now, smith f_free:
  smith(p, f_free, to_Ry_X, from_Ry_X, NULL, NULL);

  //we don't need rel_free in this list anymore since next call to smith will be on rel_free. MILDLY HACKY.
  to_Ry_X->length = to_Ry_X->length - 1;  

  int n = 0;
  mpq_t zero;
  mpq_init(zero);
  mpq_set_ui(zero, 0, 1);

  while (n < f_free->width && n < f_free->height
      && !mpq_equal(f_free->entries[n*(f_free->width+1)], zero)) {
    n++;
  }

  Matrix *rel_K = matrix_alloc();
  matrix_init(rel_K, rel_free->height - n, rel_free->width);
  set_submatrix(n, 0, 0, 0, rel_K->height, rel_K->width, rel_free, rel_K);

  MatrixArray *to_free_K = matrix_arr_alloc();
  matrix_arr_init(to_free_K, to_Ry_X->length);

  for (int i = 0; i < to_Ry_X->length; i++) {
    to_free_K->entries[i] = matrix_alloc();
    matrix_init(to_free_K->entries[i], rel_free->height - n, to_Ry_X->entries[i]->width);
    set_submatrix(n, 0, 0, 0, to_free_K->entries[i]->height, to_free_K->entries[i]->width,
        to_Ry_X->entries[i], to_free_K->entries[i]);
  }

  MatrixArray *from_free_K = matrix_arr_alloc();
  matrix_arr_init(from_free_K, from_Ry_X->length);

  for (int i = 0; i < from_Ry_X->length; i++) {
    from_free_K->entries[i] = matrix_alloc();
    matrix_init(from_free_K->entries[i],from_Ry_X->entries[i]->height, rel_free->height - n);
    set_submatrix(0, n, 0, 0, from_free_K->entries[i]->height, from_free_K->entries[i]->width,
        from_Ry_X->entries[i], from_free_K->entries[i]);
  }

  cokernel(p, rel_K, NULL, to_free_K, from_free_K, k, to_K, from_K);

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

void cokernel(const int p, const Matrix *f, const AbelianGroup *y, const MatrixArray *to_Y_0,
    const MatrixArray *from_Y_0, AbelianGroup *c, MatrixArray *to_C, MatrixArray *from_C) {

  Matrix *f_rel_y = matrix_alloc();

  int tor_rank_y =0;
  if(y!=NULL) {
    tor_rank_y=y->tor_rank;
  }
  matrix_init(f_rel_y, f->height, f->width + tor_rank_y); //build from f and rel_y.

  set_submatrix(0, 0, 0, 0, f->height, f->width, f, f_rel_y);
  if(y!=NULL) {
    set_diag_p_powers(p, 0, f->width, tor_rank_y, y->orders, f_rel_y);
  }
  //copy to_Y_0, from_Y_0.

  MatrixArray *to_Y = matrix_arr_copy(to_Y_0);
  MatrixArray *from_Y = matrix_arr_copy(from_Y_0);

  //now smith

  smith(p, f_rel_y, NULL, NULL, to_Y, from_Y);

  mpq_t one;
  mpq_init(one);
  mpq_set_ui(one, 1, 1);

  int n = 0;

  while (n < f_rel_y->height && n < f_rel_y->width
      && mpq_equal(f_rel_y->entries[n*(f_rel_y->width +1)], one)) {
    n++;
  }
  mpq_clear(one);

  //n is now the first non-1 entry.
  //so the g_factor guys come from the columns of the g beginning with n,
  //and the projection is given by the rows beginning with n.
  if (c != NULL) {
    mpq_t zero;
    mpq_init(zero);
    mpq_set_ui(zero, 0, 1);

    int tor_rank_c = 0;
    while (tor_rank_c < f_rel_y->height - n && tor_rank_c < f_rel_y->width - n
        && !mpq_equal(f_rel_y->entries[n*(f_rel_y->width +1)], zero)) {
      tor_rank_c++;
    }

    abelian_init(c, tor_rank_c, f_rel_y->height - n - tor_rank_c);

    for (int i = 0; i < c->tor_rank; i++) {

      c->orders[i] = valuation(p, f_rel_y->entries[(n+i)*(f_rel_y->width + 1)]);
    }

    mpq_clear(zero);
  }

  if (to_C != NULL) {
    matrix_arr_init(to_C, to_Y->length);

    for (int i = 0; i < to_C->length; i++) {
      to_C->entries[i] = matrix_alloc();
      matrix_init(to_C->entries[i],to_Y->entries[i]->height - n,
          to_Y->entries[i]->width);
      set_submatrix(n, 0, 0, 0, to_C->entries[i]->height, to_C->entries[i]->width,
          to_Y->entries[i], to_C->entries[i]);
    }
  }

  if (from_C != NULL) {
    matrix_arr_init(from_C, from_Y->length);

    for (int i = 0; i < from_C->length; i++) {
      from_C->entries[i] = matrix_alloc();
      matrix_init(from_C->entries[i],from_Y->entries[i]->height,
          from_Y->entries[i]->width - n);
      set_submatrix(0, n, 0, 0, from_C->entries[i]->height, from_C->entries[i]->width,
          from_Y->entries[i], from_C->entries[i]);
    }
  }

  matrix_clear(f_rel_y);
  matrix_arr_clear(from_Y);
  matrix_arr_clear(to_Y);
}

void epi_mono(const int p, const Matrix *f, const AbelianGroup *x, const AbelianGroup *y,
    AbelianGroup *img, Matrix *proj, Matrix *inc) {
  //first, compute kernel of f, tracking the inclusion K -> X.
  //then, compute the cokernel of K -> X, tracking the projection X -> C and a factorization of f over it.

  MatrixArray *from_X =matrix_arr_alloc();
  matrix_arr_init(from_X, 1);
  from_X->entries[0]=matrix_alloc();
  matrix_init(from_X->entries[0],f->width, f->width);
  set_unit(from_X->entries[0]);

  MatrixArray *from_K = matrix_arr_alloc();

  kernel(p, f, x, y, NULL, from_X, NULL, NULL, from_K); 

  MatrixArray *to_X_2 = matrix_arr_alloc();
  matrix_arr_init(to_X_2, 1);
  to_X_2->entries[0] = matrix_alloc();
  matrix_init(to_X_2->entries[0], f->width, f->width);
  set_unit(to_X_2->entries[0]);

  MatrixArray *from_X_2 = matrix_arr_alloc();
  matrix_arr_init(from_X_2, 1);
  from_X_2->entries[0] = matrix_copy(f);

  MatrixArray *to_C = matrix_arr_alloc();
  MatrixArray *from_C = matrix_arr_alloc();
  cokernel(p, from_K->entries[0], x, to_X_2, from_X_2, img, to_C, from_C);

  matrix_copy_to(to_C->entries[0],proj);
  matrix_copy_to(from_C->entries[0],inc);

  matrix_arr_clear(to_X_2);
  matrix_arr_clear(from_X);
  matrix_arr_clear(from_X_2);
  matrix_arr_clear(from_K);
  matrix_arr_clear(to_C);
  matrix_arr_clear(from_C);
}

/*void compose_diag_p_power(int p, Matrix *f, int *exponents, Matrix **res) {
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

}*/
  /*
void lift_diag_p_power(int p, Matrix *f, int *exponents, Matrix **res) {
  //    *exponents is supposed to be the normalized orders of a Z_(p) module, i.e.:
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

  / mpq_t zero;//DEBUG ONLY
   mpq_init(zero);//DEBUG ONLY
   mpq_set_ui(zero, 0, 1);//DEBUG ONLY
   /

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

  / for(int i=m; i<f->height; i++)//DEBUG ONLY
   {
   assert(mpq_equal(*entry_f,zero));
   }/

  mpq_clear(power);
  //mpq_clear(zero);//DEBUG ONLY
}*/

/*
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
}*/

void abelian_init(AbelianGroup *x, const int tor_rank, const int free_rank) {
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
