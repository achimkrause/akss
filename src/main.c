#include <stdio.h>
#include <gmp.h>
#include <assert.h>
#include <stdlib.h>

#include "smith.h"
#include "test.h"

int main() {
  /* AbelianGroup x;
   AbelianGroup y;

   abelian_init(&x, 2, 0);
   *(x.orders)=1;
   *(x.orders+1)=1;
   abelian_init(&y, 1, 0);
   *(y.orders)=2;

   Matrix *f = matrix_init(1,2);
   int val[2] = {2,2};
   fill_matrix(val, f);

   Matrix *g = matrix_init(2,1);
   int val_g[2] = {1,1};
   fill_matrix(val_g, g);

   test_kernel(2, f, x, y, g);*/

  AbelianGroup x;
  AbelianGroup y;

  abelian_init(&x, 0, 2);
  abelian_init(&y, 0, 2);

  Matrix *f = matrix_init(2, 2);
  int val[4] = { 0, 1, 0, 0 };
  fill_matrix(val, f);

  test_epi_mono(2, f, x, y);

  abelian_clear(&x);
  abelian_clear(&y);

  return 0;
}

/*
 void test_val()
 {
 mpz_t test;
 mpq_t testfrac;
 mpz_init(test);
 mpq_init(testfrac);
 mpz_set_str(test,"144", 10);
 mpq_set_str(testfrac,"-200/144",10);
 mpq_canonicalize(testfrac);
 printf("The 2-val of 144 is: %d\n",val_int(2,test));
 printf("The 3-val of 144 is: %d\n",val_int(3,test));
 printf("The 5-val of -200/144 is: %d\n",val(5,testfrac));
 printf("The 2-val of -200/144 is: %d\n",val(2,testfrac));
 mpz_clear(test);
 mpq_clear(testfrac);
 }*/
