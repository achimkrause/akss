#include "test.h"

/*int* parse_array(const char *input) {
  char *input_cpy = malloc(strlen(input)+1);
  strcpy(input_cpy, input);
  int length =1;
  char *in_ptr = input_cpy;
  char *s;
  
  while(*in_ptr != '\0') {

    if(*in_ptr == ',') {
      length++;
    }
    in_ptr++;
  }
  int* result = malloc(sizeof(int)*length);
  s=strtok(input_cpy, "{}, ");
  for(int i=0; i<length; i++) {
    sscanf(s,"%u",result+i);
    s=strtok(NULL, "{}, ");
  }
  return result;

}

int parse_matrix(const char *input, Matrix *mat) {
  char *input_cpy = malloc(strlen(input)+1);
  strcpy(input_cpy, input);
  int height =1;
  int width =1;
  char *in_ptr = input_cpy;
  char *s;

  while(*in_ptr != '}') {

    if(*in_ptr == ',') {
      width++;
    }
    in_ptr++;
  }
  while(*in_ptr != '\0') {
    if(*in_ptr == '{') {
      height++;
    }
    in_ptr++;
  }
  matrix_init(mat, height, width);
  s=strtok(input_cpy, "{}, ");

  for(int i=0; i<height*width; i++) {
    mpq_set_str(mat->entries[i], s, 10);

    s=strtok(NULL, "{}, ");
  }
  free(input_cpy);
  return 1;
}*/







void pprint_matrix(Matrix *mat) {
  mpq_t *entry = mat->entries;
  for (int i = 0; i < mat->height; i++) {
    for (int j = 0; j < mat->width; j++) {
      gmp_printf("%Qd ", *entry);
      entry++;
    }
    printf("\n");
  }
  printf("\n");
}



























/*void test_smith()
 {
 Matrix *mat;
 if(!(mat=init(4,4)))
 {
 printf("Allocation fault");
 return;
 }
 const int vals[] = {10,24,3,4,5,6,7,8,9,10,11,12,13,14,15,17};
 fill_matrix(mat, vals);

 Matrix *a;
 Matrix *b;
 Matrix *a2;
 Matrix *b2;
 a=init(4,4);
 b=init(4,4);
 a2=init(4,4);
 b2=init(4,4);
 set_unit(a);
 set_unit(a2);
 set_unit(b);
 set_unit(b2);
 printf("all set\n");
 smith(2,mat,a,b,a2,b2);
 pprint_mat(mat);
 pprint_mat(a);

 pprint_mat(b);
 pprint_mat(a2);
 pprint_mat(b2);

 }*/
/*
void pprint_abelian(int p, AbelianGroup arr) {
  int *entry = arr.orders;
  bool rep = false;
  mpz_t pow;
  mpz_init(pow);
  for (int i = 0; i < arr.tor_rank; i++) {
    if (rep) {
      printf("+");
    } else {
      rep = true;
    }
    mpz_ui_pow_ui(pow, p, *entry);
    gmp_printf("Z/%Zd", pow);
  }
  for (int i = 0; i < arr.free_rank; i++) {
    if (rep) {
      printf("+");
    } else {
      rep = true;
    }
    printf("Z");
  }
  printf("\n");
  mpz_clear(pow);
}


void fill_matrix(int *vals, Matrix *mat) {
  int *entry_vals = vals;
  mpq_t *entry_mat = mat->entries;

  for (int i = 0; i < mat->width * mat->height; i++) {
    mpq_set_ui(*entry_mat, *entry_vals, 1);
    entry_vals++;
    entry_mat++;
  }
}

void test_kernel(int p, Matrix *f, AbelianGroup x, AbelianGroup y,
    Matrix *g) {
  printf("given abelian groups X,Y:\n");
  pprint_abelian(p, x);
  pprint_abelian(p, y);
  printf("and a matrix f: X->Y: \n");
  pprint_matrix(f);

  MatrixArray from_X;
  matrix_arr_init(&from_X, 1);
  (*(from_X.entries)) = matrix_init(f->width, f->width);
  set_unit_range(0, 0, f->width, f->width, *from_X.entries);

  MatrixArray to_X;
  if (g != NULL) {
    matrix_arr_init(&to_X, 1);
    *(to_X.entries) = g;

    printf("and a matrix g to X:\n");
    pprint_matrix(*(to_X.entries));
  } else {
    matrix_arr_init(&to_X, 0);
  }

  AbelianGroup k;
  MatrixArray to_K;
  MatrixArray from_K;

  kernel(p, f, x, y, to_X, from_X, &k, &to_K, &from_K);

  printf("the kernel of f is given by:\n");
  pprint_abelian(p, k);
  printf("the inclusion of it into X is given by:\n");
  pprint_matrix(*(from_K.entries));
  if (g != NULL) {
    printf("and a lift of g through the kernel is given by:\n");
    pprint_matrix(*(to_K.entries));

  }
  matrix_arr_clear(from_X);
  free(to_X.entries);
  abelian_clear(&k);
  matrix_arr_clear(from_K);
  matrix_arr_clear(to_K);
}

void test_cokernel(int p, Matrix *f, AbelianGroup x, AbelianGroup y,
    Matrix *g) {
  printf("given abelian groups X,Y:\n");
  pprint_abelian(p, x);
  pprint_abelian(p, y);
  printf("and a matrix f: X->Y: \n");
  pprint_matrix(f);

  MatrixArray to_Y;
  matrix_arr_init(&to_Y, 1);
  *(to_Y.entries) = matrix_init(f->height, f->height);
  set_unit_range(0, 0, f->height, f->height, *to_Y.entries);

  MatrixArray from_Y;
  if (g != NULL) {
    matrix_arr_init(&from_Y, 1);
    *(from_Y.entries) = g;

    printf("and a matrix g from Y:\n");
    pprint_matrix(*(from_Y.entries));
  } else {
    matrix_arr_init(&from_Y, 0);
  }

  AbelianGroup c;
  MatrixArray to_C;
  MatrixArray from_C;

  cokernel(p, f, y, to_Y, from_Y, &c, &to_C, &from_C);

  printf("the cokernel of f is given by:\n");
  pprint_abelian(p, c);
  printf("the projection from Y to it is given by:\n");
  pprint_matrix(*(to_C.entries));
  if (g != NULL) {
    printf("and a factorization of g through it is given by:\n");
    pprint_matrix(*(from_C.entries));
  }

  matrix_arr_clear(to_Y);
  free(from_Y.entries);
  abelian_clear(&c);
  matrix_arr_clear(to_C);
  matrix_arr_clear(from_C);
}

void test_epi_mono(int p, Matrix *f, AbelianGroup x, AbelianGroup y) {
  printf("given abelian groups X,Y:\n");
  pprint_abelian(p, x);
  pprint_abelian(p, y);
  printf("and a matrix f: X->Y: \n");
  pprint_matrix(f);

  Matrix *inc;
  Matrix *proj;
  AbelianGroup img;
  epi_mono(p, f, x, y, &img, &proj, &inc);

  printf("the image is:\n");
  pprint_abelian(p, img);

  printf("the projection matrix is:\n");
  pprint_matrix(proj);

  printf("the inclusion matrix is:\n");
  pprint_matrix(inc);

  matrix_clear(proj);
  matrix_clear(inc);
  abelian_clear(&img);
}
*/