#ifndef TEST_H
#define TEST_H

#include "smith.h"
#include "matrix.h"
#include "abelian.h"
#include <stdio.h>
#include <stdbool.h>

//void test_smith(void);

void pprint_abelian(int p, AbelianGroup arr);
void pprint_matrix(Matrix *mat);

void fill_matrix(int *vals, Matrix *mat);

void test_kernel(int p, Matrix *f, AbelianGroup x, AbelianGroup y,
    Matrix *g);

void test_cokernel(int p, Matrix *f, AbelianGroup x, AbelianGroup y,
    Matrix *g);

void test_epi_mono(int p, Matrix *f, AbelianGroup x, AbelianGroup y);

#endif
