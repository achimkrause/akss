#ifndef TEST_H
#define TEST_H

#include "smith.h"
#include "matrix.h"
#include "abelian.h"
#include <stdio.h>
#include <stdbool.h>



//void test_smith(void);



void pprint_abelian(int p, struct abelian arr);
void pprint_matrix(struct matrix *mat);

void fill_matrix(int *vals, struct matrix *mat);


void test_kernel(int p, struct matrix *f, struct abelian x, struct abelian y, struct matrix *g);

void test_cokernel(int p, struct matrix *f, struct abelian x, struct abelian y, struct matrix *g);

void test_epi_mono(int p, struct matrix *f, struct abelian x, struct abelian y);

#endif
