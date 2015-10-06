#ifndef SPECTRAL_H
#define SPECTRAL_H

#include "abelian.h"
#include "matrix.h"

typedef struct {
	AbelianGroup **groups;
	Matrix **maps;
	int length;
	int current;

} GroupSequence;

typedef struct {
	GroupSequence **entries;
	int length;
} GroupSequences1D;

typedef struct {
	GroupSequences1D **entries;
	int length;
} GroupSequences2D;

typedef struct {
	GroupSequences2D **entries;
	int length;
} GroupSequences3D;

typedef struct {
  MatrixArray **entries;
  int length;
} MatrixArray2D;

typedef struct {
  MatrixArray2D **entries;
  int length;
} MatrixArray3D;

typedef struct {
  MatrixArray3D **entries;
  int length;
} MatrixArray4D;

typedef struct {
	GroupSequences3D *kernels;
	GroupSequences3D *cokernels;
	MatrixArray4D *differentials;
} SpectralSequence;

void prep_index_group_seq_3D(GroupSequences3D *gr, const int p, const int q, const int s);
void prep_index_group_seq_2D(GroupSequences2D *gr, const int q, const int s);
void prep_index_group_seq_1D(GroupSequences1D *gr, const int s);

void prep_index_matrix_arr4D(MatrixArray4D *arr, const int p,const int q,const int s,const int r);
void prep_index_matrix_arr3D(MatrixArray3D *arr, const int q,const int s,const int r);
void prep_index_matrix_arr2D(MatrixArray2D *arr, const int s,const int r);
void prep_index_matrix_arr(MatrixArray *arr, const int r);

GroupSequence* getGroupSequence(GroupSequences3D *gr, const int p, const int q, const int s);
MatrixArray* getMatrixArray(MatrixArray4D *arr, const int p, const int q, const int s);

int setE2(SpectralSequence *ss, AbelianGroup *ab, const int p, const int q, const int s);
int setDiff(const int prime, SpectralSequence *ss, Matrix *d, const int r, const int p, const int q, const int s);

int resetE2(SpectralSequence *ss, const int p, const int q, const int s);
int resetDiff(SpectralSequence *ss, const int r, const int p, const int q, const int s);

void group_sequence_append(GroupSequence *seq, Matrix *mat, AbelianGroup *ab); 
void group_sequence_diff(const int prime, const Matrix *d, GroupSequence *kernels, GroupSequence *cokernels);




#endif