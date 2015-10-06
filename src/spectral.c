#include "spectral.h"

void prep_index_group_seq_3D(GroupSequences3D *gr, const int p, const int q, const int s) {
	//makes sure gr->entries is large enough to be adressed with index p.
	if(gr->length <= p) {

		gr->entries = realloc(gr->entries, sizeof(GroupSequences2D*)*2*(p+1));

		for(int k=gr->length; k<2*(p+1); k++)
		{
			gr->entries[k] = malloc(sizeof(GroupSequences2D));
			gr->entries[k]->length = 0;
			gr->entries[k]->entries = NULL;
		}
		gr->length = 2*(p+1);
		
	}
	prep_index_group_seq_2D(gr->entries[p], q, s);
}

void prep_index_group_seq_2D(GroupSequences2D *gr, const int q, const int s) {
	//makes sure gr->entries is large enough to be adressed with index q.
	if(gr->length <= q) {
		gr->entries = realloc(gr->entries, sizeof(GroupSequences1D*)*2*(q+1));

		for(int k=gr->length; k<2*(q+1); k++)
		{
			gr->entries[k] = malloc(sizeof(GroupSequences1D));
			gr->entries[k]->length = 0;
			gr->entries[k]->entries = NULL;
		}
		gr->length = 2*(q+1);
		
	}
	prep_index_group_seq_1D(gr->entries[q], s);
}

void prep_index_group_seq_1D(GroupSequences1D *gr, const int s) {
	//makes sure gr->entries is large enough to be adressed with index s.
	if(gr->length <= s) {
		gr->entries = realloc(gr->entries, sizeof(GroupSequence*)*2*(s+1));

		for(int k=gr->length; k<2*(s+1); k++)
		{
			gr->entries[k] = malloc(sizeof(GroupSequence));
			gr->entries[k]->length = 0;
			gr->entries[k]->current = -1;
			gr->entries[k]->maps = NULL;
			gr->entries[k]->groups = NULL;
		}
		gr->length = 2*(s+1);
	}
}

void prep_index_matrix_arr4D(MatrixArray4D *arr, const int p,const int q,const int s,const int r) {
	if(arr->length <= p) {
		arr->entries = realloc(arr->entries, sizeof(MatrixArray3D*)*2*(p+1));
			
		for(int k=arr->length; k<2*(p+1); k++)
		{
			arr->entries[k] = malloc(sizeof(MatrixArray3D));
			arr->entries[k]->length = 0;
			arr->entries[k]->entries = NULL;
		}
		arr->length = 2*(p+1);
	}

	prep_index_matrix_arr3D(arr->entries[p],q,s,r);
}

void prep_index_matrix_arr3D(MatrixArray3D *arr, const int q, const int s, const int r) {
	if(arr->length <= q) {
		arr->entries = realloc(arr->entries, sizeof(MatrixArray2D*)*2*(q+1));

		for(int k=arr->length; k<2*(q+1); k++)
		{
			arr->entries[k] = malloc(sizeof(MatrixArray2D));
			arr->entries[k]->length = 0;
			arr->entries[k]->entries = NULL;
		}
		arr->length = 2*(q+1);
		
	}
	prep_index_matrix_arr2D(arr->entries[q],s,r);
}

void prep_index_matrix_arr2D(MatrixArray2D *arr, const int s, const int r) {
	if(arr->length <= s) {
		MatrixArray **temp = arr->entries;
		arr->entries = realloc(arr->entries, sizeof(MatrixArray*)*2*(s+1));
		
		for(int k=0; k<arr->length; k++)
		{
			arr->entries[k] = temp[k];
		}
		if(temp != NULL) {
			free(temp);	
		}

		for(int k=arr->length; k<2*(s+1); k++)
		{
			arr->entries[k] = malloc(sizeof(MatrixArray));
			arr->entries[k]->length = 0;
			arr->entries[k]->entries = NULL;
		}
		arr->length = 2*(s+1);
	}

	prep_index_matrix_arr(arr->entries[s], r);
}

void prep_index_matrix_arr(MatrixArray *arr, int r) {
	if(arr->length <= r) {
		Matrix **temp = arr->entries;
		arr->entries = malloc(sizeof(Matrix*)*2*(r+1));
		
		for(int k=0; k<arr->length; k++)
		{
			arr->entries[k] = temp[k];
		}
		if(temp != NULL) {
			free(temp);	
		}

		arr->length = 2*(r+1);
	}
}

GroupSequence* getGroupSequence(GroupSequences3D *gr, const int p, const int q, const int s){
	return gr->entries[p]->entries[q]->entries[s];
}

MatrixArray* getMatrixArray(MatrixArray4D *arr, const int p, const int q, const int s){
	return arr->entries[p]->entries[q]->entries[s];
}

int setE2(SpectralSequence *ss, AbelianGroup *ab, const int p, const int q, const int s) {
	//if necessary, extend array to (p,q,s).
	//then check whether E2(p,q,s) is already set. If no, extend coker array as well and set the corresponding E2.
	prep_index_group_seq_3D(ss->kernels, p,q ,s);
	prep_index_group_seq_3D(ss->cokernels, p,q ,s);
	
	GroupSequence *ker = getGroupSequence(ss->kernels, p, q, s);
	if(ker->current == -1) {

		GroupSequence *coker = getGroupSequence(ss->kernels, p, q, s);

		Matrix *id = matrix_alloc();
		matrix_init(id, ab->tor_rank + ab->free_rank, ab->tor_rank + ab->free_rank);
		set_unit(id);
		group_sequence_append(ker, id, ab);
		group_sequence_append(coker, id, ab);
		return 1;
	}
	else {
		//if E2 is already set there, return 0.
		return 0;
	}

	
}



int setDiff(const int prime, SpectralSequence *ss, Matrix *d, const int r, const int p, const int q, const int s) {
	//if either domain or codomain are not at current=r, return 0.
	prep_index_group_seq_3D(ss->kernels, p,q ,s);
	//differentials lower p-degree by r, increases q-degree by 2r-1, increases s-degree by 1.
	// (here p and r are double the usual p and r)
	if(p<r) {
		fprintf(stderr, "Differential too long.\n");
		return 0;
	}
	prep_index_group_seq_3D(ss->cokernels, p-r,q+2*r-1,s+1);
	GroupSequence *ker = getGroupSequence(ss->kernels, p, q, s);
	GroupSequence *coker = getGroupSequence(ss->cokernels, p-r, q+2*r-1, s+1);
	if(ker->current != r || coker->current != r) {
		fprintf(stderr, "Wrong r: ker is at %u, coker at %u.\n", ker->current, coker->current);
	}

	group_sequence_diff(prime, d, ker, coker);
	prep_index_matrix_arr4D(ss->differentials,p,q,s,r);
	MatrixArray *arr = getMatrixArray(ss->differentials, p,q,s);
	arr->entries[r] = d;
	return 1;
}

int resetE2(SpectralSequence *ss, const int p, const int q, const int s) {
	//if E2 is not set there, return 0.
	//otherwise, remove the group lists at (p,q,s) and all the differentials affected.
	return 0;
}
int resetDiff(SpectralSequence *ss, const int r, const int p, const int q, const int s) {
	//if, at (p,q,s), current<r, return 0.
	//otherwise, remove the differential at this position, and all later groups affected by it, and all differentials affected by them. 
	return 0;
}




void group_sequence_append(GroupSequence *seq, Matrix *mat, AbelianGroup *ab) {
	if(seq->current == seq->length - 1) {
		Matrix **temp_mat = seq->maps;
		AbelianGroup **temp_ab = seq->groups;

		seq->maps = malloc(2*(seq->length+1)*sizeof(Matrix*));
		seq->groups = malloc(2*(seq->length+1)*sizeof(AbelianGroup*));
		seq->length = 2*(seq->length+1);

		for(int i=0; i<=seq->current; i++)
		{
			seq->maps[i] = temp_mat[i];
			seq->groups[i] = temp_ab[i];
		}
		if(temp_mat != NULL)
		{
			free(temp_mat);
		}
		if(temp_ab != NULL)
		{
			free(temp_ab);	
		}
	}
	seq->current = seq->current+1;
	seq->maps[seq->current] = mat;
	seq->groups[seq->current] = ab;
}

void group_sequence_diff(const int p, const Matrix *d, GroupSequence *kernels, GroupSequence *cokernels) {
	
	MatrixArray *from_X = matrix_arr_alloc();
	matrix_arr_init(from_X, 1);
	from_X->entries[0] = kernels->maps[kernels->current];
	MatrixArray *from_K = matrix_arr_alloc();
	AbelianGroup *k = abelian_alloc();

	kernel(p, d, kernels->groups[kernels->current], 
		     cokernels->groups[cokernels->current], 
		     NULL, from_X, k, NULL, from_K);
	
	MatrixArray *to_Y = matrix_arr_alloc();
	matrix_arr_init(to_Y, 1);
	to_Y->entries[0] = cokernels->maps[cokernels->current];
	MatrixArray *to_C = matrix_arr_alloc();
	AbelianGroup *c = abelian_alloc();

	cokernel(p, d, cokernels->groups[cokernels->current],
	              to_Y, NULL, c, to_C, NULL); 

	
	group_sequence_append(kernels, from_K->entries[0], k);
	group_sequence_append(cokernels, to_C->entries[0], c);

    //only freeing these guys: all entries are still in the GroupSequences. Same for k and c. 
	free(from_X->entries);
	free(from_X);
	free(from_K->entries);
	free(from_K);
	free(to_C->entries);
	free(to_C);
	free(to_Y->entries);
	free(to_Y);
}
