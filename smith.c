#include "smith.h"
#include "test.h" //only for debug

////p-valuation helper functions////////////////////////////////////////////////////////////////////////7
int valuation_int(const int p, const mpz_t z)
{
        //fprintf(stderr, "valuation_int called on %u, %u\n", p, mpz_get_ui(z));
        mpz_t rem;
        mpz_init(rem);
        int result = 0;
        mpz_set(rem, z);
        while(mpz_divisible_ui_p(rem, p))
        {
                result++;
                mpz_divexact_ui(rem, rem, p);
        }
        mpz_clear(rem);
        return result;
}


int valuation(const int p, const mpq_t lam)
{
        mpz_t z;
        mpz_init(z);
        mpz_set(z, mpq_numref(lam));
        int result;
        result = valuation_int(p, z);
        if(result>0)
        {
                mpz_clear(z);
                return result;
        }
        else
        {
                mpz_set(z, mpq_denref(lam));
                result = -valuation_int(p,z);
                mpz_clear(z);
                return result;
        }
}


/////////////////////////////////////////////////////smith algorithm////////////////////////////////////////////////////
void smith(int p, struct matrix *mat, struct matrix_arr to_X, struct matrix_arr from_X_0,struct matrix_arr to_Y_0, struct matrix_arr from_Y)
{
        //mat represents a map X -> Y. to_X, from_X, to_Y, from_Y are representatives for maps to/from X/Y. 
        //smith brings mat into smith normal form by applying base changes to X and Y and transforming to_X and to_Y accordingly.
        mpq_t *entry;
        mpq_t min;
        mpq_init(min);
        mpq_t lam;
        mpq_init(lam);
        
        struct matrix_arr from_X;
        struct matrix_arr to_Y;

        from_X.length = from_X_0.length + 1;
        to_Y.length = to_Y_0.length + 1;
        
        from_X.entries = (struct matrix**) malloc(sizeof(struct matrix*)*from_X.length);
        to_Y.entries   = (struct matrix**) malloc(sizeof(struct matrix*)*to_Y.length);

        struct matrix **entry_from_X = from_X.entries;
        struct matrix **entry_from_X_0 = from_X_0.entries;

        for(int i=0; i<from_X_0.length; i++)
        {
                *entry_from_X = *entry_from_X_0;
                entry_from_X++;
                entry_from_X_0++;
        }
        *entry_from_X = mat;

        struct matrix **entry_to_Y = to_Y.entries;
        struct matrix **entry_to_Y_0 = to_Y_0.entries;

        for(int i=0; i<to_Y_0.length; i++)
        {
                *entry_to_Y = *entry_to_Y_0;
                entry_to_Y++;
                entry_to_Y_0++;
        }
        *entry_to_Y = mat;

        int block=0;

        while(block < mat->width && block < mat->height)
        {
                entry = mat->entries + block + block*(mat->width);    
                //find smallest in block:
                int i_min=-1;
                int j_min=-1;
                int min_val=0;
                for(int i=block; i<mat->height; i++)
                {
                        for(int j=block; j<mat->width; j++)
                        {
                                if(mpq_sgn(*entry)!=0)
                                {
                                        if(i_min==-1)
                                        {
                                                i_min=i;
                                                j_min=j;
                                                mpq_set(min, *entry);
                                                min_val = valuation(p, min);
                                        }
                                        else
                                        {
                                                if(valuation(p, *entry)<min_val)
                                                {
                                                        i_min=i;
                                                        j_min=j;
                                                        mpq_set(min, *entry);
                                                        min_val = valuation(p, min);
                                                }
                                        }

                                }
                                entry++;
                        }
                        entry += block;
                } 
                if(i_min==-1)
                {
                        mpq_clear(min);
                        mpq_clear(lam);
                        free(from_X.entries);
                        free(to_Y.entries);
                        return;
                } 

                entry = mat->entries + j_min + block*(mat->width);
                for(int i=block; i<mat->height; i++)
                {
                        if(i!=i_min)
                        {
                                mpq_div(lam, *entry, min);
                                add_base(to_Y, from_Y, i, i_min, lam);
                                //has the effect of subtracting lam times the i_min-row from the i-row.
                        }
                        entry += mat->width;
                }  
                
                entry = mat->entries + i_min*(mat->width) + block;
                for(int j=block; j<mat->width; j++)
                {
                        if(j!=j_min)
                        {
                                mpq_div(lam, *entry, min);
                                mpq_neg(lam, lam);
                                add_base(to_X, from_X, j_min, j, lam);
                                //has the effect of adding lam times the j_min col to the j-col.
                        }
                        entry++;
                }
                swap_base(to_Y,from_Y,i_min,block);
                swap_base(to_X,from_X,j_min,block);

                mpz_ui_pow_ui(mpq_numref(lam), p, min_val);
                mpz_set_ui(mpq_denref(lam), 1);
                mpq_div(lam, lam, min);
                mul_base(to_X,from_X, block, lam);

                block++;                 
        }
        mpq_clear(min);
        mpq_clear(lam);
        free(from_X.entries);
        free(to_Y.entries);

}

