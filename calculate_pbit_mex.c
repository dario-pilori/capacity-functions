/*
 * calculate_pbit_mex.c - Compute bit-wise probabilities
 *
 * Usage: Pbit = calculate_pbit_mex(Pk)
 * Pk   :=   Probability of each constellation symbol
 *
 * Compile with: mex -R2018a -lm calculate_pbit_mex.c
 * Works under 64-bit Linux. Don't know/care under other OSs.
 *
 * 2018 - Dario Pilori <dario.pilori@polito.it>
 */
#include "mex.h"
#include <math.h>

/* Prototypes */
unsigned int insert_zero(unsigned int i, unsigned int k, unsigned int nb);
void calculate_bit_probabilities(const double *Pk, double *Pb, int m, int M);

/* Gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    size_t M, m;
    double *Pb, *Pk;
    
    /* Verify input */
    if(nrhs != 1) {
        mexErrMsgIdAndTxt("DspLibrary:qam_gmi_mex:nrhs",
                "One input required.");
    }
    
    /* Verify output */
    if(nlhs > 1) {
        mexErrMsgIdAndTxt("DspLibrary:qam_gmi_mex:nlhs",
                "Max one output.");
    }
    
    /* Get sizes */
    M = mxGetM(prhs[0]);
    m = log2(M);
    
    /* Get probabilities */
    Pk = mxGetDoubles(prhs[0]);
    
    /* Allocate the output matrix */
    plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL);
    Pb = mxGetDoubles(plhs[0]);
    
    /* Call function */
    calculate_bit_probabilities(Pk, Pb, m, M);
}

// Actual function
void calculate_bit_probabilities(const double *Pk, double *Pb, int m, int M)
{
    int j, k;
    unsigned int bj;
    
    for(k=0; k<m; k++)
    {
        Pb[k] = 0.0;
        for(j=0; j<M/2; j++)
        {
            bj = insert_zero(j, k, m);
            Pb[k] += Pk[bj];
        }
    }
}

// Helper function to insert a zero inside a number
unsigned int insert_zero(unsigned int i, unsigned int k, unsigned int nb)
{
  unsigned int b0, left, right;
  
  left = (i<<1) & (  ( (1<<(nb-k)) - 1)<<(k+1) );
  right = i & ((1<<k)-1);
  b0 = left | right;
  
  return b0;
}


