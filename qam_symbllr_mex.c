/*
 * qam_symbllr_mex.c - Compute LLRs for QAM
 *
 * Usage: l = qam_symbllr_mex(C, sigma2, y, Pk)
 * C     :=   Complex constellation in Gray-mapping order
 * sigma2:=   Noise variance
 * y     :=   Received complex symbols
 * Pk    :=   Probability of each constellation symbol
 *
 * Use this function to compute symbol-wise log-likelihood-ratios
 * for M-QAM constellations assuming an AWGN channel.
 *
 * Compile with: mex -lm -R2018a qam_llr_mex.c
 * (requires MATLAB R2018a or newer versions)
 * Designed for 64-bit Linux.
 *
 * Copyright (c) 2018 - Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
 * MIT License
 */
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "capacity_functions.h"
#include "mex.h"

/* Gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    size_t M, Ns;                   /* constellation and data size */
    double complex *C, *y;          /* Data and constellation */
    double *l, *Pk;                 /* LLR, probabilities */
    double s2;                      /* AWGN variance */
    
    /* Verify input */
    if(nrhs != 4) {
        mexErrMsgIdAndTxt("DspLibrary:qam_symbllr_mex:nrhs",
                "Four inputs required.");
    }
    
    /* Verify output */
    if(nlhs > 1) {
        mexErrMsgIdAndTxt("DspLibrary:qam_symbllr_mex:nlhs",
                "Max one output.");
    }
    
    /* Get sizes */
    M = mxGetM(prhs[0]);
    Ns = mxGetM(prhs[2]);
    
    /* Verify sizes */
    if (mxGetM(prhs[3])!=M) {
        mexErrMsgIdAndTxt("DspLibrary:qam_symbllr_mex:sizes",
                "Probabilities must have same size as constellation.");
    }
    
    /* Get noise variance */
    s2 = mxGetScalar(prhs[1]);
        
    /* Get constellation and received data */
    C = (double complex *) mxGetComplexDoubles(prhs[0]);
    y = (double complex *) mxGetComplexDoubles(prhs[2]);
    
    /* Get probabilities */
    Pk = mxGetDoubles(prhs[3]);
    
    /* Allocate the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize) (Ns*M),1,mxREAL);
    
    /* get a pointer to the output matrix */
    l = mxGetDoubles(plhs[0]);
        
    /* Call function */
    qam_symbol_decode(y, Ns, C, Pk, M, s2, l);
}
