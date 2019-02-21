/*
 * qam_llr_maxlog_mex.c - Compute LLRs for QAM
 *
 * Usage: l = qam_llr_maxlog_mex(C, sigma2, y, Pk)
 * C     :=   Complex constellation in Gray-mapping order
 * sigma2:=   Noise variance per each constellation point
 * y     :=   Received complex symbols
 * Pk    :=   Probability of each constellation symbol
 *
 * Use this function to compute log-likelihood-ratios
 * for M-QAM constellations assuming an AWGN channel.
 * The output of this function is comptabible to MATLAB's qamdemod()
 * function of the Communication System Toolbox.
 * The max-log approximation is applied.
 *
 * Compile with: mex -lm -R2018a qam_llr_maxlog_mex.c
 * (requires MATLAB R2018a or newer versions)
 * Designed for 64-bit Linux.
 *
 * Copyright (c) 2019 - Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
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
    double *l, *Pk, *s2;                 
    
    /* Verify input */
    if(nrhs != 4) {
        mexErrMsgIdAndTxt("DspLibrary:qam_llr_maxlog_mex:nrhs",
                "Four inputs required.");
    }
    
    /* Verify output */
    if(nlhs > 1) {
        mexErrMsgIdAndTxt("DspLibrary:qam_llr_maxlog_mex:nlhs",
                "Max one output.");
    }
    
    /* Get sizes */
    M = mxGetM(prhs[0]);
    Ns = mxGetM(prhs[2]);
    
    /* Verify sizes */
    if ((mxGetM(prhs[1])!=M)||(mxGetM(prhs[3])!=M)) {
        mexErrMsgIdAndTxt("DspLibrary:qam_llr_maxlog_mex:sizes",
                "Probabilities and variances must have same size as constellation.");
    }
    
    /* Get noise variance */
    s2 = mxGetDoubles(prhs[1]);
        
    /* Get constellation and received data */
    C = (double complex *) mxGetComplexDoubles(prhs[0]);
    y = (double complex *) mxGetComplexDoubles(prhs[2]);
    
    /* Get probabilities */
    Pk = mxGetDoubles(prhs[3]);
    
    /* Allocate the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize) (Ns*log2(M)),1,mxREAL);
    
    /* get a pointer to the output matrix */
    l = mxGetDoubles(plhs[0]);
        
    /* Call function */
    qam_soft_decode_maxlog(y, Ns, C, Pk, M, s2, l);
}
