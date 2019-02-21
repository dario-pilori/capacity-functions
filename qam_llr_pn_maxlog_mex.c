/*
 * qam_llr_pn_maxlog_mex.c - Compute LLRs for QAM
 *
 * Usage: l = qam_llr_pn_maxlog_mex(C, Kn, Kp, y, Pk)
 * C     :=   Complex constellation in Gray-mapping order
 * Kn    :=   1/Noise variance
 * Kp    :=   ~ 1/Phase noise variance
 * y     :=   Received complex symbols
 * Pk    :=   Probability of each constellation symbol
 *
 * Use this function to compute log-likelihood-ratios
 * for M-QAM constellations assuming an AWGN channel
 * with phase noise (10.1109/TWC.2014.040714.130731). 
 * The output of this function is comptabible to MATLAB's qamdemod()
 * function of the Communication System Toolbox.
 * This function uses the max-log approximation.
 *
 * Compile with: mex -lm -R2018a qam_llr_maxlog_mex.c
 * (requires MATLAB R2018a or newer versions)
 * Designed for 64-bit Linux.
 *
 * Copyright (c) 2018 - Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
 * MIT License
 */
#include <complex.h>
#include <omp.h>
#include <math.h>
#include "capacity_functions.h"
#include "mex.h"

/* Gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    size_t M, Ns;                   /* constellation and data size */
    double complex *C, *y;          /* Data and constellation */
    double *l, *Pk;                 
    double Kn, Kp;
    
    /* Verify input */
    if(nrhs != 5) {
        mexErrMsgIdAndTxt("DspLibrary:qam_llr_maxlog_mex:nrhs",
                "Five inputs required.");
    }
    
    /* Verify output */
    if(nlhs > 1) {
        mexErrMsgIdAndTxt("DspLibrary:qam_llr_maxlog_mex:nlhs",
                "Max one output.");
    }
    
    /* Get sizes */
    M = mxGetM(prhs[0]);
    Ns = mxGetM(prhs[3]);
    
    /* Verify sizes */
    if (mxGetM(prhs[4])!=M) {
        mexErrMsgIdAndTxt("DspLibrary:qam_llr_maxlog_mex:sizes",
                "Probabilities must have same size as constellation.");
    }
    
    /* Get noise variances */
    Kn = mxGetScalar(prhs[1]);
    Kp = mxGetScalar(prhs[2]);
        
    /* Get constellation and received data */
    C = (double complex *) mxGetComplexDoubles(prhs[0]);
    y = (double complex *) mxGetComplexDoubles(prhs[3]);
    
    /* Get probabilities */
    Pk = mxGetDoubles(prhs[4]);
    
    /* Allocate the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize) (Ns*log2(M)),1,mxREAL);
    
    /* get a pointer to the output matrix */
    l = mxGetDoubles(plhs[0]);
        
    /* Call function */
    qam_soft_decode_pn_maxlog(y, Ns, C, Pk, M, Kn, Kp, l);
}
