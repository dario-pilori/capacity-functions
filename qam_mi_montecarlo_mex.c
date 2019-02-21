/*
 * qam_mi_montecarlo_mex.c - Compute MI for QAM using MonteCarlo for AWGN
 *
 * Usage: mi = qam_mi_montecarlo_mex(C, sigma2, a-y, Pk)
 * C     :=   Complex constellation in Gray-mapping order
 * sigma2:=   Noise variance (average)
 * a-y   :=   Received noise (transmitted-received symbols)
 * Pk    :=   Probability of each constellation symbol
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
#include <stdio.h>
#include "capacity_functions.h"
#include "mex.h"

/* Gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    size_t M, Ns;                   /* constellation and data size */
    double complex *C, *y;          /* Data and constellation */
    double *Pk;                 
    double s2, mi;
    
    /* Verify input */
    if(nrhs != 4) {
        mexErrMsgIdAndTxt("DspLibrary:qam_mi_montecarlo:nrhs",
                "Four inputs required.");
    }
    
    /* Verify output */
    if(nlhs > 1) {
        mexErrMsgIdAndTxt("DspLibrary:qam_mi_montecarlo:nlhs",
                "Max one output.");
    }
    
    /* Get sizes */
    M = mxGetM(prhs[0]);
    Ns = mxGetM(prhs[2]);
    
    /* Get noise variance */
    s2 = mxGetScalar(prhs[1]);
        
    /* Get constellation and received data */
    C = (double complex *) mxGetComplexDoubles(prhs[0]);
    y = (double complex *) mxGetComplexDoubles(prhs[2]);
    
    /* Get probabilities */
    Pk = mxGetDoubles(prhs[3]);
        
    /* Call function */
    mi = qam_montecarlo_mi(y, Ns, C, Pk, M, s2);
    
    /* get a pointer to the output matrix */
    plhs[0] = mxCreateDoubleScalar(mi);
}
