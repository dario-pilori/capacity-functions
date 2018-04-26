/*
 * qam_llr_mex.c - Compute LLRs for QAM
 *
 * Use this function to compute log-likelihood-ratios
 * for M-QAM constellations.
 *
 * Compile with: mex -lm -R2018a qam_llr_mex.c
 * (requires MATLAB R2018a or newer versions)
 *
 * 2018 - Dario Pilori <dario.pilori@polito.it>
 */
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "gausshermite_functions.h"
#include "mex.h"

/* Gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    size_t M, Ns;                   /* constellation and data size */
    double complex *C, *y;
    double *l;
    double snr, Es, sigma;
    
    /* Verify input */
    if(nrhs != 3) {
        mexErrMsgIdAndTxt("dsp-library:qam_gmi_mex:nrhs",
                "Three inputs required.");
    }
    
    /* Verify output */
    if(nlhs > 1) {
        mexErrMsgIdAndTxt("dsp-library:qam_gmi_mex:nlhs",
                "Max one output.");
    }
    
    /* Get sizes */
    M = mxGetM(prhs[0]);
    Ns = mxGetM(prhs[2]);
    
    /* Get SNR */
    snr = mxGetScalar(prhs[1]);
        
    /* Get constellation and received data */
    C = (double complex *) mxGetComplexDoubles(prhs[0]);
    y = (double complex *) mxGetComplexDoubles(prhs[2]);
    
    /* Allocate the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize) (Ns*log2(M)),1,mxREAL);
    
    /* get a pointer to the output matrix */
    l = mxGetDoubles(plhs[0]);
    
    /* Calculate symbol energy */
    Es = complex_symbol_energy(C, M);
    
    /* Call function */
    sigma = sqrt(Es) * pow(10,-snr/20);
    qam_soft_decode(y, Ns, C, M, sigma, l);
}

