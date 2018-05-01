/*
 * qam_llr_mex.c - Compute LLRs for QAM
 *
 * Usage: qam_llr_mex(C, SNR, y, Pk)
 * C    :=   Complex constellation in Gray-mapping order
 * SNR  :=   Vector of SNR (Es/No, dB)
 * y    :=   Received complex symbols
 * Pk   :=   Probability of each constellation symbol
 *
 * Use this function to compute log-likelihood-ratios
 * for M-QAM constellations.
 *
 * Compile with: mex -lm -R2018a qam_llr_mex.c
 * (requires MATLAB R2018a or newer versions)
 * Works under 64-bit Linux. Don't know/care under other OSs.
 *
 * 2018 - Dario Pilori <dario.pilori@polito.it>
 */
#include <math.h>
#include <complex.h>
#include "capacity_functions.h"
#include <omp.h>
#include "mex.h"

/* Gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    size_t M, Ns;                   /* constellation and data size */
    double complex *C, *y;          /* Data and constellation */
    double *l, *Pk;                 
    double snr, Es, sigma;
    
    /* Verify input */
    if(nrhs != 4) {
        mexErrMsgIdAndTxt("DspLibrary:qam_gmi_mex:nrhs",
                "Four inputs required.");
    }
    
    /* Verify output */
    if(nlhs > 1) {
        mexErrMsgIdAndTxt("DspLibrary:qam_gmi_mex:nlhs",
                "Max one output.");
    }
    
    /* Get sizes */
    M = mxGetM(prhs[0]);
    Ns = mxGetM(prhs[2]);
    
    /* Get SNR */
    snr = mxGetScalar(prhs[1]);
        
    /* Get constellation, received data and probabilities */
    C = (double complex *) mxGetComplexDoubles(prhs[0]);
    y = (double complex *) mxGetComplexDoubles(prhs[2]);
    Pk = mxGetDoubles(prhs[3]);
    
    /* Allocate the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize) (Ns*log2(M)),1,mxREAL);
    
    /* get a pointer to the output matrix */
    l = mxGetDoubles(plhs[0]);
    
    /* Calculate symbol energy */
    Es = complex_symbol_energy(C, Pk, M);
    
    /* Call function */
    sigma = sqrt(Es) * pow(10,-snr/20);
    qam_soft_decode(y, Ns, C, Pk, M, sigma, l);
}

