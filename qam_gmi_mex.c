/*
 * qam_gmi_mex.c - Compute MI and GMI for QAM
 * 
 * Usage: qam_gmi_mex(C, SNR, Pk)
 * C    :=   Complex constellation in Gray-mapping order
 * SNR  :=   Vector of SNR (Es/No, dB)
 * Pk   :=   Probability of each constellation symbol
 *
 * Use this function to compute MI and GMI for M-QAM
 * constellations in an AWGN channel using Gauss-Hermite quadrature.
 *
 * Compile with: mex -lm -R2018a CFLAGS="$CFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" qam_gmi_mex.c
 * (requires MATLAB R2018a or newer versions)
 * Works under 64-bit Linux. Don't know/care under other OS.
 *
 * 2018 - Dario Pilori <dario.pilori@polito.it>
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
    size_t M, N;                   /* constellation and SNR size */
    int i;                         /* loop index */
    double complex *C;
    double *mi, *gmi, *snr, *Pk;
    double Es, sigma;
    
    /* Verify input */
    if(nrhs != 3) {
        mexErrMsgIdAndTxt("DspLibrary:qam_gmi_mex:nrhs",
                "Three inputs required.");
    }
    
    /* Verify output */
    if(nlhs > 2) {
        mexErrMsgIdAndTxt("DspLibrary:qam_gmi_mex:nlhs",
                "Max two outputs.");
    }
    
    /* Get sizes */
    M = mxGetM(prhs[0]);
    N = mxGetM(prhs[1]);
        
    /* Get pointers */
    snr = mxGetDoubles(prhs[1]);
    C = (double complex *) mxGetComplexDoubles(prhs[0]);
    Pk = mxGetDoubles(prhs[2]);
    
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix((mwSize)N,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)N,1,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    mi = mxGetDoubles(plhs[0]);
    gmi = mxGetDoubles(plhs[1]);
    
    /* Calculate symbol energy */
    Es = complex_symbol_energy(C, Pk, M);
    
    /* Call functions */
    #pragma omp parallel for private(sigma)
    for(i=0; i<N; i++) {
        sigma = sqrt(Es) * pow(10,-snr[i]/20);
        mi[i] = qam_eval_mi(C, M, sigma, Pk);
        gmi[i] = qam_eval_gmi(C, M, sigma, Pk);
    }
}

