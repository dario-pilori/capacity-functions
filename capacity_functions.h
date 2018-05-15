/*
* Header for capacity_functions.h
*/
// Includes
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

// Definitions
#define N_GH 10

// External global variables
extern const double x[N_GH];
extern const double w[N_GH];

// Evaluate symbol energy for PAM and QAM
double symbol_energy(const double *C, const double *Pk, int M);
double complex_symbol_energy(const double complex *C, const double *Pk, int M);

// Helper function: inserts a zero inside integer i at position k. i has at most nb bits.
unsigned int insert_zero(unsigned int i, unsigned int k, unsigned int nb);

// Helper function to aid evaluation of sum
void maxxx(double *a, double b);

// Evaluation of AWGN and BMD(GMI) for PAM
double pam_eval_mi(const double *C, int M, double s, const double *Pk);
double pam_eval_gmi(const double *C, int M, double s, const double *Pk);

// Evaluation of AWGN and BMD(GMI) for QAM
double qam_eval_mi(const double complex *C, int M, double s, const double *Pk);
double qam_eval_gmi(const double complex *C, int M, double s, const double *Pk);

// Evaluation of log-likelihood ratios for PAM and QAM
void qam_soft_decode(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, const double *s2, double *l);
void pam_soft_decode(const double *y, int Ns, const double *C, 
        const double *Pk, int M, const double *s2, double *l);

// Evaluation of log-likelihood ratios for QAM with phase noise
void qam_soft_decode_pn(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, double Kn, double Kp, double B0, double *l);
void qam_soft_decode_pn_maxlog(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, double Kn, double Kp, double *l);

// Monte-Carlo evaluation of AWGN MI
double qam_montecarlo_mi(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, double s2);
double pam_montecarlo_mi(const double *y, int Ns, const double *C,
        const double *Pk, int M, double s2);