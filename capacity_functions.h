/*
* Header for capacity_functions.h
*/
// Includes
#include <math.h>
#include <complex.h>

// Definitions
#define N_GH 10

// External global variables
extern const double x[N_GH];
extern const double w[N_GH];

// Evaluate symbol energy for PAM and QAM
double symbol_energy(const double *C, const double *Pk, int M);
double complex_symbol_energy(const double complex *C, int M);

// Helper function: inserts a zero inside integer i at position k. i has at most nb bits.
unsigned int insert_zero(unsigned int i, unsigned int k, unsigned int nb);

// Evaluation of AWGN and BMD(GMI) for PAM
double pam_eval_mi(const double *C, int M, double s, const double *Pk);
double pam_eval_gmi(const double *C, int M, double s);

// Evaluation of AWGN and BMD(GMI) for QAM
double qam_eval_mi(const double complex *C, int M, double s);
double qam_eval_gmi(const double complex *C, int M, double s);

// Evaluation of log-likelihood ratios for PAM and QAM
void qam_soft_decode(const double complex *y, int Ns, const double complex *C, int M, double s, double *l);
void pam_soft_decode(const double *y, int Ns, const double *C, int M, double s, double *l);
