/*
* Header for gausshermite_functions.h
*/
// Includes
#include <math.h>
#include <complex.h>

// Definitions
#define N_GH 10

// External global variables
extern const double x[N_GH];
extern const double w[N_GH];

// Functions
double symbol_energy(const double *C, const double *Pk, int M);
double complex_symbol_energy(const double complex *C, int M);
unsigned int insert_zero(unsigned int i, unsigned int k, unsigned int nb);
double pam_eval_mi(const double *C, int M, double s, const double *Pk);
double pam_eval_gmi(const double *C, int M, double s);
double qam_eval_mi(const double complex *C, int M, double s);
double qam_eval_gmi(const double complex *C, int M, double s);
