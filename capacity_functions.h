/*! \file capacity_functions.h
 * 	\brief Information-theoretical functions
 *
 *  This file contains all the information-theoretical functions of this 
 *  project. These function can either be directly called with a small
 *  C program, or used with MATLAB's MEX library.
 *  
 *  Copyright (c) 2018 Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
 *  Released under MIT license.
*/

/* Include guard */
#ifndef CAPACITY_FUNCTIONS_H
#define CAPACITY_FUNCTIONS_H

/* Complex variables */
#include <complex.h>

/*! Number of points of Gauss-Hermite quadrature */
#define N_GH 10

/*! Roots of Hermite polynomial */
extern const double x[N_GH];
/*! Gauss-Hermite weights */
extern const double w[N_GH];

/*! Average real symbol energy
 * \param[in] C     Real constellation (e.g. PAM)
 * \param[in] Pk    Probabilities of each constellation point
 * \param[in] M     Number of points in the constellation
 */
double symbol_energy(const double *C, const double *Pk, int M);

/*! Average complex symbol energy
 * \param[in] C     Complex constellation (e.g. QAM)
 * \param[in] Pk    Probabilities of each constellation point
 * \param[in] M     Number of points in the constellation
 */
double complex_symbol_energy(const double complex *C, const double *Pk, int M);

/*! \brief Insert a zero
 *
 *  This function inserts a binary zero inside a positive binary number at a
 *  given position. This function is useful to generate the set of constellation
 *  points where a specific bit is either zero or one.
 * \param[in] i    Number over which zero is inserted 
 * \param[in] k    Position of zero insertion (0...Nb)
 * \param[in] nb   Number of bits of i
 */
unsigned int insert_zero(unsigned int i, unsigned int k, unsigned int nb);

/*! \brief Max star
 *
 *  This function helps the evaluation of log(exp(a)+exp(b)+exp(c)+...)
 *  in a numerically-stable way. This function evaluates log(exp(a)+exp(b))
 *  and inserts this value into a. Call this function multiple times
 *  to evaluate sum of more numbers.
 * \param[in,out] a    Number to sum and over which put the result 
 * \param[in]     b    Second number to sum
 */
void maxxx(double *a, double b);

/*! \brief Gauss-Hermite AWGN Mutual Information evaluation for PAM
 *
 *  This function evalautes the AWGN Mutual Information for a real-valued
 *  constellation (e.g. PAM)
 * \param[in] C     Real constellation (e.g. PAM)
 * \param[in] M     Number of points in the constellation
 * \param[in] s     Standard deviation of real-valued AWGN
 * \param[in] Pk    Probability of each constellation point
 */
double pam_eval_mi(const double *C, int M, double s, const double *Pk);

/*! \brief Gauss-Hermite BICM Mutual Information evaluation for PAM
 *
 *  This function evalautes the BICM Mutual Information
 *  (Generalized Mutual Information) for a real-valued
 *  constellation (e.g. PAM)
 * \param[in] C     Real constellation (e.g. PAM)
 * \param[in] M     Number of points in the constellation
 * \param[in] s     Standard deviation of real-valued AWGN
 * \param[in] Pk    Probability of each constellation point
 */
double pam_eval_gmi(const double *C, int M, double s, const double *Pk);

/*! \brief Gauss-Hermite AWGN Mutual Information evaluation for QAM
 *
 *  This function evalautes the AWGN Mutual Information for a complex-valued
 *  constellation (e.g. QAM)
 * \param[in] C     Complex constellation (e.g. QAM)
 * \param[in] M     Number of points in the constellation
 * \param[in] s     Standard deviation of complex-valued AWGN
 * \param[in] Pk    Probability of each constellation point
 */
double qam_eval_mi(const double complex *C, int M, double s, const double *Pk);

/*! \brief Gauss-Hermite BICM Mutual Information evaluation for QAM
 *
 *  This function evalautes the BICM Mutual Information
 *  (Generalized Mutual Information) for a complex-valued
 *  constellation (e.g. PAM)
 * \param[in] C     Real constellation (e.g. PAM)
 * \param[in] M     Number of points in the constellation
 * \param[in] s     Standard deviation of complex-valued AWGN
 * \param[in] Pk    Probability of each constellation point
 */
double qam_eval_gmi(const double complex *C, int M, double s, const double *Pk);

/*! \brief QAM BICM soft decision
 *
 *  This function calculates bit-wise log-likelihood ratios (LLRs) for a complex
 *  constellation (e.g. QAM), assuming a pragmatic (BICM) decoder over
 *  an AWGN channel.
 * \param[in]  y     Received complex samples
 * \param[in]  Ns    Number of complex samples
 * \param[in]  C     Complex constellation (e.g. QAM)
 * \param[in]  Pk    Probability of each constellation point
 * \param[in]  M     Number of constellation points
 * \param[in]  s2    Variance of complex-valued AWGN
 * \param[out] l     Bit-wise log-likelihood ratios
 */
void qam_soft_decode(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, const double *s2, double *l);

/*! \brief PAM BICM soft decision
 *
 *  This function calculates bit-wise log-likelihood ratios (LLRs) for a real
 *  constellation (e.g. PAM), assuming a pragmatic (BICM) decoder over
 *  an AWGN channel.
 * \param[in]  y     Received samples
 * \param[in]  Ns    Number of samples
 * \param[in]  C     Real constellation (e.g. PAM)
 * \param[in]  Pk    Probability of each constellation point
 * \param[in]  M     Number of constellation points
 * \param[in]  s2    Variance of real-valued AWGN
 * \param[out] l     Bit-wise log-likelihood ratios
 */
void pam_soft_decode(const double *y, int Ns, const double *C, 
        const double *Pk, int M, const double *s2, double *l);

/*! \brief QAM BICM soft decision with phase noise
 *
 *  This function calculates bit-wise log-likelihood ratios (LLRs) for a real
 *  constellation (e.g. PAM), assuming a pragmatic (BICM) decoder over
 *  an AWGN channel with Tikhonov-distributed memoryless phase noise.
 *  See: 10.1109/TWC.2014.040714.130731
 * \param[in]  y     Received samples
 * \param[in]  Ns    Number of samples
 * \param[in]  C     Real constellation (e.g. PAM)
 * \param[in]  Pk    Probability of each constellation point
 * \param[in]  M     Number of constellation points
 * \param[in]  Kn    Concentration (1/variance) of complex-valued AWGN
 * \param[in]  Kp    Concentration of Tikhonov-distributed phase noise
 * \param[out] l     Bit-wise log-likelihood ratios
 */
void qam_soft_decode_pn(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, double Kn, double Kp, double *l);

/*! \brief QAM BICM soft decision with phase noise - MAXLOG approximation
 *
 *  This function calculates bit-wise log-likelihood ratios (LLRs) for a real
 *  constellation (e.g. PAM), assuming a pragmatic (BICM) decoder over
 *  an AWGN channel with Tikhonov-distributed memoryless phase noise.
 *  LLRs are calculated with the MAX-LOG approximation, which is numerically
 *  more stable for large SNR.
 *  See: 10.1109/TWC.2014.040714.130731
 * \param[in]  y     Received samples
 * \param[in]  Ns    Number of samples
 * \param[in]  C     Real constellation (e.g. PAM)
 * \param[in]  Pk    Probability of each constellation point
 * \param[in]  M     Number of constellation points
 * \param[in]  Kn    Concentration (1/variance) of complex-valued AWGN
 * \param[in]  Kp    Concentration of Tikhonov-distributed phase noise
 * \param[out] l     Bit-wise log-likelihood ratios
 */
void qam_soft_decode_pn_maxlog(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, double Kn, double Kp, double *l);

/*! \brief QAM Monte-Carlo AWGN Mutual Information
 *
 *  This function calculates the AWGN Mutual Information for a 
 *  complex (e.g. QAM) constellation using Monte-Carlo integration. 
 * \param[in]  y     Received complex samples
 * \param[in]  Ns    Number of samples
 * \param[in]  C     Complex constellation (e.g. PAM)
 * \param[in]  Pk    Probability of each constellation point
 * \param[in]  M     Number of constellation points
 * \param[in]  s2    Variance of complex-valued AWGN
 */
double qam_montecarlo_mi(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, double s2);

/*! \brief PAM Monte-Carlo AWGN Mutual Information
 *
 *  This function calculates the AWGN Mutual Information for a 
 *  real (e.g. PAM) constellation using Monte-Carlo integration. 
 * \param[in]  y     Received samples
 * \param[in]  Ns    Number of samples
 * \param[in]  C     Real constellation (e.g. PAM)
 * \param[in]  Pk    Probability of each constellation point
 * \param[in]  M     Number of constellation points
 * \param[in]  s2    Variance of real-valued AWGN
 */
double pam_montecarlo_mi(const double *y, int Ns, const double *C,
        const double *Pk, int M, double s2);

/*! \brief QAM AWGN soft decision
 *
 *  This function calculates symbol-wise log-likelihood ratios (LLRs) for a complex
 *  constellation (e.g. QAM), assuming an AWGN channel.
 * \param[in]  y     Received complex samples
 * \param[in]  Ns    Number of complex samples
 * \param[in]  C     Complex constellation (e.g. QAM)
 * \param[in]  Pk    Probability of each constellation point
 * \param[in]  M     Number of constellation points
 * \param[in]  s2    Variance of complex-valued AWGN
 * \param[out] l     Symbol-wise log-likelihood ratios
 */
void qam_symbol_decode(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, double s2, double *l);

#endif /* CAPACITY_FUNCTIONS_H */