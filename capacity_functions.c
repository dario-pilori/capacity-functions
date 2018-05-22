/*
* capacity_functions.c
* Set of functions to evaluate capacity
*
* 2018 - Dario Pilori <dario.pilori@polito.it>
*/
// Includes
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include "capacity_functions.h"

// Gauss-Hermite parameters
const double x[N_GH] = {-3.436159118837737603327,
-2.532731674232789796409,
-1.756683649299881773451,
-1.036610829789513654178,
-0.3429013272237046087892,
0.3429013272237046087892,
1.036610829789513654178,
1.756683649299881773451,
2.532731674232789796409,	
3.436159118837737603327};	
const double w[N_GH] = {7.64043285523262062916E-6,
0.001343645746781232692202,
0.0338743944554810631362,
0.2401386110823146864165,
0.6108626337353257987836,
0.6108626337353257987836,
0.2401386110823146864165,
0.03387439445548106313617,
0.001343645746781232692202,
7.64043285523262062916E-6};

// Helper function to calculate symbol energy
double symbol_energy(const double *C, const double *Pk, int M)
{
	int i;
	double Es = 0.0;
	
	for(i=0;i<M;i++)
	{
		Es = Es + Pk[i]*pow(C[i],2.0);
	}	
	
	return Es;
}

// Helper function to insert a zero at position k of number i with nb bits
unsigned int insert_zero(unsigned int i, unsigned int k, unsigned int nb)
{
  unsigned int b0, left, right;
  
  left = (i<<1) & (  ( (1<<(nb-k)) - 1)<<(k+1) );
  right = i & ((1<<k)-1);
  b0 = left | right;
  
  return b0;
}

// Helper function to calculate symbol energy
double complex_symbol_energy(const double complex *C, const double *Pk, int M)
{
	int i;
	double Es = 0.0;
	
	for(i=0;i<M;i++)
	{
		Es += Pk[i]*pow(cabs(C[i]),2.0);
	}	
	
	return Es;
}

// Helper function to evaluate log(sum(exp(.))) in a safe way
void maxxx(double *a, double b)
{
    if(*a>b)
        *a += log(1.0 + exp(b-*a));
    else
        *a = b + log(1.0 + exp(*a-b));       
}

/* 
 * Calculate AWGN mutual information (MI) for PAM using Gauss-Hermite 
 * qadrature assuming an AWGN channel.
 */
double pam_eval_mi(const double *C, int M, double s, const double *Pk)
{
  double MI = 0.0;
  double tmp;
  int i, l, j;

  // Cycle through constellation point
  for(i=0; i<M; i++)
  {
     // Cycle through Gauss-Hermite parameters
     for(l=0; l<N_GH; l++)
     {
       tmp = 0.0;
       
       // Cycle through constellation point for the logarithm
       for(j=0; j<M; j++)
       {
         tmp += Pk[j]*exp(-(pow(C[j]-C[i],2.0) - sqrt(8.0)*x[l]*s*(C[j]-C[i]))/(2*pow(s,2.0)));
       }
       
       MI -= Pk[i]*w[l]*log2(tmp);
     }
  }
  
  MI /= sqrt(M_PI);  
  return MI;
}

/*
 * Calculate AWGN mutual information (MI) for PAM using MonteCarlo
 * integration assuming an AWGN channel.
 */
double pam_montecarlo_mi(const double *y, int Ns, const double *C,
        const double *Pk, int M, double s2)
{
    double MI = 0;
    double tmp;
    int i, l, j;
    const unsigned int m = log2(M);
    
    // For each received symbol
    #pragma omp parallel for private(tmp,i,j)
    for(l=0; l<Ns; l++)
    {
        // Cycle through constellation point
        for(i=0; i<M; i++)
        {
            tmp = 0.0;
            
            // Cycle through constellation point for the logarithm
            for(j=0; j<M; j++)
            {
                tmp += exp(-(pow(C[j]-C[i],2.0) +
                        2*y[l]*(C[j]-C[i]))/(2*s2))*Pk[j];
            }
            
            MI -= log2(tmp)*Pk[i];
        }
    }
    
    // Prepare output
    MI /= Ns;
    
    return MI;
}


/* 
 * Calculate BICM mutual information (GMI) for PAM using Gauss-Hermite 
 * qadrature assuming an AWGN channel.
 */
double pam_eval_gmi(const double *C, int M, double s, const double *Pk)
{
  // Variables
  const unsigned int m = log2(M);
  unsigned int i, l, j, k, b;
  unsigned int bi, bj;
  double GMI = 0.0;
  double tmp_num, tmp_den;
  double *Pb0;
  
  // Calculate bit-wise probabilities
  Pb0 = malloc(sizeof(double)*m);
  for(k=0; k<m; k++)
  {
    Pb0[k] = 0.0;
    for(j=0; j<M/2; j++)
    {
      bj = insert_zero(j, k, m);
      Pb0[k] += Pk[bj];
    }
  }   
  
  // Cycle through constellation bit
  for(k=0; k<m; k++)
  {  
  // Bit can be either 0 or 1
    for(b=0; b<=1; b++)
    {
      // Constellation points where k-th bit is equal to b
      for(i=0; i<M/2; i++)
      {
        bi = insert_zero(i, k, m) + (b<<k);
        
        // Cycle through Gauss-Hermite parameters
        for(l=0; l<N_GH; l++)
        {
          // Initialize numerator and denominator of the logarithm
          tmp_num = 0.0;
          tmp_den = 0.0;
       
          // Numerator of the logarithm
          for(j=0; j<M; j++)
          {
            tmp_num += exp(-(pow(C[bi]-C[j],2.0) - sqrt(8.0)*x[l]*s*(C[bi]-C[j]))/(2*pow(s,2.0)))*Pk[j];
          }
        
          // Denominator of the logarithm
          for(j=0; j<M/2; j++)
          {
            bj = insert_zero(j, k, m) + (b<<k);
            tmp_den += exp(-(pow(C[bi]-C[bj],2.0) - sqrt(8.0)*x[l]*s*(C[bi]-C[bj]))/(2*pow(s,2.0)))*Pk[bj];
          }
          
          // Apply bit probability
          if(b)
            tmp_den /= 1-Pb0[k];
          else
            tmp_den /= Pb0[k];
       
          // Evaluate GMI
          GMI -= w[l]*log2(tmp_num/tmp_den)*Pk[bi];
        } 
      }
    }
  }
  
  // Divide by sqrt(pi) according to G-H
  GMI /= sqrt(M_PI);
  
  // Add entropy of constellation
  for(j=0; j<M; j++)
  {
    GMI -= Pk[j]*log2(Pk[j]);
  }
  
  // Remove entropy of each bit
  for(k=0; k<m; k++)
  {
    GMI += Pb0[k]*log2(Pb0[k]);
    GMI += (1-Pb0[k])*log2(1-Pb0[k]);
  }
  // Free memory and return
  free(Pb0);

  // Sanity check
  if(GMI < 0)
    GMI = 0.0;
  
  return GMI;
}

/* 
 * Calculate AWGN mutual information (MI) for QAM using Gauss-Hermite 
 * qadrature assuming an AWGN channel.
 */
double qam_eval_mi(const double complex *C, int M, double s, const double *Pk)
{
  double MI = 0;
  double tmp;
  int i, l1, l2, j;
  const unsigned int m = log2(M);

  // Cycle through constellation point
  for(i=0; i<M; i++)
  {
     // Cycle through Gauss-Hermite parameters
     for(l1=0; l1<N_GH; l1++)
     {
       for(l2=0; l2<N_GH; l2++)
       {
         tmp = 0.0;
       
         // Cycle through constellation point for the logarithm
         for(j=0; j<M; j++)
         {
           tmp += exp(-( pow(cabs(C[j]-C[i]),2.0) - 
                  2*s*creal((x[l1]+I*x[l2])*(C[j]-C[i])) )/pow(s,2.0))*Pk[j];
         }
       
         MI -= w[l1]*w[l2]*log2(tmp)*Pk[i];
       }
     }
  }
  
  // Prepare output
  MI /= M_PI;
  
  return MI;
}

/*
 * Calculate AWGN mutual information (MI) for QAM using MonteCarlo
 * integration assuming an AWGN channel.
 */
double qam_montecarlo_mi(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, double s2)
{
    double MI = 0;
    double tmp;
    int i, l, j;
    const unsigned int m = log2(M);
    
    // For each received symbol
    #pragma omp parallel for private(tmp,i,j) reduction(-:MI)
    for(l=0; l<Ns; l++)
    {
        // Cycle through constellation point
        for(i=0; i<M; i++)
        {
            tmp = 0.0;
            
            // Cycle through constellation point for the logarithm
            for(j=0; j<M; j++)
            {
                tmp += exp(-(pow(cabs(C[j]-C[i]),2.0) +
                        2*creal(y[l]*(C[j]-C[i])))/s2)*Pk[j];
            }
            
            MI -= log2(tmp)*Pk[i];
        }
    }
    
    // Prepare output
    MI /= Ns;
    
    return MI;
}

/* 
 * Calculate BICM mutual information (GMI) for QAM using Gauss-Hermite 
 * qadrature assuming an AWGN channel.
 */
double qam_eval_gmi(const double complex *C, int M, double s, const double *Pk)
{
  const int m = log2(M);
  
  int i, l1, l2, j, k, b;
  int bi, bj;
  double GMI = 0.0;
  double tmp_num, tmp_den;
  double *Pb0;
  
  // Calculate bit-wise probabilities
  Pb0 = malloc(sizeof(double)*m);
  for(k=0; k<m; k++)
  {
    Pb0[k] = 0.0;
    for(j=0; j<M/2; j++)
    {
      bj = insert_zero(j, k, m);
      Pb0[k] += Pk[bj];
    }
  }   

  // Cycle through constellation bit
  for(k=0; k<m; k++)
  {  
  // Bit can be either 0 or 1
    for(b=0; b<=1; b++)
    {
      // Constellation points where k-th bit is equal to b
      for(i=0; i<M/2; i++)
      {
        bi = insert_zero(i, k, m) + (b<<k);
        
        // Cycle through Gauss-Hermite parameters
        for(l1=0; l1<N_GH; l1++)
        {
          for(l2=0; l2<N_GH; l2++)
          {
            // Initialize numerator and denominator of the logarithm
            tmp_num = 0.0;
            tmp_den = 0.0;
       
            // Numerator of the logarithm
            for(j=0; j<M; j++)
            {
              tmp_num += exp(-(pow(cabs(C[bi]-C[j]),2.0) - 
              	         2*s*creal((x[l1]+I*x[l2])*(C[bi]-C[j])))/pow(s,2.0))*Pk[j];
            }
        
            // Denominator of the logarithm
            for(j=0; j<M/2; j++)
            {
              bj = insert_zero(j, k, m) + (b<<k);
              tmp_den += exp(-(pow(cabs(C[bi]-C[bj]),2.0) - 
              	         2*s*creal((x[l1]+I*x[l2])*(C[bi]-C[bj])))/pow(s,2.0))*Pk[bj];
            }
       
          // Apply bit probability
          if(b)
            tmp_den /= 1-Pb0[k];
          else
            tmp_den /= Pb0[k];
            
          // Evaluate GMI
          GMI -= w[l1]*w[l2]*log2(tmp_num/tmp_den)*Pk[bi];
          
          }
        } 
      }
    }
  }
  
  // Prepare output
  GMI /= M_PI;
  
  // Add entropy of constellation
  for(j=0; j<M; j++)
  {
    GMI -= Pk[j]*log2(Pk[j]);
  }
  
  // Remove entropy of each bit
  for(k=0; k<m; k++)
  {
    GMI += Pb0[k]*log2(Pb0[k]);
    GMI += (1-Pb0[k])*log2(1-Pb0[k]);
  }  
  
  // Free memory
  free(Pb0);
  
  // Sanity check if GMI is negative
  if(GMI < 0)
    GMI = 0.0;
  
  return GMI;
}

/* 
 * Calculate Log-Likelihood Ratios (LLRs) for QAM assuming a BICM-AWGN channel.
 */
void qam_soft_decode(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, const double *s2, double *l)
{
  const int m = log2(M);
  int i, k, j, bj;
  double tmp_num, tmp_den;
  
  // Cycle through received symbol
  #pragma omp parallel for private(tmp_num,tmp_den,bj,k,j)
  for(i=0; i<Ns; i++)
  {
    // Cycle through constellation bit
    for(k=m-1; k>=0; k--)
    {
      // Initialize numerator and denominator of the logarithm
      tmp_num = 0.0;
      tmp_den = 0.0;
      
      // Numerator of the logarithm
      for(j=0; j<M/2; j++)
      {
        bj = insert_zero(j, k, m);
        tmp_num += exp(-pow(cabs(y[i]-C[bj]),2.0)/s2[bj])*Pk[bj];
      }
      
      // Denominator of the logarithm
      for(j=0; j<M/2; j++)
      {
        bj = insert_zero(j, k, m) + (1<<k);
        tmp_den += exp(-pow(cabs(y[i]-C[bj]),2.0)/s2[bj])*Pk[bj];
      }
      
      // Calculate LLR
      l[i*m + k] = log(tmp_num/tmp_den);
    }      
  }
}

/* 
 * Calculate Log-Likelihood Ratios (LLRs) for PAM assuming a BICM-AWGN channel.
 */
void pam_soft_decode(const double *y, int Ns, const double *C, 
        const double *Pk, int M, const double *s2, double *l)
{
  const int m = log2(M);
  int i, k, j, bj;
  double tmp_num, tmp_den;
  
  // Cycle through received symbol
  #pragma omp parallel for private(tmp_num,tmp_den,bj,k,j)
  for(i=0; i<Ns; i++)
  {
    // Cycle through constellation bit
    for(k=m-1; k>=0; k--)
    {
      // Initialize numerator and denominator of the logarithm
      tmp_num = 0.0;
      tmp_den = 0.0;
      
      // Numerator of the logarithm
      for(j=0; j<M/2; j++)
      {
        bj = insert_zero(j, k, m);
        tmp_num += exp(-pow(y[i]-C[bj],2.0)/(2*s2[bj]))*Pk[bj];
      }    
      
      // Denominator of the logarithm
      for(j=0; j<M/2; j++)
      {
        bj = insert_zero(j, k, m) + (1<<k);
        tmp_den += exp(-pow(y[i]-C[bj],2.0)/(2*s2[bj]))*Pk[bj];
      }
      
      // Calculate LLR
      l[i*m + k] = log(tmp_num/tmp_den);
    }      
  }
}


/* 
 * Calculate Log-Likelihood Ratios (LLRs) for QAM assuming a BICM-AWGN channel
 * with phase noise.
 */
void qam_soft_decode_pn(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, double Kn, double Kp, double *l)
{
    const int m = log2(M);
    int i, k, j, bj;
    double tmp_num, tmp_den;
            
    // Cycle through received symbol
    #pragma omp parallel for private(tmp_num,tmp_den,bj,k,j)
    for(i=0; i<Ns; i++)
    {        
        // Cycle through constellation bit
        for(k=m-1; k>=0; k--)
        {
            // Initialize numerator and denominator of the logarithm
            tmp_num = 0.0;
            tmp_den = 0.0;
            
            // Numerator of the logarithm
            for(j=0; j<M/2; j++)
            {
                bj = insert_zero(j, k, m);                
                
                tmp_num += exp(-Kn/2*pow(cabs(C[bj]),2)+
                        cabs(conj(y[i])*C[bj]*Kn+Kp))*Pk[bj];
            }
            
            // Denominator of the logarithm
            for(j=0; j<M/2; j++)
            {
                bj = insert_zero(j, k, m) + (1<<k);
                
                tmp_den += exp(-Kn/2*pow(cabs(C[bj]),2)+
                        cabs(conj(y[i])*C[bj]*Kn+Kp))*Pk[bj];
            }
            
            // Calculate LLR
            l[i*m + k] = log(tmp_num)-log(tmp_den);
        }
    }
}

/* 
 * Calculate Log-Likelihood Ratios (LLRs) for QAM assuming a BICM-AWGN channel
 * with phase noise. This function uses the max-log approximation
 */
void qam_soft_decode_pn_maxlog(const double complex *y, int Ns, const double complex *C,
        const double *Pk, int M, double Kn, double Kp, double *l)
{
    const int m = log2(M);
    int i, k, j, bj;
    double tmp_num, tmp_den, max_num, max_den;
    double alpha, beta;
        
    // Cycle through received symbol
    #pragma omp parallel for private(tmp_num,tmp_den,bj,k,j,alpha,beta,max_num,max_den)
    for(i=0; i<Ns; i++)
    {
        // Cycle through constellation bit
        for(k=m-1; k>=0; k--)
        {
            // Initialize numerator and denominator of the logarithm
            max_num = 0.0;
            max_den = 0.0;
            
            // Numerator of the logarithm
            for(j=0; j<M/2; j++)
            {
                bj = insert_zero(j, k, m);                
                alpha = Kp + Kn*creal(conj(y[i])*C[bj]);
                beta = -Kn*cimag(conj(y[i])*C[bj]);
                
                tmp_num = -pow(cabs(C[bj]),2)*Kn/2 +
                        sqrt(pow(alpha,2)+pow(beta,2)) + 
                        log(Pk[bj]);
                
                if(tmp_num>max_num)
                    max_num = tmp_num;
            }
            
            // Denominator of the logarithm
            for(j=0; j<M/2; j++)
            {
                bj = insert_zero(j, k, m) + (1<<k);
                alpha = Kp + Kn*creal(conj(y[i])*C[bj]);
                beta = -Kn*cimag(conj(y[i])*C[bj]);
                
                tmp_den = -pow(cabs(C[bj]),2)*Kn/2 +
                        sqrt(pow(alpha,2)+pow(beta,2)) +
                        log(Pk[bj]);
                
                if(tmp_den>max_den)
                    max_den = tmp_den;
            }
            
            // Calculate LLR
            l[i*m + k] = max_num-max_den;
        }
    }
}

