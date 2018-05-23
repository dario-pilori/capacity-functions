/*
* qam_gmi_sweep.c 	Calculate MI/GMI with Gauss-Hermite quadrature with multiple SNR
*
* 2018 - Dario Pilori <dario.pilori@polito.it>
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <omp.h>
#include "capacity_functions.h"
#define MAX_M_QAM 1024
#define CONST_DIR "/Documents/MATLAB/dsp-library/text_constellations/"

int get_constellation(double complex *C, const char *cst_name);
double *linspace(double start, double end, unsigned int num);

int main(int argc, char *argv[])
{
	// Declare variables	
	double Es, sigma_n, start_snr, end_snr, tmp_Pk, lambda;
    double *snr, *mi, *gmi, *Pk;
	double complex C[MAX_M_QAM];
	int M_qam, num_snr, i;
	char cst_name[10];
    const char *out_file = "mi_gmi_res.txt";
    FILE *fd;
	
	// Get values from command line
	if (argc >= 6) {
	  strncpy(cst_name, argv[1], 9);
	  start_snr = atof(argv[2]);
	  end_snr = atof(argv[3]);
	  num_snr = atoi(argv[4]);
      lambda = atof(argv[5]);
	} else {
	  printf("Usage: %s <constellation> <start_snr> <end_snr> <num_snr> <lambda>\n",argv[0]);
	  return 1;
	}

	// Get constellation
	M_qam = get_constellation(C, cst_name);
	if(M_qam == 0)
	  return 1;
    
    // Calculate probabilities
    tmp_Pk = 0;
    Pk = (double *) malloc(M_qam*sizeof(double));
    for(i=0; i<M_qam; i++) {
        Pk[i] = exp(-lambda*pow(cabs(C[i]),2));
        tmp_Pk += Pk[i];
    }
    
    // Normalize probability
    for(i=0; i<M_qam; i++) {
        Pk[i] /= tmp_Pk;
    }    
    
    // Calculate symbol energy
    Es = complex_symbol_energy(C, Pk, M_qam);

    // Get snr
    snr = linspace(start_snr, end_snr, num_snr);
    if(snr == NULL)
        return 1;

    // Allocate MI and GMI
    mi = (double *) malloc(num_snr*sizeof(double));
    gmi = (double *) malloc(num_snr*sizeof(double));
    
    // Parallel evaluation of MI and GMI
    #pragma omp parallel for private(sigma_n)
    for(i=0; i<num_snr; i++) {
	    sigma_n = sqrt(Es) * pow(10.0,-snr[i]/20.0);
	    gmi[i] = qam_eval_gmi(C, M_qam, sigma_n, Pk);
	    mi[i] = qam_eval_mi(C, M_qam, sigma_n, Pk);
    }
    
    // Free probabilities
    free(Pk);

    // Print results
    printf("  SNR    |    MI    |   GMI\n------------------------------\n");
    for(i=0; i<num_snr; i++) {
	    printf("%f | %f | %f\n",snr[i],mi[i],gmi[i]);	
    }

    // Save results to file
    fd = fopen(out_file,"w");
    for(i=0; i<num_snr; i++)
	    fprintf(fd,"%f %f %f\n",snr[i],mi[i],gmi[i]);	
    fclose(fd);

    // Free memory and return
    free(snr);
    free(mi);
    free(gmi);

	return 0;
}

// Helper function to allocate the SNR
double *linspace(double start, double end, unsigned int num)
{
    // Variables
    double inc;
    double *vec;
    int i;

    // Check parameters
    if(end<=start || num>100000)
        return NULL;

    // Allocate vector of SNR
    vec = (double *) malloc(num*sizeof(double));

    // Fill
    inc = (end - start)/(num - 1);
    for(i=0; i<num-1; i++)
        vec[i] = start + inc*i;
    vec[num-1] = end;

    // Return
    return vec;
}

// Helper function to get constellation from text file
int get_constellation(double complex *C, const char *cst_name)
{
  // Variable declaration
  FILE *fd;
  int M_qam = 0;
  double tmp1, tmp2;
  char *file_path;
  const char *homedir = getenv("HOME");
  
  // Build full path of file
  file_path = malloc((strlen(homedir)+strlen(cst_name)+strlen(CONST_DIR)+5)*sizeof(char));
  strcpy(file_path,homedir);
  strcat(file_path,CONST_DIR);
  strcat(file_path,cst_name);
  strcat(file_path,".txt");
  
  // Open file and erase path
  fd = fopen(file_path,"r");
  free(file_path);
  if (fd == NULL) {
    printf("Unable to open %s\n",cst_name);
    return 0;
  }
  
  // Read all constellation points
  while (fscanf(fd,"%lf %lf",&tmp1,&tmp2)==2 && M_qam<MAX_M_QAM) {
    C[M_qam] = tmp1 + I*tmp2;
    M_qam++;
  }
  
  // Close and exit
  fclose(fd);
  return M_qam;
}


