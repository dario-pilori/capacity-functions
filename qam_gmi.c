/*
* qam_gmi.c 	Simple program to calculate MI/GMI with Gauss-Hermite quadrature
*
* 2018 - Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
* MIT License
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "capacity_functions.h"
#define MAX_M_QAM 1024
#define CONST_DIR "/home/dario/Documents/MATLAB/dsp-library/text_constellations/"

int get_constellation(double complex *C, const char *cst_name);

int main(int argc, char *argv[])
{
	// Declare variables	
	double mi, gmi, snr, Es, sigma_n;
	char cst_name[10];
	int M_qam;
	double complex C[MAX_M_QAM];
	
	// Get SNR from command line
	if (argc >= 3) {
	  strncpy(cst_name, argv[1], 9);
	  snr = atof(argv[2]);
	} else {
	  printf("Usage: %s <constellation> <snr in dB>\n",argv[0]);
	  return 1;
	}

	// Get constellation
	M_qam = get_constellation(C, cst_name);
	if(M_qam == 0)
	  return 1;

	// Calculate symbol and noise energies
	Es = complex_symbol_energy(C, M_qam);
	sigma_n = sqrt(Es) * pow(10.0,-snr/20.0);
	
	// Calculate
	gmi = qam_eval_gmi(C, M_qam, sigma_n);
	mi = qam_eval_mi(C, M_qam, sigma_n);
	
	// Print results
	printf("MI=%f, GMI=%f\n",mi,gmi);	

	return 0;
}

// Helper function to get constellation from text file
int get_constellation(double complex *C, const char *cst_name)
{
  // Variable declaration
  FILE *fd;
  int M_qam = 0;
  double tmp1, tmp2;
  char *file_path;
  
  // Build full path of file
  file_path = malloc((strlen(cst_name)+strlen(CONST_DIR)+5)*sizeof(char));
  strcpy(file_path,CONST_DIR);
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


