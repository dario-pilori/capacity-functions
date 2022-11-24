/*
* example_pam_gmi.c
* Simple program to test pam_eval_mi and pam_eval_gmi functions.
*
* Copyright (c) 2018-2022 Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
* SPDX-License-Identifier: MIT
*/

#include <stdio.h>
#include <math.h>
#include "capacity_functions.h"
#define M_PAM 4

int main(int argc, char *argv[])
{
	// Indexes
	int i;
	
	// Initialize MI and GMI	
	double MI, GMI;
	
	// Constellation (natural bit mapping order)
	double C[M_PAM] = {-3.0, -1.0, 3.0, 1.0};
    // Probability for each symbol to be transmitted
	double Pk[M_PAM] = {0.5, 0.3, 0.15, 0.05};

	// Signal-to-noise ratio (dB)
	double snr = 10;
	double Es = symbol_energy(C, Pk, M_PAM);
	double  s = sqrt(Es)*pow(10.0,-snr/20.0);
	
	// Calculate
	GMI = pam_eval_gmi(C, M_PAM, s, Pk);
	MI= pam_eval_mi(C, M_PAM, s, Pk);
	
	// Print results
	printf("MI=%f, GMI=%f\n",MI,GMI);	

	return 0;
}

