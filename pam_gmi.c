#include <stdio.h>
#include <math.h>
#include "capacity_functions.h"
#define M_PAM 4

int main(int argc, char *argv[])
{
	// Initialize MI and GMI	
	double MI, GMI;
	
	// Constellation (natural bit mapping order)
	double C[M_PAM] = {-3.0, -1.0, 3.0, 1.0};
	double Pk[M_PAM] = {0.25, 0.25, 0.25, 0.25};

	// Signal-to-noise ratio (dB)
	double snr = 12;
	double Es = symbol_energy(C, Pk, M_PAM);
	double  s = sqrt(Es)*pow(10.0,-snr/20.0);
	
	// Calculate
	GMI = pam_eval_gmi(C, M_PAM, s);
	MI= pam_eval_mi(C, M_PAM, s, Pk);
	
	// Print results
	printf("MI=%f, GMI=%f\n",MI,GMI);	

	return 0;
}
