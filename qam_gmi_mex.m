%QAM_GMI_MEX  Compute MI and GMI for QAM
%  Use this function to compute MI and GMI for M-QAM
%  constellations in an AWGN channel using Gauss-Hermite quadrature.
%
%  Usage: [gmi,mi] = qam_gmi_mex(C, SNR, Pk)
%  C    :=   Complex constellation in Gray-mapping order
%  SNR  :=   Vector of SNR (Es/No, dB)
%  Pk   :=   Probability of each constellation symbol
%
%  Compile with: mex -lm -R2018a CFLAGS="$CFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" qam_gmi_mex.c
%  (requires MATLAB R2018a or newer versions)
%  Designed for 64-bit Linux.
%
%  Copyright (c) 2018 Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
%  SPDX-License-Identifier: MIT

%  MEX File function.