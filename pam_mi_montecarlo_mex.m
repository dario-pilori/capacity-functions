%PAM_MI_MONTECARLO_MEX  Compute MI for PAM using MonteCarlo for AWGN
%   Usage: pam_mi_montecarlo_mex(C, sigma2, a-y, Pk)
%   C     :=   Real constellation in Gray-mapping order
%   sigma2:=   Noise variance (average)
%   a-y   :=   Received noise (transmitted minus received symbols)
%   Pk    :=   Probability of each constellation symbol
%  
%   Compile with: mex -lm -R2018a qam_llr_mex.c
%   (requires MATLAB R2018a or newer versions)
%   Designed for 64-bit Linux.
%  
%   Copyright (c) 2018 Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
%   SPDX-License-Identifier: MIT

%   MEX File function.