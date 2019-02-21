%PAM_LLR_MEX  Compute LLRs for PAM
%   Use this function to compute log-likelihood-ratios
%   for M-PAM constellations assuming an AWGN channel.
%
%   Usage: pam_llr_mex(C, sigma2, y, Pk)
%   C     :=   Real constellation in Gray-mapping order
%   sigma2:=   Noise variance per each constellation point
%   y     :=   Received complex symbols
%   Pk    :=   Probability of each constellation symbol
%
%   Compile with: mex -lm -R2018a pam_llr_mex.c
%   (requires MATLAB R2018a or newer versions)
%   Designed for 64-bit Linux.
%
%   Copyright (c) 2018 Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
%   MIT License

%   MEX File function.