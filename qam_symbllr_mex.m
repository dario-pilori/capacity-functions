%QAM_SYMBLLR_MEX  Compute LLRs for QAM
% Use this function to compute symbol-wise log-likelihood-ratios
% for M-QAM constellations assuming an AWGN channel.
%
% Usage: l = qam_symbllr_mex(C, sigma2, y, Pk)
% C     :=   Complex constellation in Gray-mapping order
% sigma2:=   Noise variance
% y     :=   Received complex symbols
% Pk    :=   Probability of each constellation symbol
%
% Compile with: mex -lm -R2018a qam_llr_mex.c
% (requires MATLAB R2018a or newer versions)
% Designed for 64-bit Linux.
%
% Copyright (c) 2018 - Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
% MIT License

%   MEX File function.