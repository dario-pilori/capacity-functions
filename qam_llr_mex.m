%QAM_LLR_MEX  Compute LLRs for QAM
% Use this function to compute log-likelihood-ratios
% for M-QAM constellations assuming an AWGN channel.
% The output of this function is comptabible to MATLAB's qamdemod()
% function of the Communication System Toolbox.
%
% Usage: l = qam_llr_mex(C, sigma2, y, Pk)
% C     :=   Complex constellation in Gray-mapping order
% sigma2:=   Noise variance per each constellation point
% y     :=   Received complex symbols
% Pk    :=   Probability of each constellation symbol
%
% Compile with: mex -lm -R2018a qam_llr_mex.c
% (requires MATLAB R2018a or newer versions)
% Designed for 64-bit Linux.
%
% Copyright (c) 2018 - Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
% SPDX-License-Identifier: MIT

%   MEX File function.