%QAM_LLR_PN_MAXLOG_MEX  Compute LLRs for QAM
% Use this function to compute log-likelihood-ratios
% for M-QAM constellations assuming an AWGN channel
% with phase noise (10.1109/TWC.2014.040714.130731).
% The output of this function is comptabible to MATLAB's qamdemod()
% function of the Communication System Toolbox.
% This function uses the max-log approximation.
%
% Usage: l = qam_llr_pn_maxlog_mex(C, Kn, Kp, y, Pk)
% C     :=   Complex constellation in Gray-mapping order
% Kn    :=   1/Noise variance
% Kp    :=   ~ 1/Phase noise variance
% y     :=   Received complex symbols
% Pk    :=   Probability of each constellation symbol
%
% Compile with: mex -lm -R2018a qam_llr_maxlog_mex.c
% (requires MATLAB R2018a or newer versions)
% Designed for 64-bit Linux.
%
% Copyright (c) 2018 - Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
% MIT License

%   MEX File function.