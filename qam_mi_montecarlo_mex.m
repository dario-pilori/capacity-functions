%QAM_MI_MONTECARLO_MEX  Compute MI for QAM using MonteCarlo for AWGN
% Compile with: mex -lm -R2018a qam_llr_mex.c
% (requires MATLAB R2018a or newer versions)
% Designed for 64-bit Linux.
%
% Usage: mi = qam_mi_montecarlo_mex(C, sigma2, a-y, Pk)
% C     :=   Complex constellation in Gray-mapping order
% sigma2:=   Noise variance (average)
% a-y   :=   Received noise (transmitted-received symbols)
% Pk    :=   Probability of each constellation symbol
%
% Copyright (c) 2018 - Dario Pilori, Politecnico di Torino <dario.pilori@polito.it>
% MIT License

%   MEX File function.