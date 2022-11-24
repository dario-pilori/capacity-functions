function l = qam_llr(C, sigma2, y, Pk)
%QAM_LLR  Compute LLRs for QAM
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
% Copyright (c) 2018-2022 Dario Pilori
% SPDX-License-Identifier: MIT

% TODO: validate input

% Constellation size
M = size(C,1);

% Cycle through constellation bits
l = NaN(size(y,1),log2(M));

for k = 1:log2(M)
    % Get sets where the k-th bit of the constellation is 0 and 1
    c_1 = reshape((1:2^k:M) + (0:(2^(k-1)-1))',M/2,1);
    c_0 = c_1 + 2^(k-1);
    
    % Evaluate LLR assuming an AWGN channel
    l(:,k) = log(sum(exp(-abs(y-C(c_1).').^2/sigma2).*Pk(c_1).',2) ./...
        sum(exp(-abs(y-C(c_0).').^2/sigma2).*Pk(c_0).',2));
end

% Compatibility: convert to the same output format of qam_llr_mex
l = reshape(l.',[],1);