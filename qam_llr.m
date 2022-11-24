function l = qam_llr(C, sigma2, y, Pk)
%QAM_LLR  Compute LLRs
% Use this function to compute log-likelihood-ratios
% for M-QAM and M-PAM constellations assuming an AWGN channel.
%
% Usage: l = qam_llr_mex(C, sigma2, y, Pk)
% C     :=   Complex constellation in Gray-mapping order
% sigma2:=   Noise variance per each constellation point
% y     :=   Received complex symbols
% Pk    :=   Probability of each constellation symbol
%
% Copyright (c) 2018-2022 Dario Pilori <d.pilori@inrim.it>
% SPDX-License-Identifier: MIT

% Validate input
validateattributes(C,{'numeric'},{'column'},'','C',1);
if (numel(sigma2)==1)
    sigma2 = sigma2*ones(size(C));
end
validateattributes(sigma2,{'double'},{'column','nrows',size(C,1)},'','sigma2',2);
validateattributes(y,{'double'},{'column'},'','y',3);
validateattributes(Pk,{'double'},{'column','nrows',size(C,1)},'','Pk',4);

% Constellation size
M = size(C,1);
m = log2(M);

% Cycle through constellation bits
l = NaN(m, size(y,1));

for k = 1:m
    % Get sets where the k-th bit of the constellation is 0 and 1
    c_1 = reshape((1:2^k:M) + (0:(2^(k-1)-1))',M/2,1);
    c_0 = c_1 + 2^(k-1);
    
    % Evaluate LLR assuming an AWGN channel
    l(k,:) = log(sum(exp(-abs(y.'-C(c_1)).^2./sigma2(c_0)).*Pk(c_1)) ./...
        sum(exp(-abs(y.'-C(c_0)).^2./sigma2(c_0)).*Pk(c_0)));
end

% Compatibility: convert to the same output format of qam_llr_mex
l = l(:);