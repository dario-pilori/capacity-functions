function mi = montecarlo_mi(C, sigma2, y, Pk)
%MONTECARLO_MI  Compute MI using Monte-Carlo for AWGN channels
% Use this function to compute the AWGN Mutual Information from a set
% of received noisy symbols.
%
% Usage: mi = montecarlo_mi(C, sigma2, a-y, Pk)
% C     :=   Complex constellation in Gray-mapping order
% sigma2:=   Noise variance (average)
% a-y   :=   Received noise (transmitted-received symbols)
% Pk    :=   Probability of each constellation symbol
%
% Copyright (c) 2018-2022 Dario Pilori <d.pilori@inrim.it>
% SPDX-License-Identifier: MIT

% Validate input
validateattributes(C,{'numeric'},{'column'},'','C',1);
validateattributes(sigma2,{'double'},{'scalar'},'','sigma2',2);
validateattributes(y,{'double'},{'column'},'','y',3);
validateattributes(Pk,{'double'},{'column','nrows',size(C,1)},'','Pk',4);

% Get constellation size
M = size(C,1);
mi = 0;

% Calculate MI for every constellation point
for i = 1:M
    mi = mi - ...
        Pk(i)*sum(log2(sum(exp(-((C.'-C(i)).^2 + 2*real(y .* C.'))/sigma2).*Pk.',2)));
end

% Sanity check and divide by number of symbols
mi = max(0, mi/size(y,1));