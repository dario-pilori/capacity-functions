function [GMI, MI] = qam_gmi(C, SNR, Pk)
%QAM_GMI Compute MI and GMI for QAM
%  Use this function to compute MI and GMI for M-QAM
%  constellations in an AWGN channel using Gauss-Hermite quadrature.
%
%  Usage: [gmi,mi] = qam_gmi_mex(C, SNR, Pk)
%  C    :=   Complex constellation in Gray-mapping order
%  SNR  :=   Vector of SNR (Es/No, dB)
%  Pk   :=   Probability of each constellation symbol
%
%  Copyright (c) 2018-2022 Dario Pilori <dario.pilori@polito.it>
%  SPDX-License-Identifier: MIT

%% Validate input
validateattributes(C,{'numeric'},{'column'},'','C',1);
validateattributes(SNR,{'double'},{'column'},'','SNR',2);
validateattributes(Pk,{'double'},{'column','nrows',size(C,1)},'','Pk',3);

%% Preliminary computations
% Gauss-Hermite parameters
% Source: https://www.efunda.com/math/num_integration/findgausshermite.cfm
N_GH = 10;
x = [-3.436159118837737603327
    -2.532731674232789796409
    -1.756683649299881773451
    -1.036610829789513654178
    -0.3429013272237046087892
    0.3429013272237046087892
    1.036610829789513654178
    1.756683649299881773451
    2.532731674232789796409
    3.436159118837737603327];
w = [7.64043285523262062916E-6
    0.001343645746781232692202
    0.0338743944554810631362
    0.2401386110823146864165
    0.6108626337353257987836
    0.6108626337353257987836
    0.2401386110823146864165
    0.03387439445548106313617
    0.001343645746781232692202
    7.64043285523262062916E-6];

% Constellation size
M = size(C,1);
m = log2(M);

% Evaluate bit-wise probabilities
Pb = NaN(m,1);
for k = 1:m
    % Get sets where the k-th bit of the constellation is 1
    c_1 = reshape((1:2^k:M) + (0:(2^(k-1)-1))',M/2,1);

    % Compute probabilities
    Pb(k) = sum(Pk(c_1));
end

% Compute variance of AWGN
sigma2 = mean(abs(C).^2)*10.^(-SNR/10);

%% Calculate GMI
% For each SNR point
GMI = zeros(size(sigma2,1),1);
for s = 1:size(sigma2,1)
    % For each constellation bit
    for k = 1:m
        % Bit can be either zero or one
        for b = [false,true]
            % Get sets where the k-th bit of the constellation is b
            c_b = reshape((1:2^k:M) + (0:(2^(k-1)-1))',M/2,1) + 2^(k-1)*b;

            for i = 1:M/2
                % Compute numerator and denominator of the logarithm
                num = squeeze( ...
                    sum(...
                    exp(-(abs(C(c_b(i))-C).^2 - ...
                    2*sqrt(sigma2(s))*real((reshape(x,1,1,N_GH)+1j*x.').* ...
                    (C(c_b(i))-C)))/sigma2(s)).* ...
                    Pk));

                den = squeeze( ...
                    sum(...
                    exp(-(abs(C(c_b(i))-C(c_b)).^2 - ...
                    2*sqrt(sigma2(s))*real((reshape(x,1,1,N_GH)+1j*x.').* ...
                    (C(c_b(i))-C(c_b))))/sigma2(s)).* ...
                    Pk(c_b)));

                % Apply bit probability
                if b
                    den = den/(1-Pb(k));
                else
                    den = den/Pb(k);
                end

                % Compute GMI
                GMI(s) = GMI(s) - sum(sum(w .* w.' .* log2(num./den)*Pk(c_b(i))));
            end
        end
    end
end

% Add entropy of constellation and remove entropy of each bit
GMI = max(0, GMI/pi - sum(Pk.*log2(Pk)) + sum(Pb.*log2(Pb)+(1-Pb).*log2(1-Pb)));

%% Calculate MI, if needed
if nargout>1
    MI = zeros(size(sigma2,1),1);
    % For each SNR
    for s = 1:size(sigma2,1)
        % For each constellation point
        for i = 1:M
            tmp = squeeze(...
                sum(...
                exp(-(abs(C(i)-C).^2 - ...
                2*sqrt(sigma2(s))*real((reshape(x,1,1,N_GH)+1j*x.').*(C(i)-C)))/...
                sigma2(s)).*...
                Pk...
                ));
            MI(s) = MI(s) - sum(sum(w .* w.' .* log2(tmp)))*Pk(i);
        end
    end
    MI = max(0, MI/pi);
end
