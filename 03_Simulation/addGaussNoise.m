function [V] = addGaussNoise(v,C,T)
% =========================================
% =========================================
%
% Add Noise
% By: Damennick Henry and Jordan Murphy
% Date: 12/19/18
% Description: Adds Gaussian(0,C) noise to vector quantities
% Inputs
%       v - Vectors to add noise to [n x m]
%       C - Covariance matrix
%       T - Noise to state matrix
% Outpus
%       V - Noisey vectors
%
% =========================================
% =========================================

% Noise vector size
[~, l] = size(T);
% Number of vectors
[~, m] = size(v);
% Cholesky factorization of covariance matrix
Sv = chol(C,'lower');
% Add noise
V = v + T*Sv*randn(l, m);
end

