function [F, Q] = ct2dtStochastic(A, T, W, dt)
% =========================================
% =========================================
%
% [F, Q] = ct2dtStochastic(A, T, W, dt)
%
% Continuous two discrete time dynamics
% By: Damennick Henry 
% Date: 10/2/18
% Description: Van Loan's method for dicretizing system dynamics and
% process noise
% Inputs
%       A  - Dynamics matrix
%       T  - Process noise to state matrix
%       W  - Process noise covariance matrix
%       dt - Sampling time
%
% =========================================
% =========================================

% Dimensions of A
[n,~] = size(A);

% Van Loan matrix
Z = dt*[-A T*W*T'; zeros(n) A'];
% Van Loan matrix exponential
eZ = expm(Z);
% State Transition matrix
F = eZ(end-(n-1):end, end-(n-1):end)';
% Discretized process noise covariance matrix
Q = F*eZ(1:n, end-(n-1):end);

end