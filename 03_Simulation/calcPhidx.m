function [phidx] = calcPhidx(xnom, XS, YS, range)
% =========================================
% =========================================
%
% [phidx] = calcPhidx(xnom, XS, YS)
%
% Calculate Elevation Jacobian
% By: Damennick Henry and Jordan Murphy
% Date: 12/19/18
% Description: Computes the sensing matrix, expected measurement
% perturbations and Kalman gain
% Inputs
%       xnom  - Nominal state at time t
%       XS    - Station X position
%       YS    - Statin Y position
%       range - Range measurement
% Outpus
%       rhodotdx - Range rate Jacobian
% =========================================
% =========================================

% Initialize
phidx = zeros(1,4);

% Calculate elements
phidx(1) = -(xnom(3)- YS)/(range^2);

% 2 = 0

phidx(3) = 1/(range^2);

% 4 = 0;