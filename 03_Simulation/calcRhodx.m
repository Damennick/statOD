function [rhodx] = calcRhodx(xnom, XS, YS, range)
% =========================================
% =========================================
%
% [rhodx] = calcRhodx(xnom, XS, YS, range)
%
% Calculate Range Jacobian
% By: Damennick Henry and Jordan Murphy
% Date: 12/19/18
% Description: Computes the sensing matrix, expected measurement
% perturbations and Kalman gain
% Inputs
%       xnom  - Nominal state at time t
%       XS    - Station X position
%       YS    - Statin Y position
%       XSdot - Station x velocity
%       YSdot - Station y velocity
%       range - Measurement
% Outpus
%       rhodx - Range jacobian
% =========================================
% =========================================

% Initialize
rhodx = zeros(1,4);

% Compute all elements
rhodx(1) = (xnom(1)-XS)/range;

% 2 = 0

rhodx(3) = (xnom(3) - YS)/range;

% 4 = 0
end