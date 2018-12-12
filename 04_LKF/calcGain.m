function [Hf, dyhat, K] = calcGain(xnom, dxhat, Pminus, t, R)
% =========================================
% =========================================
%
% [Hf, dyhat, K] = sensMatrix(xnom, dxhat, Pminus, t, R)
%
% Calculate Gain
% By: Damennick Henry and Jordan Murphy
% Date: 12/19/18
% Description: Computes the sensing matrix, expected measurement
% perturbations and Kalman gain
% Inputs
%       dxhat  - Estimated state perturbation at time t
%       xnom   - Nominal state at time t
%       Pminus - Predicted state error covariance 
%       t      - time
% Outpus
%       Hf    - Sensing matrix
%       dyhat - Measurement perturbations
%       K     - Kalman gain
% =========================================
% =========================================

% Estimated state
xhat = xnom + dxhat;
% Station positions
[XS, YS, XSdot, YSdot, si] = getStationPositions(xhat, t);
% Number of stations visible
[~, m] = size(YS);
if m == 2
    debug = 1;
end
% Measurement perturbations
dyhat = NaN(36,1);
% Final Sensing matrix
Hf = [];
% Kalman gain matrix
K = [];
% Create sensing matrix
for jdx = 1:m
    % Range for station
    range = sqrt((xnom(1) - XS(jdx))^2 + (xnom(3) - YS(jdx))^2);
    % First row of sensing matrix for station
    H = calcRhodx(xnom,XS(jdx),YS(jdx),range);
    % Second row
    H = [H; calcRhoDotdx(xnom,XS(jdx), YS(jdx), XSdot(jdx), YSdot(jdx), range)];
    % Third row
    H = [H; calcPhidx(xnom,XS(jdx),YS(jdx),range)];
    dyhat((3*si(jdx))-2:3*si(jdx)) = H*dxhat(:);
    Hf = [Hf; H];
    if nargin > 3
        K = [K, Pminus*H'/(H*Pminus*H' + R)];
    end
end
end