function [Hf, dyhat, K] = sensingMatrix(xnom, dxhat, t,Pminus, R)
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
    H = [(xnom(1)-XS(jdx))/range, 0, (xnom(3) - YS(jdx))/range, 0];
    % Second row
    H = [H; ((xnom(2) - XSdot(jdx))*range - H(1)*((xnom(1) - ...
        XS(jdx))*(xnom(2) - XSdot(jdx)) + (xnom(3) ...
        - YS(jdx))*(xnom(4) - YSdot(jdx))))/(range^2),...
        ((xnom(4) - YSdot(jdx))*range - H(3)*((xnom(1) - ...
        XS(jdx))*(xnom(2) - XSdot(jdx)) + (xnom(3) - ...
        YS(jdx))*(xnom(4) - YSdot(jdx))))/(range^2),...
        (xnom(1) - XS(jdx))/range,...
        (xnom(3) - YS(jdx))/range];
    % Third row
    H = [H; -(xnom(3) - YS(jdx))/(range^2), 0, (xnom(1) - XS(jdx))/(range^2),0];
    dyhat((3*si(jdx))-2:3*si(jdx)) = H*dxhat(:);
    Hf = [Hf; H];
    if nargin > 3
        K = [K, Pminus*H'/(H*Pminus*H' + R)];
    end
end
end