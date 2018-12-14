function Hf = calcH(xnom, t)
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
xhat = xnom;
% Station positions
[XS, YS, XSdot, YSdot] = getStationPositions(xhat, t, false);
% Number of stations visible
[~, m] = size(YS);
if m == 2
    debug = 1;
end

Hf = [];

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
    Hf = [Hf; H];

end
end