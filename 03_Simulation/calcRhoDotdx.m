function rhodotdx = calcRhoDotdx(xnom, XS, YS, XSdot, YSdot, range)
% =========================================
% =========================================
%
% [rhodotdx] = calcRhodx(xnom, XS, YS, range)
%
% Calculate Range Rate Jacobian
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
%       rhodotdx - Range rate Jacobian
% =========================================
% =========================================

% Initialize
rhodotdx = zeros(1,4);

% Calculate elements
rhodotdx(1) = ((xnom(2) - XSdot)/range)...
    - (xnom(1) - XS)*((xnom(1) - XS)*(xnom(2) - XSdot) + (xnom(3) - YS)*(xnom(4) - YSdot))/(range^3);

rhodotdx(2) = (xnom(1) - XS)/range;

rhodotdx(3) = ((xnom(4) - YSdot)/range) -...
    (xnom(3) - YS)*((xnom(1) - XS)*(xnom(2) - XSdot) + (xnom(3) - YS)*(xnom(4) - YSdot))/(range^3);

rhodotdx(4) = (xnom(3) - YS)/range;