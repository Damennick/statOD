function F = getJac(x, u)
% =========================================
% =========================================
%
% [Hf, dyhat, K] = sensMatrix(xnom, dxhat, Pminus, t, R)
%
% Calculate Gain
% By: Damennick Henry and Jordan Murphy
% Date: 12/19/18
% Description: Computes the Jacobian Matrix 
% Inputs
%       x      - Current estimated state
%       u      - gravitational parameter
% Outputs
%       F    - Jacobian Matrix
% =========================================
% =========================================

F = [0 1 0 0;
    (-u*(x(1)^2 + x(3)^2 - 3*x(1)^2))/((x(1)^2 + x(3)^2)^(5/2)) 0 3*u*x(1)*x(3)*((x(1)^2 +x(3)^2)^(-5/2)) 0;
    0 0 0 1;
    3*u*x(1)*x(3)*((x(1)^2 +x(3)^2)^(-5/2)) 0 (-u*(x(1)^2 + x(3)^2 - 3*x(3)^2))/((x(1)^2 + x(3)^2)^(5/2)) 0];



end