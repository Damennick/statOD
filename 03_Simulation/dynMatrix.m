function A = dynMatrix(t, mu, n, r0)
% =========================================
% =========================================
% A = dynMatrix(t, mu, n, r0)
% Dynamics Matrix
% By: Damennick Henry and Jordan Murphy
% Date: 12/19/18
% Description: Computes the dynamics matrix at time t
% Inputs
%       t  - Time to evaluate dynamics matrix
%       mu - Gravitational parameter
%       n  - Mean motion for nominal trajectory
%       r0 - Nominal radius
% Outpus
%       A - Dynamics matrix
%
% =========================================
% =========================================
A = [0, 1, 0, 0; (-mu + 3*mu*((cos(n*t))^2))/(r0^3), 0 ...
       (3*mu*sin(2*n*t))/(2*(r0^3)), 0;...
       0, 0, 0, 1; (3*mu*sin(2*n*t))/(2*(r0^3))...
       0, (-mu + 3*mu*((sin(n*t))^2))/(r0^3), 0];


end