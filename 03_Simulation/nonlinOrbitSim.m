function xdot = nonlinOrbitSim(t,x,mu,wk)
% =========================================
% =========================================
%
% Nonlinear Orbit Sim
% By: Damennick Henry and Jordan Murphy
% Date: 12/4/18
% Description: Simulates linear stat OD system
%
% =========================================
% =========================================

xdot = zeros(4,1);
xdot(1) = x(2);
xdot(2) = (-mu*x(1))/(sqrt(x(1)^2 + x(3)^2)^3) + wk(1);
xdot(3) = x(4);
xdot(4) = (-mu*x(3))/(sqrt(x(1)^2 + x(3)^2)^3) + wk(2);
end