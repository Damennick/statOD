function [x, y, dx] = linOrbitSim(t, dx0, mu, r0,dt, xnonlin)
% =========================================
% =========================================
%
% Linear Orbit Sim
% By: Damennick Henry and Jordan Murphy
% Date: 12/4/18
% Description: Simulates linear stat OD system
%
% =========================================
% =========================================

% Mean motion
n = sqrt(mu/(r0^3));
% Nominal trajectory
xnom = [r0*cos(n*t); -r0*n*sin(n*t); r0*sin(n*t); r0*n*cos(n*t)];
ynom = getY(t(2:end),xnom(:,2:end), false);
dx = zeros(4,length(xnom));
dx(:,1) = dx0;
% Initialize state trajectories
x = xnom;
x(:,1) = xnom(:,1) + dx0;
% Input to state matrix
Bnom = [0 0; 1 0; 0 0; 0 1];
y = ynom;
% Simulate the system
for idx = 2:length(t)
    % Nominal dynamics matrix
    Anom = dynMatrix(t(idx-1),mu,n,r0);
    % State transition matrix
    F = eye(4) + Anom*dt;
    % Propogate state
    dx(:,idx) = F*dx(:,idx-1);
    x(:,idx) = xnom(:,idx) + dx(:,idx);
    [~, dy, ~] = sensingMatrix(xnom(:,idx),xnonlin(:,idx),dx(:,idx),t(idx));
    % Measurement
    y(:,idx-1) = ynom(:,idx-1) + dy;
end
end