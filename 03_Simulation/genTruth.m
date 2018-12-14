function [xt, yt, xNoNoise, yNoNoise] = genTruth(t, dt, mu, x0, T, Q, R)
% =========================================
% =========================================
%
% Generate truth
% By: Damennick Henry and Jordan Murphy
% Date: 12/15/18
% Description: Generates ground truth states and measurements for stat OD
% system
% Inputs
%       t  - Time vector [1 x n]
%       dt - Time step [1]
%       mu - Gravitational parameter [1]
%       x0 - Initial state
%       T  - Noise to state matrix
%       Q  - Process noise covariance
%       R  - Measurement noise covariance
% Outputs
%       xt - Ground truth states [4 x n]
%       yt - Ground truth measurements [4 x n]
%       
%
% =========================================
% =========================================

%% True States
% Initialize array of true states
xt = zeros(4,length(t));
xt(:,1) = x0;
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[~, xNoNoise] = ode45(@(t,x) nonlinOrbitSim(t,x,mu,zeros(2,1)),t,xt(:,1),options);
xNoNoise = xNoNoise';
for idx = 1:length(t)-1
    % Time span for simulation
    tspan = t(idx):dt:t(idx + 1);
    % Cholesky factorization of covariance matrix
    Sv = chol(Q,'lower');
    % Sample noise
    wk = Sv*randn(2,1);
    % Run the CT nonlinear simulation
    [~,x] = ode45(@(t,x) nonlinOrbitSim(t,x,mu,wk),tspan,xt(:,idx),options);
    % Append to state vectors
    xt(:,idx+1) = x(end,:)';
end

%% True Measuremetns
% Get measurements
yt = getY(t(2:end),xt(:,2:end), true);
yNoNoise = getY(t(2:end), xNoNoise(:,2:end),true);
% Add noise
R = repmat({R}, 1, 12);
R = blkdiag(R{:});
yt = addGaussNoise(yt, R, eye(36));
