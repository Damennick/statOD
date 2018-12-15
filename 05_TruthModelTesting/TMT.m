% =========================================
% =========================================
%
% Truth Model Testing - Main
% By: Damennick Henry and Jordan Murphy
% Date: 12/19/18
% Description: Performs Truth Model Testing
%
% =========================================
% =========================================
clear all; close all
addpath('..')
addpath('../03_Simulation')
addpath('../01_DataVis')
%addpath('../04_LKF')
addpath('../06_EKF') 

%% Parameters
% Model parameters
% ------------------------------------------------
% Earth's graviational parameter [km^3/s^2]
mu = 3.986e5;
% Nominal orbit radius [km]
r0 = 6678;
% Nominal orbit period [s]
T = 2*pi*sqrt((r0^3)/mu);
% Mean motion for nominal trajectory
n = sqrt(mu/(r0^3));
% Sampling time
dt = 10;
% Noise to state matrix
gamma = [0 0; 1 0; 0 0; 0 1];

Q_EKF = .95 * diag([1E-8 1E-7 1E-8 1E-7]);
% ------------------------------------------------

% Simulation parameters
% ------------------------------------------------
% Final simulation time
tf = 1000;
% Simulation time vector
t = 0:dt:tf;
% Nominal trajectory
xnom = [r0*cos(n*t); -r0*n*sin(n*t); r0*sin(n*t); r0*n*cos(n*t)];
% Number of trials
N = 1000;
% ------------------------------------------------

% Filter parameters
% ------------------------------------------------
% Flag that determines which filter to use: L = LKF, E = EKF
filter = 'E';
% Nominal beginning state
x0 = [6678; 0; 0; r0*sqrt(mu/(r0^3))];
% Perturbation from initial state
dx0hat = [0; 0; 0; 0];
% Initial covariance guess
P0 = diag([1e-8, 1e-9, 1e-7, 1e-8])
% Initial covarianc guess Cholesky factorization
SvP0 = chol(P0);
% Kalman Filter process noise covariance
load('orbitdeterm_finalproj_KFdata.mat')
QLKF = diag([1e-8, 1e-7, 1e-8, 1e-7]);
% ------------------------------------------------



%% TMT
allNEES = zeros(N, length(t)-1);
allNIS  = zeros(N, length(t)-1);
% figState = figure;
% figMes = figure;
for idx = 1:N
    % Generate Data
    % ------------------------------------------------
    dx0 = dx0hat + SvP0*randn(4,1);
    [xtrue, ytrue, xn, yn] = genTruth(t, dt, mu, x0 + dx0, gamma, Qtrue, Rtrue);
    dx = xtrue - xnom;
    % ------------------------------------------------
    
    % Run Filternumel(find(sampMean > r1 & sampMean < r2))
    % ------------------------------------------------
    if filter == 'L'
        [xhat, yhat, dxhat, sigmas, NEES, NIS, innov, mesSigs] = ...
            linearizedKF2(t, ytrue, xtrue, dx0hat, P0, mu, r0, dt, QLKF, gamma, Rtrue);
%         plotinnov(t(2:end),innov,mesSigs, figMes)
%         shg
%         plotStateErrors(t, dx - dxhat,sigmas,figState);
%         shg
        allNEES(idx,:) =  NEES;
        allNIS(idx,:) = NIS;
    elseif filter == 'E'
        [xhat, sigmas, NEES, NIS, inn, mesSigmas] = ...
            extendedKF(t, xtrue, ytrue, dx0hat, P0, mu, r0, dt, Q_EKF, gamma, Rtrue);
        allNEES(idx,:) =  NEES;
        allNIS(idx,:) = NIS;
    end
    % ------------------------------------------------
end

%% Perform chi-square tests
% Perform NEES testing
performNEESTest(allNEES,4,0.05);
% Perform NIS test
performNISTest(allNIS, ytrue, 0.05);