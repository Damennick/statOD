% =========================================
% =========================================
%
% Linearized Kalman Filter - Main
% By: Damennick Henry and Jordan Murphy
% Date: 12/15/18
% Description: Runs a linearized Kalman filter on a stat od system
%
% =========================================
% =========================================
clear all; close all
addpath('..')
addpath('../03_Simulation')
addpath('../01_DataVis')
addpath('../05_TruthModelTesting')
% rng(23)
%% Parameters
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
% Final simulation time [s]
tf = 1000;
% Simulation time vector
t = 0:dt:tf;
% Initial state
x0 = [6678; 0; 0; r0*sqrt(mu/(r0^3))];
% Perturbation from initial state
dx0hat = [0; 0; 0; 0];
% Noise to state matrix
gamma = [0 0; 1 0; 0 0; 0 1];
% Initial covariance guess
P0 = diag([100, 1, 100, 1])/10000000000
Sv = chol(P0,'lower');
load('orbitdeterm_finalproj_KFdata.mat')


%% Generate data
dx0 = dx0hat + Sv*randn(4,1);
xnom = [r0*cos(n*t); -r0*n*sin(n*t); r0*sin(n*t); r0*n*cos(n*t)];
[xtrue, ytrue, xn, yn] = genTruth(t, dt, mu, x0 + dx0, gamma, Qtrue, Rtrue);
dx = xtrue - xnom;
%% Implement LKF

QLKF = diag([1e-5, 1e-7, 1e-5, 1e-7]);
[xhat, yhat, dxhat, sigmas] = linearizedKF2(t, ytrue, xtrue, dx0hat, P0, mu, r0, dt, QLKF, gamma, Rtrue);

%% Data Visualization
plotMeasurements(ytrue,t(2:end));
plotStateErrors(t, dx - dxhat,sigmas);
plotStates(dxhat,t);
compareStates(xhat,xtrue,t,t);