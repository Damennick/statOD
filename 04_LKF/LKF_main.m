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
rng(23)
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
% Final simulation time
tf = 100;
% Simulation time vector
t = 0:dt:tf;
% Initial state
x0 = [6678; 0; 0; r0*sqrt(mu/(r0^3))];
% Perturbation from initial state
dx0 = [0; 0; 0; 0];
% Noise to state matrix
gamma = [0 0; 1 0; 0 0; 0 1];
load('orbitdeterm_finalproj_KFdata.mat')


%% Generate data
xnom = [r0*cos(n*t); -r0*n*sin(n*t); r0*sin(n*t); r0*n*cos(n*t)];
[xtrue, ytrue] = genTruth(t, dt, mu, x0 + dx0, gamma, Qtrue, Rtrue);
dx = xtrue - xnom;
%% Implement LKF
% Initial covariance guess
P0 = diag([10, 1, 10, 1]);
dx0hat = dx0;
[xhat, yhat, dxhat, sigmas] = linearizedKF2(t, ytrue, dx0hat, P0, mu, r0, dt, Qtrue, gamma, Rtrue);

%% Data Visualization
plotMeasurements(yhat,t(2:end))
plotStateErrors(t, dx - dxhat,sigmas);
plotStates(dxhat,t)
compareStates(xhat,xtrue,t,t);