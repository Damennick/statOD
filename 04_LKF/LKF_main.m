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
P0 = diag([1e-8, 1e-9, 1e-7, 1e-8]);
Sv = chol(P0,'lower');
load('orbitdeterm_finalproj_KFdata.mat')


%% Generate data
% dx0 = dx0hat + Sv*randn(4,1);
% xnom = [r0*cos(n*t); -r0*n*sin(n*t); r0*sin(n*t); r0*n*cos(n*t)];
% [xtrue, ytrue, xn, yn] = genTruth(t, dt, mu, x0 + dx0, gamma, Qtrue, Rtrue);
% dx = xtrue - xnom;
%% Implement LKF
Mes = [];
QLKF = diag([1e-5, 1e-7, 1e-5, 1e-7]);
for idx = 2:length(ydata)
    mes = ydata{idx};
    Mes = [Mes, mes];
    [~, nMes] = size(mes);
    yt = nan(36,1);
    for iMes = 1:nMes
        station = mes(4,iMes);
        yt((3*station)-2:3*station,1) = mes(1:3, iMes);
    end
    ytrue(:, idx-1) = yt;
end
xtrue = zeros(4, length(ytrue) + 1);
t = 0:dt:length(ytrue)*dt;
[xhat, yhat, dxhat, sigmas] = linearizedKF2(t, ytrue, xtrue, dx0hat, P0, mu, r0, dt, QLKF, gamma, Rtrue);

%% Data Visualization
% plotMeasurements(ytrue,t(2:end));
% fig = plotStateErrors(t, dx - dxhat,sigmas);
% subplot(4,1,1)
% legend('Error', '2$$\sigma$$ bounds')
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 4 4];
% figure(2);
% print('../07_Images/02_LKFAnalysis/stateErrors1000s', '-dpng', '-r200')
fig = plotStates(xhat,t, 'Estimated States');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 4];
print('../07_Images/02_LKFAnalysis/EstimatedStatesRealData', '-dpng', '-r200')
fig = plotStates(2*sigmas, t, 'Estimated 2$$\sigma$$');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 4];
print('../07_Images/02_LKFAnalysis/EstimatedSigmasRealData', '-dpng', '-r200')
% compareStates(xhat,xtrue,t,t);