% =========================================
% =========================================
%
% Linear vs Nonlinear Orbit Sim
% By: Damennick Henry and Jordan Murphy
% Date: 12/4/18
% Description: Compares linear and nonlinear stat OD system
%
% =========================================
% =========================================
clear all; close all
addpath('../03_Simulation')
addpath('../01_DataVis')
%% Parameters
% Earth's graviational parameter [km^3/s^2]
mu = 3.986e5;
% Nominal orbit radius [km]
r0 = 6678;
% Nominal orbit period [s]
T = 2*pi*sqrt((r0^3)/mu);
% Sampling time
dt = 10;
% Final simulation time
tf = 14000;
% Simulation time vector
tlin = 0:dt:tf;
% Initial state
x0 = [6678; 0; 0; r0*sqrt(mu/(r0^3))];
% Perturbation from initial state
dx0 = [0; 0.075; 0; -0.021];

%% Nonlinear system
tspan = tlin;
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[tnonlin,xnonlin] = ode45(@(t,x) nonlinOrbitSim(t,x,mu,zeros(2,1)),tspan,x0+dx0,options);
xnonlin = xnonlin';
ynonlin = getY(tnonlin,xnonlin, true);

%% Linear system
[xlin, ylin, dx] = linOrbitSim(tlin,dx0,mu,r0,dt,xnonlin);

%% Data Visualization
% Compare the states
fig = compareStates(xlin, xnonlin, tlin, tnonlin);
subplot(4,1,1)
legend('Nonlinear', 'Linear')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 4];
print('../07_Images/01_BasicSystemAnalysis/NLvL', '-dpng', '-r200')

% Compare measurements
fig = compareMes(ylin, ynonlin, tlin(2:end), tnonlin);
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 4];
print('../07_Images/01_BasicSystemAnalysis/Mescomp', '-dpng', '-r200')

% Plot perturbations
fig = plotStates(dx,tlin, 'Linearized Perturbations');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4 4];
print('../07_Images/01_BasicSystemAnalysis/perts', '-dpng', '-r200')

% Plot nonlinear measurements
fig = plotMeasurements(ynonlin,tnonlin, 'Nonlinear Measurements');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4];
print('../07_Images/01_BasicSystemAnalysis/nonlinMes', '-dpng', '-r200')

% Plot linearized measurements
fig = plotMeasurements(ylin,tlin(2:end), 'Linearized Measurements');
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4];
print('../07_Images/01_BasicSystemAnalysis/linMes', '-dpng', '-r200')