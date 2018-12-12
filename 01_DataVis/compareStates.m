function fig = compareStates(xl, xn, tl, tn)
% =========================================
% =========================================
%
% Compare States
% By: Damennick Henry and Jordan Murphy
% Date: 12/4/18
% Description: Plots linear and nonlinear states on the same graphs for
% comparison
% Inputs
%       xl - Linear states (4 x n)
%       xn - Nonlinear state (4 x m)
%       tl - Time vector for linear states (1 x n)
%       tn - Time vector for nonlinear states (1 x m)
% Outputs
%       fig - Handle for figure
%
% =========================================
% =========================================

% Plot first state
fig = figure;
subplot(4,1,1)
hold on
box on
plot(tl,xl(1,:),'r--')
plot(tn,xn(1,:),'k')
ylabel('$$X$$')
legend

% Plot second state
subplot(4,1,2)
hold on
box on
plot(tl,xl(2,:),'r--')
plot(tn,xn(2,:),'k')
ylabel('$$\dot{X}$$')

% Plot third state
subplot(4,1,3)
hold on
box on
plot(tl,xl(3,:),'r--')
plot(tn,xn(3,:),'k')
ylabel('$$Y$$')

% Plot fourth state
subplot(4,1,4)
hold on
box on
plot(tl,xl(4,:),'r--')
plot(tn,xn(4,:),'k')
ylabel('$$\dot{Y}$$')


end