function fig = plotStates(x, t, titletext)
% =========================================
% =========================================
%
% Plot Perturbations
% By: Damennick Henry and Jordan Murphy
% Date: 12/4/18
% Description: Plots the linearized perturbations
% Inputs
%       dx - Linear perturbations (4 x n)
%       t  - Time vector for linear measurements (1 x n)
% Outputs
%       fig - Handle for figure
%
% =========================================
% =========================================

% Plot first perturbation
fig = figure;
subplot(4,1,1)
hold on
box on
plot(t,x(1,:),'k')
ylabel('$$X$$ (km)')

% Plot second perturbation
subplot(4,1,2)
hold on
box on
plot(t,x(2,:),'k')
ylabel('$$ \dot{X}$$ (km/s)')

% Plot third perturbation
subplot(4,1,3)
hold on
box on
plot(t,x(3,:),'k')
ylabel('$$ Y$$ (km)')

% Plot fourth perturbation
subplot(4,1,4)
hold on
box on
plot(t,x(4,:),'k')
ylabel('$$ \dot{Y}$$ (km/s)')
xlabel('Time (s)')

if nargin > 2
    h = suptitle('');
    h.Interpreter = 'latex';
    h.String = titletext;
end
end