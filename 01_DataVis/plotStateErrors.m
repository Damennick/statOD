function fig = plotStateErrors(t, e, sigmas, fig)
% =========================================
% =========================================
%
% Plot State Errors
% By: Damennick Henry and Jordan Murphy
% Date: 12/4/18
% Description: Plots linear and nonlinear states on the same graphs for
% comparison
% Inputs
%       t - Time vector for states (1 x n)
%       e -  State error (4 x n)
%       sigmas - Standard deviation of state
% Outputs
%       fig - Handle for figure
%
% =========================================
% =========================================

% Plot first state
if nargin > 3
    figure(fig.Number)
    clf
else
    fig = figure;
end
subplot(4,1,1)
hold on
box on
plot(t,e(1,:),'k')
plot(t,2*sigmas(1,:),'r--')
plot(t,-2*sigmas(1,:),'r--', 'HandleVisibility','off')
ylabel('$$X$$ (km)')
legend

% Plot second state
subplot(4,1,2)
hold on
box on
plot(t,e(2,:),'k')
plot(t,2*sigmas(2,:),'r--')
plot(t,-2*sigmas(2,:),'r--', 'HandleVisibility','off')
ylabel('$$\dot{X}$$ (km/s)')

% Plot third state
subplot(4,1,3)
hold on
box on
plot(t,e(3,:),'k')
plot(t,2*sigmas(3,:),'r--')
plot(t,-2*sigmas(3,:),'r--', 'HandleVisibility','off')
ylabel('$$Y$$ (km)')

% Plot fourth state
subplot(4,1,4)
hold on
box on
plot(t,e(4,:),'k')
plot(t,2*sigmas(4,:),'r--')
plot(t,-2*sigmas(4,:),'r--', 'HandleVisibility','off')
ylabel('$$\dot{Y}$$ (km/s)')
xlabel('Time (s)')

end