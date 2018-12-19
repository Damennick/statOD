function fig = plotMeasurements(y, t, titletext)
% =========================================
% =========================================
%
% Plot Measurements
% By: Damennick Henry and Jordan Murphy
% Date: 12/4/18
% Description: Plots the measurements
% Inputs
%       y  - Vector of measuremetns for each station (12 x n)
%       t  - Time vector for measurements (1 x n)
% Outputs
%       fig - Handle for figure
%
% =========================================
% =========================================

fig = figure;
fig.Units = 'inches';
fig.Position = [3, 3, 8, 7];
hold on
map = colormap(lines(12));
for idx = 1:12
    subplot(3,1,1)
    ylabel('Range (km)')
    hold on
    box on
    plot(t,y((3*idx - 2),:),'x','Color', map(idx, :))
    
    subplot (3,1,2)
    ylabel('Range Rate (km/s)')
    hold on
    box on
    plot(t,y((3*idx - 1),:),'x','Color', map(idx, :))
    
    subplot (3,1,3)
    ylabel('Elevation (rad)')
    xlabel('Time (s)')
    hold on
    box on
    plot(t,y((3*idx),:),'x','Color', map(idx, :))
end

if nargin > 2
    suptitle(titletext) 
end
end