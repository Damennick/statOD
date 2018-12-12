function fig = plotMeasurements(y,t)
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
    hold on
    plot(t,y((3*idx - 2),:),'o','Color', map(idx, :))
    
    subplot (3,1,2)
    hold on
    plot(t,y((3*idx - 1),:),'o','Color', map(idx, :))
    
    subplot (3,1,3)
    hold on
    plot(t,y((3*idx),:),'o','Color', map(idx, :))
end

end