function fig = compareMes(yl, yn, tl, tn)
% =========================================
% =========================================
%
% Compare Measurements
% By: Damennick Henry and Jordan Murphy
% Date: 12/4/18
% Description: Plots linear and nonlinear measurements of first three
% stations on the same graphs for comparison
% Inputs
%       yl - Linear measurement (2 x n)
%       yn - Nonlinear measurement (2 x m)
%       tl - Time vector for linear measurements (1 x n)
%       tn - Time vector for nonlinear measurements (1 x m)
% Outputs
%       fig - Handle for figure
%
% =========================================
% =========================================

fig =figure;
fig.Units = 'inches';
fig.Position = [3, 3, 7, 7];
stcount = 0;
% Loo
for idx = 1:9
   subplot(3,3,idx)
   hold on
   box on
   if mod(idx,3) == 1
       stcount = stcount + 1;
      ylabel(['Station ' num2str(stcount)])
   end
   plot(tl,yl(idx,:),'rx')
   plot(tn,yn(idx,:),'k*')
   if idx == 1
       title('Range [km]')
       legend('Linear', 'Nonlinear')
   end
   if idx == 2
       title('Range Rate [km/s]')
   end
   if idx == 3
       title('Elevation [rad]')
   end
end

