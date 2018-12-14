function plotinnov(t, in, sigs, fig)

if nargin > 3
   figure(fig.Number) 
else
    fig = figure;
end
% set(gcf, 'Position', get(0, 'Screensize'));

for idx = 1:6
    subplot(6,1,idx)
    hold on
    plot(t, in(idx,:), 'k')
    plot(t, 2*sigs(idx,:), 'r--')
    plot(t, -2*sigs(idx,:), 'r--')
end
end