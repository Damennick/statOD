function fig = performNISTest(NIS,y,alpha)
% =========================================
% =========================================
%
% fig = performNEESTest(chiVar, n, alpha)
%
% Perform NEES test
% By: Damennick Henry and Jordan Murphy
% Date: 12/19/18
% Description: Performs a chi squared test by plotting chi-squared samples
% with corresponding significance level bounds
% Inputs
%       NIS - Normalized Innovation Squared Error and the size of each NIS
%       value
%       alpha  - significance level
% Outputs
%       fig - Figure chi-squared test is performed on
%
% =========================================
% =========================================

fig = figure;
hold on
% Number of trials and number of chivars in a run
[N,K] = size(NIS);
% Sample mean
sampMean = mean(NIS);
% Collect measurement size
[~, nMeas] = size(y);
mesSize = zeros(1,nMeas);
for idx = 1:nMeas()
    yvec = y(:,idx);
    mesSize(idx) = length(yvec(~isnan(yvec)));
end
[sizes, ~, iMes] = unique(mesSize);
% Number of different vector sizes
nSize = length(sizes);
% Colormap for all sizes
map = colormap(lines(nSize));
inBounds = 0;
for idx = 1:nSize
    sizeIdx = find(iMes == idx)';
    plot(sizeIdx, sampMean(sizeIdx),'x', 'Color', map(idx,:))
    % Chi square bounds
    r1 = chi2inv(alpha/2, N*sizes(idx))/N;
    r2 = chi2inv(1 - (alpha/2), N*sizes(idx))/N;
    plot(r1*ones(size(1:K)), '--', 'Color', map(idx,:), 'LineWidth', 2, ...
        'HandleVisibility', 'off')
    plot(r2*ones(size(1:K)), '--', 'Color', map(idx,:), 'LineWidth', 2, ...
        'HandleVisibility', 'off')
    inBounds = inBounds + numel(find(sampMean(sizeIdx) > r1 & sampMean(sizeIdx) < r2));
end
legend('Measurement Size = 3', 'Measurement Size = 6')
disp(['In bounds: ' num2str(inBounds/K)])
end