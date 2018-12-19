function fig = performNEESTest(chiVar, n, alpha)
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
%       chiVar - Chi-squared samples
%       n      - Size of error
%       alpha  - significance level
% Outputs
%       fig - Figure chi-squared test is performed on
%
% =========================================
% =========================================

fig = figure;
hold on
box on
% Number of trials and number of chivars in a run
[N,K] = size(chiVar);
% Sample mean
sampMean = mean(chiVar);
plot(1:K,sampMean, 'kx')
% Chi square bounds
r1 = chi2inv(alpha/2, N*n)/N;
r2 = chi2inv(1 - (alpha/2), N*n)/N;
plot(r1*ones(size(1:K)), 'r--')
plot(r2*ones(size(1:K)), 'r--')
ylim([r1 - r1/2, r2 + r2/2])
disp(['In bounds: ' num2str(numel(find(sampMean > r1 & sampMean < r2))/K)])
xlabel('Time Step')
ylabel('NEES')
end