function [NEES, NIS] = calcNEESNIS(ex, ey, Pk, Sk)
% =========================================
% =========================================
%
% [NEES, NIS] = calcNEESNIS(xt, xh, yt, yh, Pk, R, H)
%
% Calculate Normalized Estimation Error Squared
% By: Damennick Henry and Jordan Murphy
% Date: 12/19/18
% Description: Calculates NEES
% Inputs
%       xt - True states
%       xh - Estimated states
%       yt - True measurements
%       yh - Estimated measurements
%       Pk - Cell array of error covariances
%       R  - Measurement covariance
%       H  - Sensing matrix
% Outputs
%       NEES - Normalized Estimaiton Error Squared
%
% =========================================
% =========================================

% NEES
NEES = ex'*(Pk\ex);

if isempty(ey)
    NIS = nan(1);
else
    % NIS
    NIS = ey'*(Sk\ey);
end

end