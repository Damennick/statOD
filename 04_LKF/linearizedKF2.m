function [xhat, yhat, dxhat, sigmas, NEES, NIS, inn, mesSigmas] = linearizedKF2(t, y, x, dx0, P0, mu, r0, dt, Q, T, R)
% =========================================
% =========================================
%
% Linearized KF 
% By: Damennick Henry and Jordan Murphy
% Date: 12/19/18
% Description: Implements and LKF for the stat OD system
% Inputs:
%       t   - Time Vector
%       y   - Vector of measurements
%       dx0 - Initial perturbation guess
%       P0  - Initial covariance guess
%       mu  - Gravitational parmaeter
%       r0  - Orbit radius
%       dt  - Time step
%       Q   - Process noise covariance
%       T   - Noise to state matrix
%       R   - Measurement noise covariance
% Outputs:
%       xhat   - Estimated state
%       yhat   - Estimated Measurements
%       dx     - Estimated perturbations
%       sigmas - Standard deviations of each state
% =========================================
% =========================================

%% Initialization
% Mean motion for nominal trajectory
n = sqrt(mu/(r0^3));
% Process noise to state matrix
Omega = dt*T;
% Nominal state
xnom = [r0*cos(n*t); -r0*n*sin(n*t); r0*sin(n*t); r0*n*cos(n*t)];
% Peturbed state prediction
dxhat = zeros(4, length(t));
dxhat(:,1) = dx0;
% State estimate
xhat = dxhat + xnom;
% True peturbed state
dx = x - xnom;
% Nominal Measurements
ynom = getY(t(2:end),xnom(:,2:end), false);
% Peturbed measurements
dy = y - ynom;
% Estimated measurements
yhat = zeros(36,length(t)-1);
% Standard deviations
sigmas = dxhat;
% Measurment standard deviations
mesSigmas = nan(6,length(t)-1);
% Initial sigmas
PPrev = P0;
Pk = cell(1,length(t));
Pk{1} = P0;
sigmas(:,1) = sqrt(diag(PPrev));
NEES = zeros(1,length(t)-1);
NIS  = zeros(1,length(t)-1);
inn = nan(6,length(t)-1);
%% Filter
for idx = 2:length(t)
    % Prediction Update
    % ------------------------------------------------------
    % Compute dynamics matrix
    Anom = dynMatrix(t(idx), mu, n, r0);
    % State transition matrix
    F = eye(4) + Anom*dt;
    % State perturbation prediction
    dxminus = F*dxhat(:,idx-1);
    Pminus = F*PPrev*F' + Q;
    % ------------------------------------------------------
    
    % Measurement Update
    % ------------------------------------------------------
    if t(idx) == 4670
       debug = 1; 
    end
    % Compute sensing matrix, estimated measurements, and Kalman Gain
    [H, dyhat, K] = calcGain(xnom(:,idx),dxminus,Pminus,t(idx),R);
    yhat(:,idx-1) = ynom(:,idx-1) + dyhat;
    % Measurement perturbation error
    ey = dy(:,idx-1) - dyhat;
    ey = ey(~isnan(ey));
    inn(1:length(ey),idx-1) = ey;
    % Fix Kalman gain so that it corresponds to to stations that see
    % satellite
    K = K(:,~isnan(dy(:,idx-1)));
    % Same for the sensing matrix
    H = H(~isnan(dy(:,idx-1)),:);
    if isempty(ey)
        K = zeros(4,3);
        H = zeros(3,4);
        ey = zeros(3,1);
        disp(['No measurement at t = ' num2str(t(idx))])
    else
        debug = 1;
    end
    % State measurement update
    dxhat(:,idx) = dxminus + K*ey;
    Rblk = repmat({R}, 1, length(ey)/3);
    Rblk = blkdiag(Rblk{:});
    % Covariance measurement update
    PPrev = (eye(4) - K*H)*Pminus*(eye(4) - K*H)' + K*Rblk*K';
    if ~isreal(sqrt(diag(PPrev)))
        debug = 1;
    end
    % Standard deviation
    sigmas(:,idx) = sqrt(diag(PPrev));
    Pk{idx} = PPrev;
    xhat(:,idx) = xnom(:,idx) + dxhat(:,idx);
    % S matrix
    S = H*Pminus*H' + Rblk;
    % Measurement Sigmas
    mesSigmas(1:length(ey),idx-1) = sqrt(diag(S));
    [NEES(idx-1), NIS(1,idx-1)] = calcNEESNIS(dx(:,idx) - dxhat(:,idx), ...
        ey, PPrev, S);
    % ------------------------------------------------------
end
end