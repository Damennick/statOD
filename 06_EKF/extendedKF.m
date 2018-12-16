function [xhatplus, sigmas, NEES, NIS, inn, mesSigmas] = extendedKF(t, xtrue, y, dx0, P0, mu, r0, dt, Q, T, R)
% =========================================
% =========================================
%
% Extended KF 
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
%CT noise for state
wk = [0;0];
% Nominal state
xnom = [r0*cos(n*t); -r0*n*sin(n*t); r0*sin(n*t); r0*n*cos(n*t)];
% Peturbed state prediction
xhatplus = zeros(4, length(t));
xhatplus(:,1) = dx0 + xnom(:,1);
% Estimated measurements
yhat = zeros(36,length(t)-1);
% Standard deviations
sigmas = xhatplus;
% Initial sigmas
Pplus = P0;
sigmas(:,1) = sqrt(diag(Pplus));
% Measurment standard deviations
mesSigmas = nan(6,length(t)-1);

NEES = zeros(1,length(t)-1);
NIS  = zeros(1,length(t)-1);
inn = nan(6,length(t)-1);

%% Filter
for idx = 2:length(t)
    % Prediction Update
    % ------------------------------------------------------
    % Update xhatminus using integration
    [~,xhatminus] = ode45(@(t,x) nonlinOrbitSim(t,x,mu,wk), [0 dt], xhatplus(:,idx-1));
    xhatminus = xhatminus(end,:)';
    % State transition matrix
    F = eye(4) + dt * getJac(xhatplus(:,idx-1), mu);
    % Covariance prediction
    Pminus = F*Pplus*F' + Q;

    % ------------------------------------------------------
    
    % Measurement Update
    % ------------------------------------------------------
    % Compute yhatminus
    yhatminus = getY(t(idx), xhatminus, false);
    % Compute Sensing Matrix
    H = calcH(xhatminus, t(idx));
    % Fix Sensing Matrix so that it corresponds to when both measurements are
    % seen
    H = H(~isnan(y(:,idx-1)),:);
    % Measurement perturbation error
    ey = y(:,idx-1) - yhatminus;
    ey = ey(~isnan(ey));
    if isempty(ey)
        K = zeros(4,3);
        H = zeros(3,4);
        ey = zeros(3,1);
    end
    %Gain Matrix
    Rblk = repmat({R}, 1, length(ey)/3);
    Rblk = blkdiag(Rblk{:});
    K = Pminus * H' / (H * Pminus * H' + Rblk);
   
    % State measurement update
    xhatplus(:,idx) = xhatminus + K * ey;

    % Covariance measurement update
    Pplus = (eye(4) - K*H) * Pminus * (eye(4) - K*H)' + K*Rblk*K';
    % Standard deviation
    sigmas(:,idx) = sqrt(diag(Pplus));
    % S matrix
    S = H*Pminus*H' + Rblk;
    % Measurement Sigmas
    mesSigmas(1:length(ey),idx-1) = sqrt(diag(S));
    [NEES(idx-1), NIS(1,idx-1)] = calcNEESNIS(xtrue(:,idx) - xhatplus(:,idx), ...
        ey, Pplus, S);
    % ------------------------------------------------------
end
end