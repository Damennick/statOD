function [xhat, yhat, dx, sigmas] = linearizedKF(t, y, dx0, P0, mu, r0, dt, Q, T, R)
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

T = dt*T;
% Mean motion
n = sqrt(mu/(r0^3));
% Nominal trajectory
xnom = [r0*cos(n*t); -r0*n*sin(n*t); r0*sin(n*t); r0*n*cos(n*t)];
ynom = getY(t(2:end),xnom(:,2:end));
dx = zeros(4,length(xnom));
dx(:,1) = dx0;
% Initialize state estimates
xhat = xnom;
xhat(:,1) = xnom(:,1) + dx0;
% Initialize measurement estimates
yhat = ynom;
% Measurement residuals
dy = y - ynom;
% Initial error covariance
PPrev = P0;
% Standard deviations
sigmas = dx;
sigmas(:,1) = sqrt(diag(P0));
% Simulate the system
for idx = 2:length(t)
    % Prediction Step
    % --------------------------------------------
    % Nominal dynamics matrix
    Anom = [0, 1, 0, 0; (-mu + 3*mu*((cos(n*t(idx-1)))^2))/(r0^3), 0 ...
       (3*mu*sin(2*n*t(idx-1)))/(2*(r0^3)), 0;...
       0, 0, 0, 1; (3*mu*sin(2*n*t(idx-1)))/(2*(r0^3))...
       0, (-mu + 3*mu*((sin(n*t(idx-1)))^2))/(r0^3), 0];
    % State transition matrix
    F = eye(4) + Anom*dt;
    % Propogate state
    dxminus = F*dx(:,idx-1);
    xhat(:,idx) = xnom(:,idx) + dx(:,idx);
    % Predict covariance
    Pminus = F*PPrev*F' + T*Q*T';
    % --------------------------------------------
    
    % Meaurement Update
    % --------------------------------------------
    % Create sensing matrix and Kalman Gain matrix
    K = [];
    Ht = [];
    % Station positions
    [XS, YS, XSdot, YSdot, si] = getStationPositions(xhat(:,idx-1), t(idx-1));
    % Number of stations visible
    [~, m] = size(YS);
    if m == 2
        debug = 1;
    end
    dyhat = NaN(36,1);
    for jdx = 1:m
        % Range for station
        range = sqrt((xnom(1,idx-1) - XS(jdx))^2 + (xnom(3,idx-1) - YS(jdx))^2);
        % First row of sensing matrix for station
        Hi = [(xnom(1,idx-1)-XS(jdx))/range, 0, (xnom(3,idx-1) - YS(jdx))/range, 0];
        % Second row
        Hi = [Hi; ((xnom(2,idx-1) - XSdot(jdx))*range - Hi(1)*((xnom(1,idx-1) - ...
            XS(jdx))*(xnom(2,idx-1) - XSdot(jdx)) + (xnom(3,idx-1) ...
            - YS(jdx))*(xnom(4,idx-1) - YSdot(jdx))))/(range^2),...
            (xnom(1,idx-1) - XS(jdx))/range,...
            ((xnom(4,idx-1) - YSdot(jdx))*range - Hi(3)*((xnom(1,idx-1) - ...
            XS(jdx))*(xnom(2,idx-1) - XSdot(jdx)) + (xnom(3,idx-1) - ...
            YS(jdx))*(xnom(4,idx-1) - YSdot(jdx))))/(range^2),...
            (xnom(3,idx-1) - YS(jdx))/range];
        % Third row
        Hi = [Hi; -(xnom(3,idx-1) - YS(jdx))/(range^2), 0, (xnom(1,idx-1) - XS(jdx))/(range^2),0];
        dyhat((3*si(jdx))-2:3*si(jdx)) = Hi*dx(:,idx-1);
        Ht = [Ht;Hi];
        % Kalman Gain
        K = [K, (Pminus*Hi')/(Hi*Pminus*Hi' + R)];
    end
    % Measurement error
    ey = dy(:,idx-1) - dyhat;
    ey = ey(~isnan(ey));
    % Only use measurements present in truth and expected dy
    K = K(:,1:length(ey));
    Ht = Ht(1:length(ey),:);
    % State measurement update
    dx(:,idx) = dxminus + K*ey;
    % Covariance measurement update
    PPrev = (eye(4) - K*Ht)*Pminus;
    % Ensure matrix is symmetric
    if ~issymmetric(PPrev)
        PPrev = (PPrev + PPrev')/2;
    end
    % --------------------------------------------
    % Measurement estimate
    yhat(:,idx-1) = ynom(:,idx-1) + dyhat;
    % Standard deviation
    sigmas(:,idx) = sqrt(diag(PPrev));
    if ~isreal(sigmas(:,idx))
        debug = 1;
    end
end
end