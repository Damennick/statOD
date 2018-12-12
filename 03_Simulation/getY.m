function y = getY(t, x)
% =========================================
% =========================================
%
% Get Measurements
% By: Damennick Henry and Jordan Murphy
% Date: 12/4/18
% Description: Calculates nonlinear measurements corresponding to the state
% x
%
% =========================================
% =========================================
for idx = 1:length(t)
    [XS, YS, XSdot, YSdot,satidx] = getStationPositions(x(:,idx),t(idx));
    % Number of stations visible
    [~, m] = size(YS);
    y(:,idx) = NaN(36,1);
    % Loop through satellite positions (jdx)
    for jdx = 1:m
        % Range for station
        range = sqrt((x(1,idx) - XS(jdx))^2 + (x(3,idx) - YS(jdx))^2);
        % Range rate
        rdot = (x(1,idx)-XS(jdx))*(x(2,idx) -XSdot(jdx)) + (x(3,idx)-YS(jdx))*(x(4,idx) -YSdot(jdx));
        rdot = rdot/range;
        % Elevation
        EL = atan2(x(3,idx) - YS(jdx), x(1,idx) - XS(jdx));
        y(3*satidx(jdx)-2:3*satidx(jdx),idx) = [range; rdot; EL];
    end
end
end