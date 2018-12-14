function [XS, YS, XSdot, YSdot,stationidx] = getStationPositions(x,t,checkVis)
% =========================================
% =========================================
%
% Get Station Positions
% By: Damennick Henry and Jordan Murphy
% Date: 12/4/18
% Description: Calculates the position of the visible stations
% Inputs
%       x        - State of spacecraft
%       t        - Time of state
%       checkVis - Flag for whether or not to check visibility
%
% =========================================
% =========================================

% Radius of Earth
Re = 6378;
% Earth rotation rate
we = 2*pi/86400;
% Initial angular positions
theta0 = ((1:12) -1)*pi/6;
% Station positions
XS = Re*cos(we*t + theta0);
YS = Re*sin(we*t + theta0);
XSdot = -Re*we*sin(we*t + theta0);
YSdot = Re*we*cos(we*t + theta0);
% Elevation angles
EL = atan2((x(3) - YS),(x(1) - XS));
% Angle limits
thetai = atan2(YS,XS);
EL = EL - thetai;
% Return only visible satellites
if checkVis
    stationidx = find(abs(EL) <= pi/2 | abs(EL) >= 3*pi/2);
    % if isempty(satidx)
    %    debug = 1; 
    %    figure
    %    plot(XS, YS)
    %    hold on
    %    plot([XS(4) x(1)],[YS(4) x(3)],'r.-')
    % end
    XS = XS(stationidx);
    YS = YS(stationidx);
    XSdot = XSdot(stationidx);
    YSdot = YSdot(stationidx);
end
end