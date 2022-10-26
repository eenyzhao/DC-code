function  calcuTau = id(refQ,fs)
% inverse dynamics
%% subject anatomial data
bodyMass = 78; %kg
bodyHeight = 178;% cm
lofHand = 0.08;%m

m = bodyMass*0.0065;   % Kg
g = 9.81;
%CoM : l = 4.11 + 0.026*bodyweight(kg) + 0.033*height(cm) from 3rd dactylion
l = 0.08;   % hand length

% Inertia parameters
J = (-13.68 + 0.088*bodyMass + 0.092*bodyHeight)/10000;
J = J + m*lofHand^2;
B = 0.3; % damping parameter

dt = 1/fs;
dtheta  = diff(refQ)/dt;
ddtheta = diff(dtheta)/dt;

% spline match the motion with EMG
%data cubic spline
duration = (length(refQ)-1)/fs;
t = 0:1/fs:duration;

t_dtheta = dt:dt:duration;
t_ddtheta = dt*2:dt:duration;

dtheta= spline(t_dtheta , dtheta , t)';
ddtheta = spline(t_ddtheta , ddtheta , t)';

calcuTau = J*ddtheta + 0.3*dtheta + 0.6 * refQ;

end