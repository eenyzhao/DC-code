function [Lmt,Lm,ddtheta,fv]= thetaPred(theta,dtheta,a,para,vmt)
%% function used to compute the tau at current time step
% fv=zeros(1,5);
%% Declare hand parameters;
bodyWeight = 88;             % Kg
bodyHeight = 178;            % cm
% m = bodyWeight*0.0065;       % --> hand mass (Kg)
m = 0.589;
g = 9.81;                    % --> gravity coefficient (m/s^2)
%CoM : l = 4.11 + 0.026*bodyweight(kg) + 0.033*height(cm) from 3rd dactylion
l = 0.08;                    % --> hand length (m)
J = (-13.68 + 0.088*bodyWeight + 0.092*bodyHeight)/10000;% --> inertia of hand
inertia = J + m*l^2;         % --> inertia around joint rotation center
%% allocate the parameter
lmo = para(1,1:5);          % optimal muscle fibre length
fmo = para(1,6:10);         % maximum isometric force
lts = para(1,11:15);        % tendon slack length
phi = [0.05 0.2 0.01 0.16 0.06];  % Pennation angle
mtScal = para(1,16:20);    % mtScale
% dampingC = parameter(22);% damping coefficient
C = 0.3;    %damping coefficient
K = 0.6;   %stiffness coefficient
%% Compute the joint torque using musculoskeletal model
%---> Step one: with initial state of Q, calculate the moment arm and change
% of muscle-tendon length
[Lmt,MA] = getMTUandMA(theta);
Lmt = Lmt.*mtScal;
MA  = MA.*mtScal;
% --->Step two: Calculate the current muscle fibre length
Lm = sqrt((lmo.*sin(phi)).^2 + (Lmt - lts).^2);
% --->Step Three: Calculate the joint moment
[fmt(1),fv(1)]  = updateT(Lm(1),lmo(1),phi(1),fmo(1),vmt(1),mtScal(1),a(1));
[fmt(2),fv(2)]  = updateT(Lm(2),lmo(2),phi(2),fmo(2),vmt(3),mtScal(2),a(2));
[fmt(3),fv(3)]  = updateT(Lm(3),lmo(3),phi(3),fmo(3),vmt(3),mtScal(3),a(3));
[fmt(4),fv(4)]  = updateT(Lm(4),lmo(4),phi(4),fmo(4),vmt(4),mtScal(4),a(4));
[fmt(5),fv(5)]  = updateT(Lm(5),lmo(5),phi(5),fmo(5),vmt(5),mtScal(5),a(5));
%---> Step Four: Calculate the joint moment
moment = fmt.*MA;
%---> Step Five: Calculate the jiont torque
tau = moment(1) + moment(2) + moment(3) + moment(4)+ moment(5);
% tau = (abs(moment(1))+ abs(moment(2))) - (abs(moment(3))+abs(moment(4))+abs(moment(5)));
%---> Step Six: Calculate the acceleration and integrated to theta
% ddtheta = tau/J;
ddtheta = (tau-C*dtheta- K*theta)/inertia;
% Debug purpose
if isnan(ddtheta)
    disp(tau);
    disp(ddtheta);
    disp(J);
end
end