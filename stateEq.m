function  F = stateEq(theta,dtheta,ddtheta,a,parameter)
%% implicit state equation for direction collocation method
%% Declare hand parameters for each subjects;
bodyWeight = 88;             % Kg -- subject's body mass kg
bodyHeight = 178;            % cm -- subject's height cm
m = bodyWeight*0.0065;       % --> hand mass (Kg)
g = 9.81;                    % --> gravity coefficient (m/s^2)

l = 0.08;                    % -- hand length (measured)
J = (-13.68 + 0.088*bodyWeight + 0.092*bodyHeight)/10000;% --> inertia of hand
inertia = J + m*l^2;         % --> inertia around joint rotation center

%% allocate the parameter
lmo = parameter(1,1:5);             % optimal muscle fibre length
fmo = parameter(1,6:10);            % maximum isometric force
lts = parameter(1,11:15);           % tendon slack length
mtScal = parameter(1,16:20);        % tendon slack length
phi = [0.05 0.2 0.01 0.16 0.06];    % Pennation angle
B = 0.3;                            % damping parameters
K = 0.6;                            % stiffness parameters
%% Compute the joint torque using musculoskeletal model
%---> Step one: with state of Q, calculate the moment arm and muscle-tendon length
[Lmt,MA] = getMTUandMA(theta);
Lmt = Lmt.*mtScal;
MA  = MA.*mtScal;
% --->Step two: Calculate the current muscle fibre length
% Lm2  = updateLm(lmo,lts,phi,Lmt);
Lm = sqrt((lmo.*sin(phi)).^2 + (Lmt - lts).^2);
% --->Step Three: Calculate the joint moment
fmt(1) = updateT(Lm(1),lmo(1),phi(1),fmo(1),mtScal(1),a(1));
fmt(2) = updateT(Lm(2),lmo(2),phi(2),fmo(2),mtScal(2),a(2));
fmt(3) = updateT(Lm(3),lmo(3),phi(3),fmo(3),mtScal(3),a(3));
fmt(4) = updateT(Lm(4),lmo(4),phi(4),fmo(4),mtScal(4),a(4));
fmt(5) = updateT(Lm(5),lmo(5),phi(5),fmo(5),mtScal(5),a(5));
%---> Step Four: Calculate the joint moment
moment = fmt.*MA;
%---> Step Five: Calculate the jiont torque
tau = moment(1) + moment(2) + (moment(3) + moment(4)+ moment(5));
%---> Step Six: Calculate the acceleration and integrated to theta;
%% ddtheta*J + dampingC*dtheta+ m*g*l*sin(theta) - tau = 0;
F = ddtheta*inertia + B*dtheta + K*theta - tau;
%
end