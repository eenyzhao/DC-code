%% Function to get the theta
function [thetaP,a] = simulation(p,act,refTheta,time)

afactor = p(21);
threshold = 0; % threshold is set to remove the artefact
a = activation(afactor,act,threshold);

%pre_allocate space fore
dtheta = zeros(length(time),1);
ddtheta = zeros(length(time),1);
theta = zeros(length(time),1);
vmt = zeros(length(time),5);
Lmt = zeros(length(time),5);
Lm  = zeros(length(time),5);
MA  = zeros(length(time),5);
% Initial condition
theta(1) = refTheta(1);
dtheta(1) = 0;
ddtheta(1) = 0;

%% Runge kutta 4th order 4th order
for i = 1:length(time)-1
    dt = time(i+1) - time(i);
    [Lmt(i,:),Lm(i,:),ddtheta(i),fv(i,:)] = thetaPred(theta(i),dtheta(i),a(i,:),p,vmt(i,:));
    
    % kv1 & ks1
    kv1 = ddtheta(i)*dt;
    ks1 = dtheta(i)*dt;
    
    % kv2 & ks2
    kv2 = (ddtheta(i)+ ks1*0.5)*dt;
    ks2 = (dtheta(i) + 0.5*kv1)*dt;
    
    % kv3 & ks3
    kv3 = (ddtheta(i)+ ks2*0.5)*dt;
    ks3 = (dtheta(i) + 0.5*kv2)*dt;
    
    % kv4 & ks4
    kv4 = (ddtheta(i)+ ks3*0.5)*dt;
    ks4 = (dtheta(i) + kv3)*dt;
    
    % Wrist flexion/extension integration
    dtheta(i+1) = dtheta(i) + 1/6*(kv1+2*kv2+2*kv3+kv4);
    theta(i+1)  = theta(i)  + 1/6*(ks1+2*ks2+2*ks3+ks4);
    
    [Lmt(i+1,:),MA(i+1,:)] = getMTUandMA(theta(i+1));
    vmt(i+1,:) = (Lmt(i+1,:) - Lmt(i,:))/dt;    
end
thetaP = theta;
end