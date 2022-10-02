function C = constraints(x,aux)
% This script is used to impose the Hill model constraints
% And the contrainst setting for EMG-driven model
%% --> Get measured theta from aux structure
gridN      = aux.gridN;
meas_theta = aux.reftheta;
h          = aux.h;
Nstate     = aux.Nstate;
Nctrl      = aux.Nctrl;
Npara      = aux.Npara;

Ncon   = (gridN -1)*(Nstate+Npara)+7; % The number of constraints
                                    % Without task contraints
tact   = 0.015;  % activation time constant   - 15ms
tdeact = 0.050;  % deactivation time constant - 50ms
%% --> Get the state and control out of the vector(Column = state/control; row = grid point)
theta  = zeros(gridN,1);
dtheta = zeros(gridN,1);
act    = zeros(gridN,5);
emg    = zeros(gridN,5);
para   = zeros(gridN,21);

for i = 1:gridN
    % index information
    indtheta = (i-1)*(Nstate+Nctrl) + 1;
    indtheta_dot = (i-1)*(Nstate+Nctrl) + 2;
    indact    = (i-1)*(Nstate+Nctrl) + 3 : (i-1)*(Nstate+Nctrl) + 7;
    indemg    = (i-1)*(Nstate+Nctrl) + 8 : (i-1)*(Nstate+Nctrl) + 12;
    indpara   = (i-1)*(Nstate+Nctrl) + 13: (i-1)*(Nstate+Nctrl) + 33;
    theta(i)  = x(indtheta);     % angles
    dtheta(i) = x(indtheta_dot); % velocity
    act(i,:) = x(indact);       % activation
    emg(i,:)  = x(indemg);       % emg
    para(i,:) = x(indpara);      % parameters    
end

%% By discretizing the control and state variable
%  At each grid point, the defect calculated by euler method and
%  the differenrtiation of the state with respect time should be zero.
%  Second methods: Using the implicit formulation to compute the
%  constraints;
%  In this script, the mid-point rule is used;
%% Equality constraints
% Pre-alloacte space for algebric equality constraints
C = zeros(1,Ncon);
c_Temp = zeros(1,Nstate+Npara);
act_dot = zeros(1,size(act,2));
a_non = zeros(1,size(act,2));
u = zeros(1,size(emg,2));
a = zeros(1,size(act,2));

for k = 1:gridN-1
    % Constraints index
    indC = ((k-1)*(Nstate+Npara) + 1) : ((Nstate+Npara)*k);
    % state derivative
    theta_dot  = (theta(k+1)  - theta(k))/h;   % velocity 
    dtheta_dot = (dtheta(k+1) - dtheta(k))/h;  % acceleration
    act_dot    = (act(k+1,:)  - act(k,:))/h;   
    
    % state and control at current grid
    q = (theta(k)    + theta(k+1))  /2;
    v = (dtheta(k)   + dtheta(k+1)) /2;
    u = (emg(k,:)    + emg(k+1,:))  /2;
    a = (act(k+1,:)  + act(k,:))    /2;
    p = (para(k+1,:) + para(k,:))  /2;
    % Compute the path constraints using implicited formulation of
    % musculoskeletal model
    % Equation of motion
    c_Temp(1) = theta_dot - v;    
    
    % MSK dynamics constraints
    a_non= (( exp( p(21)* a) - 1)./(exp(p(21))- 1 ) ); % non-linear
              
    c_Temp(2) = stateEq(q,v,dtheta_dot,a_non,p);   
    % muscle activation constraints
    for kk = 1:size(emg,2)
        c_Temp(kk+2) = act_dot(kk) - (u(kk)/tact + ((1-u(kk))/tdeact))*(u(kk)-a(kk));
    end
    
    for num = 1:size(para,2)
        c_Temp(num+Nstate) = para(k+1,num) - para(k,num);
    end    
    C(indC) = c_Temp;
end

%% task contraints
% At the initial condition
% --> state at initial condition
C((gridN -1)*(Nstate+Npara)+1) = theta(1) - meas_theta(1);        % start position
C((gridN -1)*(Nstate+Npara)+2) = dtheta(1);                       % initial joint velocity == 0;

% initial activation state
C((gridN -1)*(Nstate+Npara)+3)  = act(1,1)- emg(1,1);
C((gridN -1)*(Nstate+Npara)+4)  = act(1,2)- emg(1,2);
C((gridN -1)*(Nstate+Npara)+5)  = act(1,3)- emg(1,3);
C((gridN -1)*(Nstate+Npara)+6)  = act(1,4)- emg(1,4);
C((gridN -1)*(Nstate+Npara)+7)  = act(1,5)- emg(1,5);
end