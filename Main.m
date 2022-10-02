%% Matlab script to optimize the parameter using direct collocation methods
clear all; close all; clc;

%% load subject EMG data and reference trajectory
% ECRL & ECRB share has the same EMG activation
% FCR; FCU; ECRL; ECRB; ECU
saveFlag = 0;
%% Optimzated parameter for prediction 
%% load subject experiment data

load('eg1.mat');

% Data reollocation
time  = redacted_wristData.time;
eFCR  = redacted_wristData.FCR;
eFCU  = redacted_wristData.FCU;
eECRL = redacted_wristData.ECRL;
eECRB = redacted_wristData.ECRB;
eECU  = redacted_wristData.ECU;
refQ  = redacted_wristData.angle;

sampleF = 1000; % sample frequency
% re-arrange time start from zero
tspan = time(end) - time(1);
dt = tspan/(length(time)-1);
t = 0:dt:tspan;
t = t';

e = [eFCR,eFCU,eECRL,eECRB,eECU];

%% Direct collocation methods setup
gridN = 201;                % number of grid point
duration = t(end)-t(1);     % time in sec
h = duration/(gridN-1);     % time inverval between nodes
dc_time = t(1):h:t(end);    % list of time points
dc_time = dc_time';

%% allocate the reference theta into grid point
meas_theta = interp1(t,refQ,dc_time);
uFCR       = interp1(t,eFCR,dc_time);
uFCU       = interp1(t,eFCU,dc_time);
uECRL      = interp1(t,eECRL,dc_time);
uECRB      = interp1(t,eECRB,dc_time);
uECU       = interp1(t,eECU,dc_time);
ss         = [uFCR,uFCU,uECRL,uECRB,uECU];

% use a struct to store the auxiliary information to pass in optimizer
aux.gridN = gridN;
aux.duration = duration;
aux.h = h;
aux.time = dc_time;
aux.reftheta = meas_theta;
aux.uFCR = uFCR;
aux.uFCU = uFCU;
aux.uECRL = uECRL;
aux.uECRB = uECRB;
aux.uECU = uECU;
%% static parameter initialization for ipopt optimization
% 1st initial guess of parameters
Lm0 =   [0.062 0.051 0.081  0.058 0.062];   % Lm0
Fm0 =   [407 479 337 252 192];              % Fm0 scale down
Lts =   [0.24  0.26   0.24  0.22  0.2285];  % Lts
Phi =   [0.05 0.2 0.01 0.16 0.06];          % Pennation angle (not Optimised)
MtScale = [1,1,1,1,1];
afactor  = -2;                              % Non-linear  factor;
% Combine parameter for optimization 
initialpara = [Lm0 Fm0 Lts MtScale afactor];


Npara = length(initialpara);
%% allocate the control variable into grid point
% add the static paramters into control variables
Nctrl = size(e,2)+ Npara;           % number of control variable with added static parameters
ctrl_temp = zeros(gridN,Nctrl);     % temporal control variable
% ctrl_u = zeros(size(e,1),size(e,2));

for j = 1:size(e,2)
    ctrl_temp(:,j)= interp1(t,e(:,j),dc_time);
end
for row = 1:gridN
    ctrl_temp(row,size(e,2)+1:end) = initialpara;
end

%% initial guess of state variable = [angle; angle velocity; activation]
theta  = (70*pi/180 - (-70*pi/180))*rand(length(t),1) +(-70*pi/180) ;       % joint angle 
dtheta = randn(length(t),1);         % joint velocity
a_state = 1*ones(length(t),5);  % muscle activation a factor
% a_state = e;
state_var = [theta dtheta a_state];

Nstate = 7; % number of state variable
state_temp  = zeros(gridN,Nstate); % temporal state variable
for k = 1 : Nstate
    state_temp(:,k) = interp1(t,state_var(:,k),dc_time);
end
%% For direct collocation method, create the initial parameter vector store the state, control, static parameters
%% in the form of X0 = [x1,u1,p,x2,u2,p,x3,u3,p ...xn,un,p];
% store the state and control number into aux
aux.Nstate = Nstate;
aux.Nctrl  = Nctrl;
aux.Npara  = Npara;

Ncon = (gridN -1)*(Nstate+Npara)+7; % # of constraints
Nx   = (gridN*(Nstate+Nctrl)); % # of variables;

x0 = zeros(1,Nx);
% Create the initial guess
for i = 1:gridN
    ix_first = (i-1)*(Nstate+Nctrl) + 1; 
    ix_last  = ix_first + (Nstate+Nctrl) - 1;
    ix = ix_first:ix_last;
    x = state_temp(i,:);
    u = ctrl_temp(i,:);
    x0(ix) = horzcat(x,u);
end

%% Boundary condition for state/control/parameters

[LB,UB] = Varbounds(aux,ctrl_temp);
% ---> ipopt setting
% solve the NLP with IPOPT
% Ipopt options
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.recalc_y_feas_tol = 1e-02;

% options.ipopt.max_cpu_time    = 43200;
options.ipopt.max_iter        = 1000;
options.ipopt.tol             = 1e-02;
options.ipopt.dual_inf_tol    = 1e-02;
options.ipopt.constr_viol_tol = 1e-03;
options.ipopt.compl_inf_tol   = 1e-03;
options.ipopt.acceptable_tol  = 1e-03;
options.ipopt.acceptable_iter = 15;

options.ipopt.bound_frac = 0.0001;  %0.0001;
options.ipopt.bound_push = 0.0001;  %0.0001;
options.ipopt.ma57_pre_alloc        = 2.0;


% Boundary Conditions
options.lb = LB;
options.ub = UB;
options.cl = zeros(1,Ncon);
options.cu = zeros(1,Ncon);
options.auxdata = aux;

% IPOPT functions
funcs.objective   = @objective;
funcs.gradient    = @objGrad;
funcs.constraints = @constraints;
funcs.jacobian    = @conJacobian;
funcs.jacobianstructure = @conJacobianstructure;

%% Debug purpose

tic; % start timer
[xopt, info] = ipopt_auxdata(x0,funcs,options);
iterInfo = info.iter;
executiontime = toc; % stop timer

% Load xopt
% load('D:\WFH\NMS_IPOPT_Ver2\xopt_ipopt\S8_xopt_dc.mat')
% Display time 
fprintf('Finished in %f seconds\n', executiontime);

% Plot optimized results
pOpt= xopt((Nstate+Nctrl) + 13 : (Nstate+Nctrl) + 33);  % static parameter
plotresults(xopt,aux);
figure
plot(dc_time,ss)
lmo_opt = pOpt(1:5);                    % Optimal muscle fibre length
fmo_opt = pOpt(6:10);                   % maximum isometric force
lts_opt = pOpt(11:15);                  % tendon slack length
mtscaler_opt = pOpt(16:20);             % scale for Lmt and MA
afactor_opt  = pOpt(21);                % muscle activation non linear shape factor


if saveFlag == 1 
    muscleInd  = {'FCR';'FCU';'ECRL';'ECRB';'ECU';};
    lmo_ind    =  lmo_opt';
    fmo_ind    =  fmo_opt';
    lts_ind    =  lts_opt';
    phi_ind    =  Phi';
    mtScal_ind  = mtscaler_opt';
    afactor_ind = [0;0;0;0;afactor_opt];
    varNames    = {'muscleInd','Lmo','Fmo','Lts','Phi','MtScal','afator'};
    outTable = table(muscleInd,lmo_ind,fmo_ind,lts_ind,...
                     phi_ind,mtScal_ind,afactor_ind,...
                     'VariableNames',varNames);
    saveFileName = strcat(subject,'_dc','_optiPara.xlsx');
    writetable(outTable,saveFileName);
    char = [subject,'_','pOpt_dc'];
    save(char,'pOpt','executiontime','iterInfo');
end

%% debug the MSK parameter maintain constant 
for num = 1:gridN
    indpara = (num-1)*(Nstate+Nctrl) + 13 : (num-1)*(Nstate+Nctrl) + 33;
    para_opt(num,:) = xopt(indpara);
end

% Constraints violation
c = constraints(xopt,aux);
figure
plot(c); 
title('Constraints violation');


