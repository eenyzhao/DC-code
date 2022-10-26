% load subject EMG data and reference trajectory
% ECRL & ECRB share has the same EMG activation
% FCR; FCU; ECRL; ECRB; ECU
clear all
close all
clc

%% Evoluntionary algorithms flags
%% 1 - GA; 2- SA; 3 - PSO ;
evoFlag   = 0;  %flag to turn on optimization
saveFlag  = 1; % save data


%% load subject experiment data
subject         = 'S6';
motion          = 'WFE';
Opimisation     = 'Opti';
fileDir         = '...\ExampleData';
filePath    = strcat(fileDir,'\',subject,'\','_',Opimisation,'.mat');

load(filePath);

time  = redacted_wristData.time;
eFCR  = redacted_wristData.FCR;
eFCU  = redacted_wristData.FCU;
eECRL = redacted_wristData.ECRL;
eECRB = redacted_wristData.ECRB;
eECU  = redacted_wristData.ECU;
refQ  = redacted_wristData.angle;


%% Reorginse the time span
fs = 1000;
dt = 1/fs;
tspan = (length(refQ)-1)/fs;
span = 0:dt:tspan;
span = span';
num = length(span);

% ODE45 solver for muscle activation
tic
[t1,aFCR]  = ode45(@(t,a) odeFunACT(t,a,eFCR,fs),span,eFCR(1));
[t2,aFCU]  = ode45(@(t,a) odeFunACT(t,a,eFCU,fs),span,eFCU(1));
[t3,aECRL] = ode45(@(t,a) odeFunACT(t,a,eECRL,fs),span,eECRL(1));
[t4,aECRB] = ode45(@(t,a) odeFunACT(t,a,eECRB,fs),span,eECRB(1));
[t5,aECU]  = ode45(@(t,a) odeFunACT(t,a,eECU,fs),span,eECU(1));
toc

act = [aFCR,aFCU,aECRL,aECRB,aECU];
%Remove negative value
for col = 1:size(act,2)
    for row = 1:size(act,1)
        if act(row,col) < 0.00 % lower bound for the tendon compliance computation
            act(row,col) = 0.00;
        end
    end
end
% debug

% Initialize the MSK parameters
Lm0 =   [0.062 0.051 0.081  0.058 0.062];   % Lm0
Fm0 =   [407 479 337 252 192];              % Fm0 scale down
Lts =   [0.24  0.26   0.24  0.22  0.2285];  % Lts
Phi =   [0.05 0.2 0.01 0.16 0.06];          % Pennation angle
MtScale  = [1,1,1,1,1];
afactor  = -2;                              % A factor;
% Combine parameter for optimization
initialpara = [Lm0 Fm0 Lts MtScale afactor];

%% Optimization setting

% Lower bound
LB = [0.062*0.85 0.051*0.85 0.081*0.85  0.058*0.85 0.062*0.85...
    407*0.5    479*0.5 337*0.5 252*0.5 192*0.5...
    0.24*0.85  0.26*0.85   0.24*0.85  0.22*0.85  0.2285*0.85...
    0.9 0.9 0.9 0.9 0.9 ...
    -3];
% upper bound
UB =  [0.062*1.15 0.051*1.15 0.081*1.15  0.058*1.15 0.062*1.15...
    407*1.5    479*1.5 337*1.5 252*1.5 192*1.5...
    0.24*1.15  0.26*1.15   0.24*1.15  0.22*1.15  0.2285*1.15...
    1.1 1.1 1.1 1.1 1.1 ......
    0.0001];


%% Optimzation through Different evoluntionary algorithms
%% 1 - GA; 2 - SA; 3 - PSO ;  4 - GSO
if evoFlag == 1 % Genetic algorithm
    nvar = 21;    
    tic
    % no linear constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nlcon = [];
    options = optimoptions('ga','TolFun',1e-2,...
        'PlotFcn',@gaplotbestf,...
        'Display','iter',...
        'MaxGenerations',1000,...
        'PopulationSize',100,...
        'StallGenLimit',51)
    [pOpt,fval,exitflag,output,population,scores] = ga(@(initialpara) objective(initialpara,act,refQ,span),nvar,A,b,Aeq,beq,LB,UB,nlcon,options);
    runTime = toc;
    iterInfo = output.generations;
    disp(['Optimization running time is :', num2str(runTime)]);
    disp(['Optimization is done with in', num2str(iterInfo), ' generations:']);   
elseif evoFlag == 2 % simulated annealing 
    tic
    options = optimoptions('simulannealbnd',...
        'TolFun',1e-2,'PlotFcn',@saplotbestf,...
        'Display','iter',...
        'ReannealInterval',100,...
        'MaxIterations',10000,...
        'StallIterLimit',2000)
    [pOpt,fval,exitflag,output] = simulannealbnd(@(initialpara) objective(initialpara,act,refQ,span),initialpara,LB,UB,options);
    runTime  = toc;
    iterInfo = output.iterations;
    disp(['Optimization running time is :',num2str(runTime)]);
    disp(['Optimization is done with in', num2str(iterInfo), ' interations:']);
elseif evoFlag == 3 % particle swarm optimisation
    nvar = 21;
    tic
    options = optimoptions('particleswarm',...
        'TolFun',1e-2,...
        'MaxIterations',1000,...
        'PlotFcn',@pswplotbestf,...
        'Display','iter',...
        'SwarmSize',100,...
        'MaxStallIterations',20)
    [pOpt,fval,exitflag,output] = particleswarm(@(initialpara) objective(initialpara,act,refQ,span),nvar,LB,UB,options);
    runTime = toc;
    iterInfo = output.iterations;
    disp(['Optimization running time is :',num2str(runTime)]);
    disp(['Optimization is done with in', num2str(iterInfo), ' interations:']);
end

% save the optimised parameters
if evoFlag == 1 == 1 && saveFlag == 1
    char = [subject,'_','pOpt_ga'];
    save(char,'pOpt','runTime','iterInfo');
elseif evoFlag == 2 && saveFlag == 1
    char = [subject,'_','pOpt_sa'];
    save(char,'pOpt','runTime','iterInfo');
elseif evoFlag == 3 && saveFlag ==1
    char = [subject,'_','pOpt_pso'];
    save(char,'pOpt','runTime','iterInfo');
end

[preQ,a] = simulation(pOpt,act,refQ,span);

h(1) = figure(1);
plot(span,preQ,'r:','LineWidth',1.2)
hold on
plot(span,refQ,'k-','LineWidth',1.2)
legend('Estimated angle','Ground Truth')
ylabel('joint angle(deg)')
xlabel('t(s)')
%% Plot and Save Results
% Display optimized muscule-tendon parameters
lmo_opt = pOpt(1:5);                    % Optimal muscle fibre length
fmo_opt = pOpt(6:10);                   % maximum isometric force
lts_opt = pOpt(11:15);                  % tendon slack length
% phi = [0.05 0.2 0.01 0.16 0.06];      % phi
mtscaler_opt = pOpt(16:20);             % scale for Lmt and MA
afactor_opt  = pOpt(21);                % muscle activation non linear shape factor
disp('The optimized EMG-driven model parameters are:');
disp(['Lm0 are:', num2str(lmo_opt)]);
disp(['Fm0 are:', num2str(fmo_opt)]);
disp(['Lts are:', num2str(lts_opt)]);
disp(['Phi are:', num2str(Phi)]);
disp(['MTLscaler are:', num2str(mtscaler_opt)]);
disp(['afactor is:', num2str(afactor_opt)]);

condition = @(x) (x - mean(x))./std(x);
CC = condition(preQ)'*condition(refQ)/sum(condition(preQ).^2);
[r2,rmse] = rsquare(refQ,preQ);
disp(['CC is:',num2str(CC)]);
disp(['R^2 is:',num2str(r2)]);
disp(['RMSE is:',num2str(rmse)]);
