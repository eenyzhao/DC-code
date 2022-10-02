function [LB,UB] = Varbounds(aux,ctrl_temp)
gridN = aux.gridN;
meas_theta = aux.reftheta;
h = aux.h;
Nstate = aux.Nstate;
Nctrl = aux.Nctrl;
Ncon = (gridN -1)*7+10; % The number of constraints
Nx   = (gridN*(Nstate+Nctrl)); % # of total variables
%% Boundary condition for state/control/parameters
% Boundary conditions for MSK parameters at each node
% lower bound --> muscle tendon parameter
LB_para = [0.062*0.85 0.051*0.85 0.081*0.85  0.058*0.85 0.062*0.85...
           407*0.5    479*0.5 337*0.5 252*0.5 192*0.5...
           0.24*0.85  0.26*0.85   0.24*0.85  0.22*0.85  0.2285*0.85...
           0.9 0.9 0.9 0.9 0.9 ...
           -3]; 
% upper bound --> muscle tendon parameter
UB_para = [0.062*1.15 0.051*1.15 0.081*1.15  0.058*1.15 0.062*1.15... 
           407*1.5    479*1.5 337*1.5 252*1.5 192*1.5...
           0.24*1.15  0.26*1.15   0.24*1.15  0.22*1.15  0.2285*1.15...
           1.1 1.1 1.1 1.1 1.1 ......
           0.0001]; 

%% re-arrange the boundary condition in the form of
%% [x1,u1,x2,u2,x3,u3,....xn,un];
LB = zeros(1,Nx);
UB = zeros(1,Nx);
ctrl_emg = zeros(gridN,5);

% reallocation the EMG and MSK parameters
for row = 1:gridN
    ctrl_emg(row,:) = ctrl_temp(row,1:5);
end



for k = 1:gridN
    ind_first = (k-1)*(Nstate+Nctrl) + 1; 
    ind_last  = ind_first + (Nstate+Nctrl) - 1;
    ind_bound = ind_first:ind_last;
    
    % boundary for states
    LB_theta  = -(70*pi/180);               % Kinematic
    LB_dtheta = -inf;                       % Kinematic
    LB_act    = 0.0001*ones(1,size(ctrl_emg,2));       % Activation state
    
    UB_theta  = (70*pi/180);                % Kinematic
    UB_dtheta = inf;                        % Kinematic
    UB_act    = 0.999*ones(1,size(ctrl_emg,2));        % Activation state
    
    % boundary for the control   
    LB_ctrl   = ctrl_emg(k,:);   % Control emg
    UB_ctrl   = ctrl_emg(k,:);    % Control EMG with be no change  
    
    LB(ind_bound) = horzcat(LB_theta,LB_dtheta,LB_act,LB_ctrl,LB_para);
    UB(ind_bound) = horzcat(UB_theta,UB_dtheta,UB_act,UB_ctrl,UB_para);
end

% LB(gridN * (Nstate+Nctrl) + 1:end) = LB_parameter;
% UB(gridN * (Nstate+Nctrl) + 1:end) = UB_parameter;

end
