function obj = objective(x,aux)
%% Define the objective function
% --> Get measured theta from aux structure
gridN = aux.gridN;
meas_theta = aux.reftheta';
h = aux.h;
Nstate = aux.Nstate;
Nctrl = aux.Nctrl;
Npara = aux.Npara;
uFCR = aux.uFCR;
uFCU = aux.uFCU;
uECRL = aux.uECRL;
uECRB = aux.uECRB;
uECU = aux.uECU;

e = [uFCR,uFCU,uECRL,uECRB,uECU];
% --> Get the state and control out of the vector(Column = state/control; row = grid point)
for i = 1:gridN
    % Index information
    % Index information       
    indtheta      = (i-1)*(Nstate+Nctrl) + 1;                              % theta
    indtheta_dot  = (i-1)*(Nstate+Nctrl) + 2;                              % dtheta
    indact        = (i-1)*(Nstate+Nctrl) + 3 :  (i-1)*(Nstate+Nctrl) + 7;   % act
    indemg        = (i-1)*(Nstate+Nctrl) + 8 :  (i-1)*(Nstate+Nctrl) + 12;  % u
    indpara       = (i-1)*(Nstate+Nctrl) + 13 : (i-1)*(Nstate+Nctrl) + 33; % msk parameters
    
    theta(i)  = x(indtheta);
    dtheta(i) = x(indtheta_dot);
    act(i,:)  = x(indact);
    emg(i,:)  = x(indemg);
    para(i,:) = x(indpara);
end

% --> Calculate the objective value
% --> Minimize the sum of square error between esitmated theta and measured
% --> theta
Err_emg = e - emg;
EMG_sumErr = 0;
for i = 1:size(emg,2)  % penalty function
    EMG_sumErr = EMG_sumErr + (h * sum(Err_emg(:,i).^2));
end

% Err_theta = 1000* h * sum((meas_theta - theta).^2) ;
weight = 1000; % add weight if sometimes it convergence too fast
Err_theta = sum((meas_theta - theta).^2)/length(meas_theta) ;

% obj = h * sum(Err_theta.^2);
obj = Err_theta;

% debug
% disp(['The objective value is ' num2str(obj)]);

end
