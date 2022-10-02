function G = objGrad(x,aux)
%% gradient vector is computed with respect to objective function
%% --> Get measured theta from aux structure
gridN = aux.gridN;
meas_theta = aux.reftheta;
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
% some  variables
Ncon  = (gridN -1)*(Nstate+Npara)+7; % The number of constraints
Nx    = (gridN*(Nstate+Nctrl));       % # of variables


% Size of gradient matrix
G = zeros(1,Nx);
%% --> Get the state and control out of the vector(Column = state/control; row = grid point)
for i = 1:gridN
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
        
%     G(indtheta) = -2*1000*h*(meas_theta(i)-theta(i));   % gradient with respect to theta
    G(indtheta) = (-2*1000*(meas_theta(i)-theta(i)))/length(meas_theta);   % gradient with respect to theta
%     G(indemg)   = -1000*h*2*(e(i,:)-emg(i,:));        % gradient with respect to activation
%     G(indact)   = 1000*h*2*(act(i,:) - emg(i,:));  % gradient with respect to activation
%     G(indemg)   = -1000*h*2*(act(i,:) - emg(i,:)); % gradient with respect to activation
end
end