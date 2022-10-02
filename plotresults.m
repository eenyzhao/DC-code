function [] = plotresults(x,aux)

%% --> Get measured theta from aux structure
gridN = aux.gridN;
meas_theta = aux.reftheta;
h = aux.h;
Nstate = aux.Nstate;
Nctrl = aux.Nctrl;
duration = aux.duration;
dc_time = aux.time;
Npara  = aux.Npara;

theta = zeros(gridN,1);
%% Get the state and control out of the vector
for i = 1:gridN
    % index information
    indtheta      = (i-1)*(Nstate+Nctrl) + 1;                              % theta
    indtheta_dot  = (i-1)*(Nstate+Nctrl) + 2;                              % dtheta
    indact        = (i-1)*(Nstate+Nctrl) + 3 : (i-1)*(Nstate+Nctrl) + 7;   % act
    indemg        = (i-1)*(Nstate+Nctrl) + 8 : (i-1)*(Nstate+Nctrl) + 12;  % u
    indpara       = (i-1)*(Nstate+Nctrl) + 13 : (i-1)*(Nstate+Nctrl) + 33; % msk parameter
    theta(i) = x(indtheta);
    dtheta(i) = x(indtheta_dot);
    act (i,:) = x(indact);
    emg(i,:) = x(indemg);
    para(i,:) = x(indpara);
end

%% The optimized parameters value
lmo_opt = para(1,1:5);                  % optimal muscle fibre length
fmo_opt = para(1,6:10);                 % maximum isometric force
lts_opt = para(1,11:15);               % tendon slack length
phi_opt = [0.05 0.2 0.01 0.16 0.06];    % pennation angle
mtScal = para(1,16:20);             % tendon slack length
% non_afactor = para(21);
afactor_opt  = para(1,21);  % muscle activation non linear shape factor
disp('The EMG-driven model parameters are:');
disp(['Lm0 are:', num2str(lmo_opt)]);
disp(['Fm0 are:', num2str(fmo_opt)]);
disp(['Lts are:', num2str(lts_opt)]);
disp(['Phi are:', num2str(phi_opt)]);
disp(['MtScale are:', num2str(mtScal)]);
disp(['afactor is :',num2str(afactor_opt)]);

%% The optimized RMSE and Correlation
% --> Caculate RMSE
error = meas_theta - theta;
rmse = sqrt(sum(error.^2)/length(theta));
disp(['The RMSE is' num2str(rmse)])

% --> Calculate Correlation
condition = @(z) (z - mean(z))./std(z);
RR = condition(theta)'*condition(meas_theta)/sum(condition(theta).^2);
% R = corrcoef(theta,refTheta);
disp(['The CC  is ' num2str(RR)])
%% show muscle activation
figure()
subplot(2,1,1)
plot(dc_time,act,'LineWidth',1)
legend('muscle activation')
subplot(2,1,2)
plot(dc_time,emg,'LineWidth',1)
legend('Neural excitation')

%% Show optimal result
figure()
% subplot(2,1,1)
plot(dc_time,theta,'r-.','LineWidth',1)
hold on
plot(dc_time,meas_theta,'k:','LineWidth',1)
% set(gca,'XTick',[0:5:20]) % Change the intervel
% set(gca,'YTick',[-1.5:0.5:1.5]) % Change the intervel
legend('Estimated Joint Angle','Measured Joint Angle','location','NorthOutside');
xlabel('time(s)','FontWeight','bold')
ylabel('Joint Angle(rad)','FontWeight','bold')
hold off
end