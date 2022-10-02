%% Hessian matrix H, gradient vector G and Jacobian matrix J
function J = conJacobian(x,aux)
% The jacobian matrix of the constraints
gridN = aux.gridN;
meas_theta = aux.reftheta;
h = aux.h;
Nstate = aux.Nstate;
Nctrl = aux.Nctrl;
Npara = aux.Npara;

% some constant variables
Ncon   = (gridN -1)*(Nstate+Npara)+7; % The number of constraints
Nx     = (gridN*(Nstate+Nctrl)); % # of variables the +20 is the muscle-tendon parameters

tact   = 0.015;  % activation time constant   - 15ms
tdeact = 0.050;  % deactivation time constant - 50ms

% Some constant in musculoskeletal model
%% Declare hand parameters for each subjects;
bodyWeight = 88;             % Kg -- subject's body mass kg
bodyHeight = 178;            % cm -- subject's height cm
m = bodyWeight*0.0065;       % --> hand mass (Kg)
g = 9.81;                    % --> gravity coefficient (m/s^2)

l = 0.08;                    % -- hand length (measured)
J = (-13.68 + 0.088*bodyWeight + 0.092*bodyHeight)/10000;% --> inertia of hand
inertia = J + m*l^2;         % --> inertia around joint rotation center
B = 0.3;                            % damping parameters
K = 0.6;                            % stiffness parameters

% --> Get the state and control out of the vector(Column = state/control; row = grid point)
% Size of Jacobian
J = zeros(Ncon,Nx);
%p = x(gridN * (Nstate+Nctrl) + 1:end); % static parameter
%% allocate the parameter
% Lmo = p(1,1:5);       % optimal muscle fibre length
% fmo = p(1,6:10);      % maximum isometric force
% lts = p(1,11:15);     % tendon slack length
% phi =   [0.05 0.2 0.01 0.16 0.06];
% pre-allocate space for lmt and ma
% Lmt  = zeros(1,5);
% MA   = zeros(1,5);
% Lm   = zeros(1,5);

for n = 1:gridN-1
    % index info for state, control and parameters
    % Index information
    % current node N:
    indtheta      = (n-1)*(Nstate+Nctrl) + 1;                              % theta
    indtheta_dot  = (n-1)*(Nstate+Nctrl) + 2;                              % dtheta
    indact        = (n-1)*(Nstate+Nctrl) + 3 : (n-1)*(Nstate+Nctrl) + 7;   % act
    indemg        = (n-1)*(Nstate+Nctrl) + 8 : (n-1)*(Nstate+Nctrl) + 12;  % u
    indpara       = (n-1)*(Nstate+Nctrl) + 13 : (n-1)*(Nstate+Nctrl) + 33; % msk parameters
    
    % node N+1:
    indtheta2     = (n)*(Nstate+Nctrl) + 1;                             % theta
    indtheta_dot2 = (n)*(Nstate+Nctrl) + 2;                             % dtheta
    indact2       = (n)*(Nstate+Nctrl) + 3 : (n)*(Nstate+Nctrl) + 7;    % act
    indemg2       = (n)*(Nstate+Nctrl) + 8 : (n)*(Nstate+Nctrl) + 12;   % u
    indpara2      = (n)*(Nstate+Nctrl) + 13 :(n)*(Nstate+Nctrl) + 33;  % msk parameters
    
    % midpoint rule
    % state at current grid
    q1 =  x(indtheta);       % theta at current N
    q2 =  x(indtheta2);      % theta at N+1;
    q  = (q1+q2)/2;          % theta_dot at current N
    v1 = x(indtheta_dot);    % theta_dot at current N
    v2 = x(indtheta_dot2);   % theta_dot at N+1
    v  = (v1+v2)/2;
    a1 = x(indact);          % act at current N
    a2 = x(indact2);         % act at current N+1
    a  = (a1 + a2)/2;        % state at curret node
    u1 = x(indemg);          % control at current node N
    u2 = x(indemg2);         % control at next node N+1
    u  = (u1+u2)/2;          % control at current node
    
    para1 = x(indpara);     % MT parameters at current N
    para2 = x(indpara2);    % MT parameters at current N+1
    para  = (para1+para2)/2;
    
%     % re-allocate parameters;
%     Lmo1 = para1(1,1:5);                  % optimal muscle fibre length
%     fmo1 = para1(1,6:10);                 % maximum isometric force
%     lts1 = para1(1,11:15);                % tendon slack length  
%     mtScal1 = para1(1,16:20);             % tendon slack length
%     
%     % re-allocate parameters;
%     Lmo2 = para2(1,1:5);                  % optimal muscle fibre length
%     fmo2 = para2(1,6:10);                 % maximum isometric force
%     lts2 = para2(1,11:15);                % tendon slack length
%     mtScal1 = para2(1,16:20);             % tendon slack length
    
    
    % re-allocate parameters;
    Lmo = para(1,1:5);                  % optimal muscle fibre length
    fmo = para(1,6:10);                 % maximum isometric force
    lts = para(1,11:15);                % tendon slack length
    phi = [0.05 0.2 0.01 0.16 0.06];    % pennation angle
    mtScal = para(1,16:20);             % tendon slack length
    % non_afactor = para(21);
    
    % some calculation --- don not needed
%     [Lmt2,MA2] = getMTUandMA(q2); % Lmt at n+1
%     [Lmt1,MA1] = getMTUandMA(q1); % Lmt at n
%     vmt      = (Lmt2  - Lmt1)/h;   % muscle velocity
    
    % Lmt = Lmt.*mtScal;
    % MA  = MA.*mtScal;
    % % --->Step two: Calculate the current muscle fibre length
    % Lm = sqrt((Lmo.*sin(phi)).^2 + (Lmt - lts).^2);
    
    % find the Non-zero elements;   
    %% Equation of motion dfdx
    J((Nstate+Npara)*(n-1)+1,indtheta)        = -1/h;
    J((Nstate+Npara)*(n-1)+1,indtheta2)       =  1/h;
    J((Nstate+Npara)*(n-1)+1,indtheta_dot)    =  -1/2;
    J((Nstate+Npara)*(n-1)+1,indtheta_dot2)   =  -1/2;
    
    %% Musculoskeletal model constraints;
    % I*ddtheta + B*detheta + K*theta - tau = 0
    %
    J((Nstate+Npara)*(n-1)+2,indtheta)                = 1/2*K - dtaudq(q1,q2,Lmo,phi,fmo,lts,mtScal,a,para(21));     % theta n
    J((Nstate+Npara)*(n-1)+2,indtheta2)               = 1/2*K - dtaudq(q1,q2,Lmo,phi,fmo,lts,mtScal,a,para(21));     % thata n+1
    
    J((Nstate+Npara)*(n-1)+2,indtheta_dot)            = inertia*(-1/h) + B*1/2;     % theta_dot;
    J((Nstate+Npara)*(n-1)+2,indtheta_dot2)           = inertia*(1/h)  + B*1/2;     % theta_dot at n +1;
    
    % w.s.t muscle activation
    dcda = dtauda(q,Lmo,lts,phi,fmo,mtScal,para(21),a);
   
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) +3) = dcda(1);            % a1
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl)   +3) = dcda(1);            % a1 at n+1
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) +4) = dcda(2); % a2
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl)   +4) = dcda(2); % a2 at n+1
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) +5) = dcda(3); % a3
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl)   +5) = dcda(3); % a3 at n+1
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) +6) = dcda(4); % a4
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl)   +6) = dcda(4); % a4 at n+1
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) +7) = dcda(5); % a5
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl)   +7) = dcda(5); % a5 at n+1
       
    % Partical derivative w.r.t lmo
    dcdlm = lm0_pd(q,Lmo,lts,phi,fmo,mtScal,a,para(21));    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 13) = dcdlm(1);     %lm0_FCR;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 13) = dcdlm(1);       %lm0_FCR at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 14) = dcdlm(2);     %lm0_FCU;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 14) = dcdlm(2);       %lm0_FCU at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 15) = dcdlm(3);     %lm0_ECRL;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 15) = dcdlm(3);       %lm0_ECRL at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 16) = dcdlm(4);     %lm0_ECRB;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 16) = dcdlm(4);       %lm0_ECRB at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 17) = dcdlm(5);     %lm0_ECU;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 17) = dcdlm(5);       %lm0_ECU at N+1;
    
    % w.r.t fmo
    dcdfmo = fm0_pd(q,Lmo,phi,lts,mtScal,a,para(21));
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 18) = dcdfmo(1);     %fm0_FCR;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 18) = dcdfmo(1);       %fm0_FCR at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 19) = dcdfmo(2);     %fm0_FCU;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 19) = dcdfmo(2);       %fm0_FCUat N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 20) = dcdfmo(3);     %fm0_ECRL;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 20) = dcdfmo(3);       %fm0_ECRL at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 21) = dcdfmo(4);     %fm0_ECRB;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 21) = dcdfmo(4);       %fm0_ECRB at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 22) = dcdfmo(5);     %fm0_ECU;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 22) = dcdfmo(5);       %fm0_ECU at N+1;
    
    % w.r.t lts
    dcdlts = lts_pd(q,Lmo,lts,fmo,phi,mtScal,a,para(21));
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 23) = dcdlts(1);     %lts_FCR;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 23) = dcdlts(1);       %lts_FCR at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 24) = dcdlts(2);     %lts_FCU;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 24) = dcdlts(2);       %lts_FCUat N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 25) = dcdlts(3);     %lts_ECRL;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 25) = dcdlts(3);       %lts_ECRL at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 26) = dcdlts(4);     %lts_ECRB;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 26) = dcdlts(4);       %lts_ECRB at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 27) = dcdlts(5);     %lts_ECU;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 27) = dcdlts(5);       %lts_ECU at N+1;
    
   
    % w.r.t mtScal
    dcdmtScal = mtScal_pd(q,Lmo,lts,fmo,phi,a,mtScal,para(21));

    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 28) = dcdmtScal(1);     %mtScal_FCR;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 28) = dcdmtScal(1);       %mtScal_FCR at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 29) = dcdmtScal(2);     %mtScal_FCU;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 29) = dcdmtScal(2);       %mtScal_FCUat N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 30) = dcdmtScal(3);     %mtScal_ECRL;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 30) = dcdmtScal(3);       %mtScal_ECRL at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 31) = dcdmtScal(4);     %mtScal_ECRB;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 31) = dcdmtScal(4);       %mtScal_ECRB at N+1;
    
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 32) = dcdmtScal(5);     %mtScal_ECU;
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 32) = dcdmtScal(5);       %mtScal_ECU at N+1;
    
    % w.r.t afactor
    dfdafactor = afactor_pd(q,Lmo,lts,fmo,phi,mtScal,a,para(21));
                 
    J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 33) = dfdafactor;         %afactro at N
    J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 33) = dfdafactor;         %afactor at N+1
    
    % Contraction dynamics; dfdx
    J((Nstate+Npara)*(n-1)+3,(n-1)*(Nstate+Nctrl) + 3)  =  -1/h + 1/2*(u(1)/tact + (1-u(1))/tdeact); %a1
    J((Nstate+Npara)*(n-1)+3,(n)*(Nstate+Nctrl)   + 3)  =  1/h  + 1/2*(u(1)/tact + (1-u(1))/tdeact); %a1 at N+1
    
    J((Nstate+Npara)*(n-1)+4,(n-1)*(Nstate+Nctrl) + 4)  =  -1/h + 1/2*(u(2)/tact + (1-u(2))/tdeact); %a2
    J((Nstate+Npara)*(n-1)+4,(n)*(Nstate+Nctrl)   + 4)  =  1/h  + 1/2*(u(2)/tact + (1-u(2))/tdeact); %a2 at N+1
    
    J((Nstate+Npara)*(n-1)+5,(n-1)*(Nstate+Nctrl) + 5)  =  -1/h + 1/2*(u(3)/tact + (1-u(3))/tdeact); %a3
    J((Nstate+Npara)*(n-1)+5,(n)*(Nstate+Nctrl)   + 5)  =   1/h + 1/2*(u(3)/tact + (1-u(3))/tdeact); %a3 at N+1
    
    J((Nstate+Npara)*(n-1)+6,(n-1)*(Nstate+Nctrl) + 6)  =  -1/h + 1/2*(u(4)/tact + (1-u(4))/tdeact); %a4
    J((Nstate+Npara)*(n-1)+6,(n)*(Nstate+Nctrl)   + 6)  =  1/h  + 1/2*(u(4)/tact + (1-u(4))/tdeact); %a4 at N+1
    
    J((Nstate+Npara)*(n-1)+7,(n-1)*(Nstate+Nctrl) + 7)  =  -1/h + 1/2*(u(5)/tact + (1-u(5))/tdeact); %a5
    J((Nstate+Npara)*(n-1)+7,(n)*(Nstate+Nctrl)   + 7)  =  1/h  + 1/2*(u(5)/tact + (1-u(5))/tdeact); %a5 at N+1
    
    %% Muscle contraction dynamics;
    J((Nstate+Npara)*(n-1)+3,(n-1)*(Nstate+Nctrl) + 8)  = ((u(1)-1)/(2*tdeact) - u(1)/(2*tact)) + (1/(2*tact) - 1/(2*tdeact))*(a(1)-u(1)); % u1 at N
    J((Nstate+Npara)*(n-1)+3,(n)*(Nstate+Nctrl) + 8)    = ((u(1)-1)/(2*tdeact) - u(1)/(2*tact)) + (1/(2*tact) - 1/(2*tdeact))*(a(1)-u(1));% u1 at N +1
    
    J((Nstate+Npara)*(n-1)+4,(n-1)*(Nstate+Nctrl) + 9)  = ((u(2)-1)/(2*tdeact) - u(2)/(2*tact)) + (1/(2*tact) - 1/(2*tdeact))*(a(2)-u(2)); % u2 at N
    J((Nstate+Npara)*(n-1)+4,(n)*(Nstate+Nctrl) + 9)    = ((u(2)-1)/(2*tdeact) - u(2)/(2*tact)) + (1/(2*tact) - 1/(2*tdeact))*(a(2)-u(2)); % u2 at N+1
    
    J((Nstate+Npara)*(n-1)+5,(n-1)*(Nstate+Nctrl) + 10) = ((u(3)-1)/(2*tdeact) - u(3)/(2*tact)) + (1/(2*tact) - 1/(2*tdeact))*(a(3)-u(3)); % u3 at N
    J((Nstate+Npara)*(n-1)+5,(n)*(Nstate+Nctrl) + 10)   = ((u(3)-1)/(2*tdeact) - u(3)/(2*tact)) + (1/(2*tact) - 1/(2*tdeact))*(a(3)-u(3)); % u3 at N+1
    
    J((Nstate+Npara)*(n-1)+6,(n-1)*(Nstate+Nctrl) + 11) = ((u(4)-1)/(2*tdeact) - u(4)/(2*tact)) + (1/(2*tact) - 1/(2*tdeact))*(a(4)-u(4)); % u4 at N
    J((Nstate+Npara)*(n-1)+6,(n)*(Nstate+Nctrl) + 11)   = ((u(4)-1)/(2*tdeact) - u(4)/(2*tact)) + (1/(2*tact) - 1/(2*tdeact))*(a(4)-u(4)); % u4 at N+1
    
    J((Nstate+Npara)*(n-1)+7,(n-1)*(Nstate+Nctrl) + 12) = ((u(5)-1)/(2*tdeact) - u(5)/(2*tact)) + (1/(2*tact) - 1/(2*tdeact))*(a(5)-u(5)); % u5 at N
    J((Nstate+Npara)*(n-1)+7,(n)*(Nstate+Nctrl) + 12)   = ((u(5)-1)/(2*tdeact) - u(5)/(2*tact)) + (1/(2*tact) - 1/(2*tdeact))*(a(5)-u(5)); % u5 at N+1
    
    %% Constraints w.r.t MT parameters themself to maintain they are constant all the time
    for num = 1:length(para)
        J((Nstate+Npara)*(n-1)+7+num,(n-1)*(Nstate+Nctrl) + 12+num) = -1;
        J((Nstate+Npara)*(n-1)+7+num,(n)*(Nstate+Nctrl) + 12+num ) = 1;
    end
    
end

% Task constraints
J((gridN-1)*(Nstate+Npara)+1,1) = 1;                                 % initial position
J((gridN-1)*(Nstate+Npara)+2,2) = 1;                                 % initial joint velocity == 0;


% initial activation state
J((gridN-1)*(Nstate+Npara)+3,3)   =  1;
J((gridN-1)*(Nstate+Npara)+4,4)   =  1;
J((gridN-1)*(Nstate+Npara)+5,5)   =  1;
J((gridN-1)*(Nstate+Npara)+6,6)   =  1;
J((gridN-1)*(Nstate+Npara)+7,7)  =  1;

J((gridN-1)*(Nstate+Npara)+3,8)   = -1;
J((gridN-1)*(Nstate+Npara)+4,9)   = -1;
J((gridN-1)*(Nstate+Npara)+5,10)  = -1;
J((gridN-1)*(Nstate+Npara)+6,11)  = -1;
J((gridN-1)*(Nstate+Npara)+7,12) = -1;
J = sparse(J);
end