function J = conJacobianstructure(aux)
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

 %Size of Jacobian
J = zeros(Ncon,Nx);
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
indpara2      = (n)*(Nstate+Nctrl) + 13 : (n)*(Nstate+Nctrl) + 33;  % msk parameters 

% midpoint rule

% find the Non-zero elements;

%% Equation of motion dfdx
J((Nstate+Npara)*(n-1)+1,indtheta)        = 1;
J((Nstate+Npara)*(n-1)+1,indtheta2)       = 1;
J((Nstate+Npara)*(n-1)+1,indtheta_dot)    = 1;
J((Nstate+Npara)*(n-1)+1,indtheta_dot2)   = 1;

%% Musculoskeletal model constraints; 
% I*ddtheta + B*detheta + K*theta - tau = 0
% 
J((Nstate+Npara)*(n-1)+2,indtheta)        = 1;     % theta n
J((Nstate+Npara)*(n-1)+2,indtheta2)       = 1;     % thata n+1
                                    
J((Nstate+Npara)*(n-1)+2,indtheta_dot)    = 1;  % theta_dot;
J((Nstate+Npara)*(n-1)+2,indtheta_dot2)   = 1;  % theta_dot at n +1;
 
% w.s.t muscle activation
% dcda = dtauda(q,Lmo,lts,phi,fmo,mtScal,para(21),a1,a2);
J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) +3) = 1; % a1
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl)   +3) = 1; % a1 at n+1 

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) +4) = 1; % a2
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl)   +4) = 1; % a2 at n+1

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) +5) = 1; % a3
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl)   +5) = 1; % a3 at n+1

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) +6) = 1; % a4
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl)   +6) = 1; % a4 at n+1

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) +7) = 1; % a5
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl)   +7) = 1; % a5 at n+1


% Partical derivative w.r.t lmo
% dcdlm = lm0_pd(q,Lmo,lts,phi,fmo,mtScal,a,para(21));

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 13) = 1;     %lm0_FCR;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 13) = 1;       %lm0_FCR at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 14) = 1;     %lm0_FCU;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 14) = 1;       %lm0_FCU at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 15) = 1;     %lm0_ECRL;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 15) = 1;       %lm0_ECRL at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 16) = 1;     %lm0_ECRB;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 16) = 1;       %lm0_ECRB at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 17) = 1;     %lm0_ECU;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 17) = 1;       %lm0_ECU at N+1;

% w.r.t fmo
% dcdfmo = fm0_pd(q,Lmo,phi,mtScal,a,para(21));
J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 18) = 1;     %fm0_FCR;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 18) = 1;       %fm0_FCR at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 19) = 1;     %fm0_FCU;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 19) = 1;       %fm0_FCUat N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 20) = 1;     %fm0_ECRL;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 20) = 1;       %fm0_ECRL at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 21) = 1;     %fm0_ECRB;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 21) = 1;       %fm0_ECRB at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 22) = 1;     %fm0_ECU;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 22) = 1;       %fm0_ECU at N+1;

% w.r.t lts
% dcdlts = lts_pd(q,Lmo,lts,fmo,phi,mtScal,a,para(21));
J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 23) = 1;     %lts_FCR;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 23) = 1;       %lts_FCR at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 24) = 1;     %lts_FCU;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 24) = 1;       %lts_FCUat N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 25) = 1;     %lts_ECRL;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 25) = 1;       %lts_ECRL at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 26) = 1;     %lts_ECRB;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 26) = 1;       %lts_ECRB at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 27) = 1;     %lts_ECU;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 27) = 1;       %lts_ECU at N+1;


% w.r.t mtScal
% dcdmtScal = mtScal_pd(q,Lmo,lts,fmo,phi,a,mtScal,para(21));
J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 28) = 1;     %mtScal_FCR;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 28) = 1;       %mtScal_FCR at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 29) = 1;     %mtScal_FCU;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 29) = 1;       %mtScal_FCUat N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 30) = 1;     %mtScal_ECRL;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 30) = 1;       %mtScal_ECRL at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 31) = 1;     %mtScal_ECRB;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 31) = 1;       %mtScal_ECRB at N+1;

J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 32) = 1;     %mtScal_ECU;
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 32) = 1;       %mtScal_ECU at N+1;

% w.r.t afactor
J((Nstate+Npara)*(n-1)+2,(n-1)*(Nstate+Nctrl) + 33) = 1;       %afactro at N
J((Nstate+Npara)*(n-1)+2,(n)*(Nstate+Nctrl) + 33) = 1;       %afactor at N+1

% Contraction dynamics; dfdx
J((Nstate+Npara)*(n-1)+3,(n-1)*(Nstate+Nctrl) + 3)  =  1; %a1 
J((Nstate+Npara)*(n-1)+3,(n)*(Nstate+Nctrl)   + 3)  =  1; %a1 at N+1

J((Nstate+Npara)*(n-1)+4,(n-1)*(Nstate+Nctrl) + 4)  =  1; %a2
J((Nstate+Npara)*(n-1)+4,(n)*(Nstate+Nctrl)   + 4)  =  1; %a2 at N+1

J((Nstate+Npara)*(n-1)+5,(n-1)*(Nstate+Nctrl) + 5)  =  1; %a3
J((Nstate+Npara)*(n-1)+5,(n)*(Nstate+Nctrl)   + 5)  =  1; %a3 at N+1

J((Nstate+Npara)*(n-1)+6,(n-1)*(Nstate+Nctrl) + 6)  =  1; %a4
J((Nstate+Npara)*(n-1)+6,(n)*(Nstate+Nctrl)   + 6)  =  1; %a4 at N+1

J((Nstate+Npara)*(n-1)+7,(n-1)*(Nstate+Nctrl) + 7)  =  1; %a5
J((Nstate+Npara)*(n-1)+7,(n)*(Nstate+Nctrl)   + 7)  =  1; %a5 at N+1

%% Muscle contraction dynamics;
J((Nstate+Npara)*(n-1)+3,(n-1)*(Nstate+Nctrl) + 8)  = 1; % u1 at N                                    
J((Nstate+Npara)*(n-1)+3,(n)*(Nstate+Nctrl) + 8)    = 1;% u1 at N +1

J((Nstate+Npara)*(n-1)+4,(n-1)*(Nstate+Nctrl) + 9)  = 1; % u2 at N
J((Nstate+Npara)*(n-1)+4,(n)*(Nstate+Nctrl) + 9)    = 1; % u2 at N+1

J((Nstate+Npara)*(n-1)+5,(n-1)*(Nstate+Nctrl) + 10) = 1; % u3 at N
J((Nstate+Npara)*(n-1)+5,(n)*(Nstate+Nctrl) + 10)   = 1; % u3 at N+1

J((Nstate+Npara)*(n-1)+6,(n-1)*(Nstate+Nctrl) + 11) = 1; % u4 at N
J((Nstate+Npara)*(n-1)+6,(n)*(Nstate+Nctrl) + 11)   = 1; % u4 at N+1

J((Nstate+Npara)*(n-1)+7,(n-1)*(Nstate+Nctrl) + 12) = 1; % u5 at N
J((Nstate+Npara)*(n-1)+7,(n)*(Nstate+Nctrl) + 12)   = 1; % u5 at N+1

%% Constraints w.r.t MT parameters themself to maintain they are constant all the time
for num = 1:Npara
    J((Nstate+Npara)*(n-1)+7+num,(n-1)*(Nstate+Nctrl) + 12+num) = 1;
    J((Nstate+Npara)*(n-1)+7+num,(n)*(Nstate+Nctrl) + 12+num ) = 1;
end

end
% 
% Task constraints
J((gridN-1)*(Nstate+Npara)+1,1) = 1;                                 % initial position
J((gridN-1)*(Nstate+Npara)+2,2) = 1;                                 % initial joint velocity == 0;

% initial activation state 
J((gridN-1)*(Nstate+Npara)+3,3)  =  1;
J((gridN-1)*(Nstate+Npara)+4,4)  =  1;
J((gridN-1)*(Nstate+Npara)+5,5)  =  1;
J((gridN-1)*(Nstate+Npara)+6,6)  =  1;
J((gridN-1)*(Nstate+Npara)+7,7) =  1;

J((gridN-1)*(Nstate+Npara)+3,8)   = 1;
J((gridN-1)*(Nstate+Npara)+4,9)   = 1;
J((gridN-1)*(Nstate+Npara)+5,10)  = 1;
J((gridN-1)*(Nstate+Npara)+6,11)  = 1;
J((gridN-1)*(Nstate+Npara)+7,12) = 1;
J = sparse(J);
end