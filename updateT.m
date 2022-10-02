function fmt  = updateT(Lm,Lmo,phi,fmo,mtScal,a) 
%% Step One --> Normalise lm
lambda   = 0.15;  % constant in force length relationship with respect to activation
normLm   = Lm/Lmo;
normLm_A = Lm/(Lmo*(lambda*(1-a)+1));

% update pennation angle
phi_num   = Lmo.*sin(phi);
phi_denum = Lm;
phi_t     = asin(phi_num./phi_denum);


%% step two -->  calcuate the active force and passive force
% % active force length relationship
w = 0.25;
FLa_norm =  exp(-((normLm_A-1)/w)^2);
% FLa_norm = exp(-(normLm_A-1)^2*(0.45^-1));  

% FLa_norm = exp(-(normLm_A./(0.56*(Lmo.*(lambda.*(1-a)+1))).^2));
FLp_norm= exp(10 * (normLm - 1 ))/exp(5);
% kpe = 4;
% e0 = 0.6;
% FLp_norm = exp(((4*(normLm-1))-1)./(exp(normLm)-1));
%% Step four ---> force-velocity curve
% vm = mtScal*vmt*cos(phi_t);
% vm_norm = vm/(10*Lmo);
% d1 = -0.318;
% d2 = -8.149;
% d3 = -0.374;
% d4 = 0.886;
% fv = d1*log((d2*vm_norm+d3)+sqrt((d2*vm_norm+d3)^2+1))+d4; 
fv =1;
%% Step seven --> calculate the muscle force
fmt = (FLa_norm*a+FLp_norm)*cos(phi_t)*fmo;
end 