function fmt  = updateT(Lm,Lmo,phi,fmo,mtScal,a) 

lambda   = 0.15;  % constant in force length relationship with respect to activation
normLm   = Lm/Lmo;
normLm_A = Lm/(Lmo*(lambda*(1-a)+1));

% update pennation angle
phi_num   = Lmo.*sin(phi);
phi_denum = Lm;
phi_t     = asin(phi_num./phi_denum);



% % active force length relationship
w = 0.25;
FLa_norm =  exp(-((normLm_A-1)/w)^2);

FLp_norm= exp(10 * (normLm - 1 ))/exp(5);

fv =1;
%% Step seven --> calculate the muscle force
fmt = (FLa_norm*a+FLp_norm)*cos(phi_t)*fmo;
end 
