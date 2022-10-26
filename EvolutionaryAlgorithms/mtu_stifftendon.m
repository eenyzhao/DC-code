function [fmt,lm] = mtu_stifftendon(a,Lmt,para,fs)
% scirpt for the muscle-tendon force computation
Lmo = para(1)';   % optimal muscle fibre length
Fmo = para(2)';   % maximum isometric force
Lts = para(3)';   % tendon slack length
phi = para(4)';   % phi

dt = 1/fs;
vmt = diff(Lmt)/dt;

for i = 1:length(Lmt)-2
    % muscle fibre length
    lm(i) = sqrt((Lmo*sin(phi))^2 + (Lmt(i) - Lts)^2) ;
    % Normalzation
    lambda   = 0.15;  % constant in force length relationship with respect to activation
    normLm(i)   = lm(i)/Lmo;
    normLm_A(i) = lm(i)/(Lmo*(lambda*(1-a(i))+1));
    
    % update pennation angle
    phi_num(i)   = Lmo*sin(phi);
    phi_denum (i)= lm(i);
    phi_t(i)     = asin(phi_num(i)/phi_denum(i));
    
    % tendon strain
    % fa
    % fa
    w = 0.3;
    fa(i) = exp(-((normLm_A(i)-1)/w)^2);
    %fp
    fp(i) = exp(10 * (normLm(i) - 1 ))/exp(5);
    %fv
    vm(i) = vmt(i)*cos(phi_t(i));
    
    vmax = 10*Lmo;
    c3 = (vmax*0.25*(1.5-1))/(0.25+1);
    
    if vm(i) <= 0
        fv(i) = (vmax  + vm(i)/vmax)/(vmax - 4*vm(i)/vmax);
    elseif vm(i) > 0
        fv(i) = (1.5*vm(i)/vmax + c3)/(vm(i)/vmax+c3);
    end
        
    fmt(i) = (fa(i)* fv(i)*a(i) + fp(i))*Fmo*cos(phi_t(i));
    
end
    fmt = fmt';
    lm =lm';
end