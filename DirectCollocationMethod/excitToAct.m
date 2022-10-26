function a = excitToAct(t,u)
 tact = 0.015;  % activation time constant
 tdeact = 0.050; % deactivation time constant
 h = (t(end)-t(1))/(length(t)-1);

 a = zeros(length(t),size(u,2));
 dadt = zeros(length(t)-1,size(u,2));
 
 for ch = 1:size(u,2)
     for iter = 1:length(t)-1
        dadt(iter,ch) = (u(iter,ch)/tact + ((1-u(iter,ch))/tdeact))*(u(iter,ch)-a(iter,ch));
        a(iter+1,ch) = a(iter,ch)+dadt(iter,ch)*h;
     end
 end


end
