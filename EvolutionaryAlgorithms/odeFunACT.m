function dadt = odeFunACT(t,a,e,fs)
    % ode function for muscle activation
    k= floor(1+t*fs);   % iteraction number corresponding to current t 
    tact = 0.015;       % activation time constant
    tdeact = 0.050;     % deactivation time constant  
    dadt = (e(k)/tact + ((1-e(k))/tdeact))*(e(k)-a);
end