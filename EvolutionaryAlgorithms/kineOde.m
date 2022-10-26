function dy = kineOde(t,y,tau,bodyMass,bodyHeight,lofHand,fs) 
% hand mass regression
% Y = -0.1165 + 0.0036*bodyweight(kg) + 0.00175*height(cm);

m = bodyMass*0.0065;   % Kg
g = 9.81;
%CoM : l = 4.11 + 0.026*bodyweight(kg) + 0.033*height(cm) from 3rd dactylion
l = 0.08;   % hand length

% Inertia parameters
J = (-13.68 + 0.088*bodyMass + 0.092*bodyHeight)/10000;
J = J + m*lofHand^2;
% J = m *(l^2)*1/3 + 0.0134;
B = 0.3; % damping parameter

Fs = 1000;
k = floor(1+t*fs);

% y,dy/dt => y(1),y(2)
dy = zeros(2,1);
dy(1) = y(2);
dy(2) = (tau(k)- B*y(2)- 0.6*y(1))/J;
end