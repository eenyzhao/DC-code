function dfda = dtauda(q,Lmo,lts,phi,fmo,mtScal,afactor,a) 
%% Step One --> Normlaise lm
lambda   = 0.15;  % constant in force length relationship with respect to activation

% FCR
S(1) = -fmo(1)*mtScal(1)*((afactor*exp(afactor*(a(1)))*exp(-((4*(Lmo(1)^2*sin(phi(1))^2 + (lts(1) - mtScal(1)*((1423*cos((1937*q)/5000))/100000 - (1927*sin((1937*q)/5000))/50000 + 1491/5000))^2)^(1/2))/(Lmo(1)*(lambda*((exp(afactor*(a(1))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)^2))/(2*(exp(afactor) - 1)) + (4*afactor*lambda*exp(afactor*(a(1)))*exp(-((4*(Lmo(1)^2*sin(phi(1))^2 + (lts(1) - mtScal(1)*((1423*cos((1937*q)/5000))/100000 - (1927*sin((1937*q)/5000))/50000 + 1491/5000))^2)^(1/2))/(Lmo(1)*(lambda*((exp(afactor*(a(1))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)^2)*((4*(Lmo(1)^2*sin(phi(1))^2 + (lts(1) - mtScal(1)*((1423*cos((1937*q)/5000))/100000 - (1927*sin((1937*q)/5000))/50000 + 1491/5000))^2)^(1/2))/(Lmo(1)*(lambda*((exp(afactor*(a(1))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)*(Lmo(1)^2*sin(phi(1))^2 + (lts(1) - mtScal(1)*((1423*cos((1937*q)/5000))/100000 - (1927*sin((1937*q)/5000))/50000 + 1491/5000))^2)^(1/2)*(exp(afactor*(a(1))) - 1))/(Lmo(1)*(exp(afactor) - 1)^2*(lambda*((exp(afactor*(a(1))) - 1)/(exp(afactor) - 1) - 1) - 1)^2))*(1 - (Lmo(1)^2*sin(phi(1))^2)/(Lmo(1)^2*sin(phi(1))^2 + (lts(1) - mtScal(1)*((1423*cos((1937*q)/5000))/100000 - (1927*sin((1937*q)/5000))/50000 + 1491/5000))^2))^(1/2)*((3732599*cos((1937*q)/5000))/250000000 + (2756351*sin((1937*q)/5000))/500000000);

% FCU
S(2) = -fmo(2)*mtScal(2)*((afactor*exp(afactor*(a(2)))*exp(-((4*(Lmo(2)^2*sin(phi(2))^2 + (lts(2) - mtScal(2)*((2491*cos((1227*q)/5000))/50000 - (1219*sin((1227*q)/5000))/20000 + 657/2500))^2)^(1/2))/(Lmo(2)*(lambda*((exp(afactor*(a(2))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)^2))/(2*(exp(afactor) - 1)) + (4*afactor*lambda*exp(afactor*(a(2)))*exp(-((4*(Lmo(2)^2*sin(phi(2))^2 + (lts(2) - mtScal(2)*((2491*cos((1227*q)/5000))/50000 - (1219*sin((1227*q)/5000))/20000 + 657/2500))^2)^(1/2))/(Lmo(2)*(lambda*((exp(afactor*(a(2))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)^2)*((4*(Lmo(2)^2*sin(phi(2))^2 + (lts(2) - mtScal(2)*((2491*cos((1227*q)/5000))/50000 - (1219*sin((1227*q)/5000))/20000 + 657/2500))^2)^(1/2))/(Lmo(2)*(lambda*((exp(afactor*(a(2))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)*(Lmo(2)^2*sin(phi(2))^2 + (lts(2) - mtScal(2)*((2491*cos((1227*q)/5000))/50000 - (1219*sin((1227*q)/5000))/20000 + 657/2500))^2)^(1/2)*(exp(afactor*(a(2))) - 1))/(Lmo(2)*(exp(afactor) - 1)^2*(lambda*((exp(afactor*(a(2))) - 1)/(exp(afactor) - 1) - 1) - 1)^2))*(1 - (Lmo(2)^2*sin(phi(2))^2)/(Lmo(2)^2*sin(phi(2))^2 + (lts(2) - mtScal(2)*((2491*cos((1227*q)/5000))/50000 - (1219*sin((1227*q)/5000))/20000 + 657/2500))^2))^(1/2)*((1495713*cos((49*q)/200))/100000000 + (3056457*sin((1227*q)/5000))/250000000);

% ECRL
S(3) = fmo(3)*mtScal(3)*((10071157*cos((7499*q)/100000))/1000000000 - (34952839*sin((7499*q)/100000))/1000000000)*((afactor*exp(afactor*(a(3)))*exp(-((4*(Lmo(3)^2*sin(phi(3))^2 + (lts(3) - mtScal(3)*((4661*cos((7499*q)/100000))/10000 + (1343*sin((7499*q)/100000))/10000 - 1317/10000))^2)^(1/2))/(Lmo(3)*(lambda*((exp(afactor*(a(3))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)^2))/(2*(exp(afactor) - 1)) + (4*afactor*lambda*exp(afactor*(a(3)))*exp(-((4*(Lmo(3)^2*sin(phi(3))^2 + (lts(3) - mtScal(3)*((4661*cos((7499*q)/100000))/10000 + (1343*sin((7499*q)/100000))/10000 - 1317/10000))^2)^(1/2))/(Lmo(3)*(lambda*((exp(afactor*(a(3))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)^2)*((4*(Lmo(3)^2*sin(phi(3))^2 + (lts(3) - mtScal(3)*((4661*cos((7499*q)/100000))/10000 + (1343*sin((7499*q)/100000))/10000 - 1317/10000))^2)^(1/2))/(Lmo(3)*(lambda*((exp(afactor*(a(3))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)*(Lmo(3)^2*sin(phi(3))^2 + (lts(3) - mtScal(3)*((4661*cos((7499*q)/100000))/10000 + (1343*sin((7499*q)/100000))/10000 - 1317/10000))^2)^(1/2)*(exp(afactor*(a(3))) - 1))/(Lmo(3)*(exp(afactor) - 1)^2*(lambda*((exp(afactor*(a(3))) - 1)/(exp(afactor) - 1) - 1) - 1)^2))*(1 - (Lmo(3)^2*sin(phi(3))^2)/(Lmo(3)^2*sin(phi(3))^2 + (lts(3) - mtScal(3)*((4661*cos((7499*q)/100000))/10000 + (1343*sin((7499*q)/100000))/10000 - 1317/10000))^2))^(1/2);

% ECRB
S(4) = fmo(4)*mtScal(4)*((afactor*exp(afactor*(a(4)))*exp(-((4*(Lmo(4)^2*sin(phi(4))^2 + (lts(4) - mtScal(4)*((4893*cos((1069*q)/5000))/100000 + (311*sin((1069*q)/5000))/5000 + 2371/10000))^2)^(1/2))/(Lmo(4)*(lambda*((exp(afactor*(a(4))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)^2))/(2*(exp(afactor) - 1)) + (4*afactor*lambda*exp(afactor*(a(4)))*exp(-((4*(Lmo(4)^2*sin(phi(4))^2 + (lts(4) - mtScal(4)*((4893*cos((1069*q)/5000))/100000 + (311*sin((1069*q)/5000))/5000 + 2371/10000))^2)^(1/2))/(Lmo(4)*(lambda*((exp(afactor*(a(4))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)^2)*((4*(Lmo(4)^2*sin(phi(4))^2 + (lts(4) - mtScal(4)*((4893*cos((1069*q)/5000))/100000 + (311*sin((1069*q)/5000))/5000 + 2371/10000))^2)^(1/2))/(Lmo(4)*(lambda*((exp(afactor*(a(4))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)*(Lmo(4)^2*sin(phi(4))^2 + (lts(4) - mtScal(4)*((4893*cos((1069*q)/5000))/100000 + (311*sin((1069*q)/5000))/5000 + 2371/10000))^2)^(1/2)*(exp(afactor*(a(4))) - 1))/(Lmo(4)*(exp(afactor) - 1)^2*(lambda*((exp(afactor*(a(4))) - 1)/(exp(afactor) - 1) - 1) - 1)^2))*(1 - (Lmo(4)^2*sin(phi(4))^2)/(Lmo(4)^2*sin(phi(4))^2 + (lts(4) - mtScal(4)*((4893*cos((1069*q)/5000))/100000 + (311*sin((1069*q)/5000))/5000 + 2371/10000))^2))^(1/2)*((332459*cos((1069*q)/5000))/25000000 - (5230617*sin((1069*q)/5000))/500000000);

% ECU
S(5) = fmo(5)*mtScal(5)*(1 - (Lmo(5)^2*sin(phi(5))^2)/(Lmo(5)^2*sin(phi(5))^2 + (lts(5) - mtScal(5)*((3160734304879671*sin((1213*q)/1000))/576460752303423488 - (646097211181677*cos((1213*q)/1000))/4611686018427387904 + 581/2000))^2))^(1/2)*((afactor*exp(afactor*(a(5)))*exp(-((4*(Lmo(5)^2*sin(phi(5))^2 + (lts(5) - mtScal(5)*((3160734304879671*sin((1213*q)/1000))/576460752303423488 - (646097211181677*cos((1213*q)/1000))/4611686018427387904 + 581/2000))^2)^(1/2))/(Lmo(5)*(lambda*((exp(afactor*(a(5))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)^2))/(2*(exp(afactor) - 1)) + (4*afactor*lambda*exp(afactor*(a(5)))*exp(-((4*(Lmo(5)^2*sin(phi(5))^2 + (lts(5) - mtScal(5)*((3160734304879671*sin((1213*q)/1000))/576460752303423488 - (646097211181677*cos((1213*q)/1000))/4611686018427387904 + 581/2000))^2)^(1/2))/(Lmo(5)*(lambda*((exp(afactor*(a(5))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)^2)*((4*(Lmo(5)^2*sin(phi(5))^2 + (lts(5) - mtScal(5)*((3160734304879671*sin((1213*q)/1000))/576460752303423488 - (646097211181677*cos((1213*q)/1000))/4611686018427387904 + 581/2000))^2)^(1/2))/(Lmo(5)*(lambda*((exp(afactor*(a(5))) - 1)/(exp(afactor) - 1) - 1) - 1)) + 4)*(Lmo(5)^2*sin(phi(5))^2 + (lts(5) - mtScal(5)*((3160734304879671*sin((1213*q)/1000))/576460752303423488 - (646097211181677*cos((1213*q)/1000))/4611686018427387904 + 581/2000))^2)^(1/2)*(exp(afactor*(a(5))) - 1))/(Lmo(5)*(exp(afactor) - 1)^2*(lambda*((exp(afactor*(a(5))) - 1)/(exp(afactor) - 1) - 1) - 1)^2))*((1872056011630391*cos((1213*q)/1000))/281474976710656000 + (6122780602838861*sin((1213*q)/1000))/36028797018963968000);


dfda = S;
end
