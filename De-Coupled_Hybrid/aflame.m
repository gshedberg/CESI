function [Tout,afr] = aflame(TXNinAir,TXNinFuel)
erxn1 = 1;
erxn2 = 1;
erxn3 = 1;
h = enthalpy(1000);
h_f = enthalpy(298);
TinAir = TXNinAir(:,1);
XinAir = TXNinAir(:,2:8);
NinAir = TXNinAir(:,9);
TinFuel = TXNinFuel(:,1);
XinFuel = TXNinFuel(:,2:8);
NinFuel = TXNinFuel(:,9);
[~,hinFuel] = enthalpy(TinFuel,XinFuel,NinFuel);
[~,hinair] = enthalpy(TinAir,XinAir,NinAir);

% LHVCO = 303000; %Lower Heating Value of CO
% LHVH2 = 240420; %Lower HEating Value of H2
% LHVCH4 = ones(length(TXNinAir(:,1)),1).*8e5;
% LHVfuel = zeros(length(TXNinFuel(:,1)),1)+LHVCH4.*TXNinFuel(:,2)+LHVCO.*TXNinFuel(:,3)+LHVH2.*TXNinFuel(:,5);%LHV of anode off products

%Reactions occuring in combustion%
hrxn1 = 2*h(:,5)+h(:,2)-h(:,1)-1.5*h(:,7); %Ch4 + 1.5 O2 --> CO + 2 H2O
hrxn2 = h(:,3)-h(:,2)-.5*h(:,7); %CO + .5 O2 --> CO2 
hrxn3 = h(:,5)-h(:,4)-.5*h(:,7); %H2 + .5 O2 -->  H2O
%Combustion Components
R1 =(XinFuel(:,1).*(NinFuel)*erxn1);
R2 = (XinFuel(:,2).*(NinFuel)+R1)*erxn2;
R3 =(XinFuel(:,4).*(NinFuel))*erxn3;

%Total Enthalpy of Stoich Combustion
hrxn = (R1.*hrxn1)-(R2.*hrxn2)-(R3.*hrxn3);



T_guess = ones(length(TinAir),1).*1000;
T_error = 100;
while abs(T_error) > 1
    [~,h_guess] = enthalpy(T_guess,XinFuel,NinFuel);
    Cp = SpecHeat(T_guess,XinAir);
    T_error = (hrxn-h_guess)./(Cp.*NinFuel);
    T_guess = T_guess + T_error;  
end

N_guess = NinAir;
N2_error = 100;
while max(abs(N2_error)) > .00001
    [~,h_lean] = enthalpy(T_guess,XinAir,N_guess);
    Cp2 = SpecHeat(T_guess,XinAir);
    N2_error = (hrxn-h_guess)./(Cp2.*N_guess);
    N_guess = N_guess+N2_error;
    T2_error = (h_lean-h_guess)./(Cp2.*N_guess);
    T_out = T_guess-T2_error;
end




