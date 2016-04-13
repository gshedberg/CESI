function [Tout_stoich] = aflame(TXNinAir,TXNinFuel)
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
[~,H_anode] = enthalpy(TinFuel,XinFuel,NinFuel);
[~,H_resid] = enthalpy(TinAir,XinAir,NinAir);

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

%Total heat of Stoich Combustion
hrxn = (R1.*hrxn1)-(R2.*hrxn2)-(R3.*hrxn3);

%Total Reaction Occuring in Combustor
Nout = (NinFuel + NinAir + .5*R1 - .5*R2 - .5*R3);

Xout(:,1) = ((XinFuel(:,1).*NinFuel + XinAir(:,1).*NinAir) - R1)./Nout;
Xout(:,2) = ((XinFuel(:,2).*NinFuel + XinAir(:,2).*NinAir) + R1-R2)./Nout;
Xout(:,3) = ((XinFuel(:,3).*NinFuel + XinAir(:,3).*NinAir) + R2)./Nout;
Xout(:,4) = ((XinFuel(:,4).*NinFuel + XinAir(:,4).*NinAir) - R3)./Nout;
Xout(:,5) = ((XinFuel(:,5).*NinFuel + XinAir(:,5).*NinAir) + R3 + 2*R1)./Nout;
Xout(:,6) = ((XinFuel(:,6).*NinFuel + XinAir(:,6).*NinAir))./Nout;
Xout(:,7) = ((XinFuel(:,7).*NinFuel + XinAir(:,7).*NinAir) - 1.5*R1 - .5*R2 - .5*R3)./Nout;

%Stoichiometric Reaction
N_O2_Stoich = NinFuel.*(XinFuel(:,4)*.5+XinFuel(:,2)*.5+XinFuel(:,1)*1.5);
N_res_stoich = N_O2_Stoich./XinAir(:,7);
S_F = N_res_stoich./NinAir;
Stoich_flow = S_F.*NinAir;

N_stoich = (NinFuel + Stoich_flow + .5*R1 - .5*R2 - .5*R3);

Xout_St(:,1) = ((XinFuel(:,1).*NinFuel + XinAir(:,1).*(Stoich_flow)) - R1)./N_stoich;
Xout_St(:,2) = ((XinFuel(:,2).*NinFuel + XinAir(:,2).*(Stoich_flow)) + R1-R2)./N_stoich;
Xout_St(:,3) = ((XinFuel(:,3).*NinFuel + XinAir(:,3).*(Stoich_flow)) + R2)./N_stoich;
Xout_St(:,4) = ((XinFuel(:,4).*NinFuel + XinAir(:,4).*(Stoich_flow)) - R3)./N_stoich;
Xout_St(:,5) = ((XinFuel(:,5).*NinFuel + XinAir(:,5).*(Stoich_flow)) + R3 + 2*R1)./N_stoich;
Xout_St(:,6) = ((XinFuel(:,6).*NinFuel + XinAir(:,6).*(Stoich_flow)))./N_stoich;
Xout_St(:,7) = ((XinFuel(:,7).*NinFuel + XinAir(:,7).*(Stoich_flow)) - 1.5*R1 - .5*R2 - .5*R3)./N_stoich;

H_stoich = H_anode+S_F.*H_resid+hrxn;

T_guess = ones(length(TinAir),1).*1000;
error = 100;
while abs(max(error)) > .01
    [~,H_guess] = enthalpy(T_guess,Xout,N_stoich);
    Cp = SpecHeat(T_guess,Xout);
    T_error = (H_stoich-H_guess)./(Cp.*N_stoich);
    T_guess = T_guess + T_error;
    error = H_stoich-H_guess;
end
%Stoich Temperature Flame
Tout_stoich = T_guess;

%Lean Reaction (more air)
N_remain = Nout - N_stoich;

Xout_lean(:,1) = ((XinFuel(:,1).*NinFuel + XinAir(:,1).*(NinAir-Stoich_flow))./N_remain)-(Xout_St(:,1).*N_stoich);
Xout_lean(:,2) = ((XinFuel(:,2).*NinFuel + XinAir(:,2).*(NinAir-Stoich_flow))./N_remain)-(Xout_St(:,2).*N_stoich);
Xout_lean(:,3) = ((XinFuel(:,3).*NinFuel + XinAir(:,3).*(NinAir-Stoich_flow))./N_remain)-(Xout_St(:,3).*N_stoich);
Xout_lean(:,4) = ((XinFuel(:,4).*NinFuel + XinAir(:,4).*(NinAir-Stoich_flow))./N_remain)-(Xout_St(:,4).*N_stoich);
Xout_lean(:,5) = ((XinFuel(:,5).*NinFuel + XinAir(:,5).*(NinAir-Stoich_flow))./N_remain)-(Xout_St(:,5).*N_stoich);
Xout_lean(:,6) = ((XinFuel(:,6).*NinFuel + XinAir(:,6).*(NinAir-Stoich_flow))./N_remain)-(Xout_St(:,6).*N_stoich);
Xout_lean(:,7) = ((XinFuel(:,7).*NinFuel + XinAir(:,7).*(NinAir-Stoich_flow))./N_remain)-(Xout_St(:,7).*N_stoich);

%Total Enthalpy exiting combustor
H_total = H_anode+H_resid+hrxn;
%Remaining enthalpy
H_remain = H_total - H_stoich;
%Lean Temperature
T_guess_2 = ones(length(TinAir),1).*1000;
error_2 = 100;
while abs(max(error_2)) > .01
    [~,H_guess] = enthalpy(T_guess_2,Xout_lean,N_remain);
    Cp = SpecHeat(T_guess_2,Xout_lean);
    T_error_2 = (H_guess - H_remain)./(Cp.*N_remain);
    T_guess_2 = T_guess_2 + T_error_2;
    error_2 = H_guess - H_remain;
end

