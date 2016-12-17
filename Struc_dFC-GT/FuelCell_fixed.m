function [V] = FuelCell_fixed(T,ASR,e2,S2C,L,W,n,Pr,i_den,util)
Ru = 8.314;
F = 96485; %Faraday constant C/mol
J = i_den*W*L;
Fuel.CH4 = J/(2000*F)/(4*util);
V = .85;
%%Initial guesses
r = .5;
h = enthalpy2(298);
hrxn1 = 2*h.H2O-2*h.H2-h.O2; %2H2 + O2 -->  2H2O
hrxn2 = h.CO2+h.H2-h.CO-h.H2O; %CO + H20 --> CO2 + H2 
hrxn3 = 3*h.H2+h.CO-h.CH4-h.H2O; %CH4+H2O --> CO + 3H2
Qgen = -hrxn1*J/(4000*F) - J*V/1000; %total heat released in kW/cell
Qcool = hrxn3 + e2*hrxn2; %cooling in kJ/kmol of CH4

h = enthalpy2(T);
s = entropy2(T);
E0 = -((h.H2O-s.H2O.*T)-(h.H2-s.H2.*T)-.5*(h.O2-s.O2.*T))/(2*F); %reference voltage
%solve for distribution
i = (ones(n,1)).*(J./(L.*W)); %Initial current distribution per cell
J_int = zeros(n,1);

r = solveRecirc(J,e2,Fuel,S2C,r);
for k=1:1:n
    J_int(k) = sum(i(1:k))*(W*L/n); %integral of current density as a function of length, total current thus far
end
H2flow_reformed = (3+e2)*Fuel.CH4;
H2used = J/(2000*F);
XH2end = (H2flow_reformed-H2used)/(3*Fuel.CH4);
XO2 = [.21;.2;.19;.18;.17;.16;.15;.14;.13;.12];
COflow_reformed = (1-e2)*Fuel.CH4;
XCOend = COflow_reformed/(3*Fuel.CH4);
X_H2 = 1+e2/3-J/(2000*F)*r/(3*Fuel.CH4)-(1-r)/(6000*F*Fuel.CH4)*J_int;
X_H2O = J/(2000*F)*r/(3*Fuel.CH4) - (1+e2)/3 + (1-r)/(6000*F*Fuel.CH4)*J_int;
X_CO2_L = e2*(1-r)/3;
X_CO_L = 1 - X_CO2_L - X_H2(end) - X_H2O(end);
X_CO_L = (1-e2)/3;
X_CO2_L = 1 - X_CO_L - X_H2(end) - X_H2O(end);
% E = E0 + Ru*T/(2*F)*log((1/Pr).*(X_H2./X_H2O.*XO2.^.5));%Nernst Potential
E = E0 + Ru*T/(2*F)*log(X_H2./X_H2O.*XO2.^.5.*Pr^.5);%Nernst Potential
error2 = 1;
count=0;
V = mean(E - i*ASR);

while abs(mean(error2))>(J*1e-6) || count < 2
    i = max(0,(E-V)/ASR);%new current distribution
    error2 = (sum(i)/n*L*W) - J;%error in total current
    V = V + .5*(error2/(L*W)*ASR) ;
    count = count + 1;
end


function r = solveRecirc(J,e2,Fuel,S2C,r)
F = 96485;
error = 1;
while abs(error)>1e-2
    S2Cguess = (J/(2000*F)-(1+e2)*Fuel.CH4)*r/(1-r)/Fuel.CH4;
    error = S2C - S2Cguess;
    r = r + .05*error;
end