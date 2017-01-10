function [i,r,FuelFlow,FlowOut,V,P,Utilization] = FuelCell_struc(T,ASR,e2,S2C,Oxidant,L,W,n,Cells,Pr)
Ru = 8.314;
F = 96485; %Faraday constant C/mol
Oxidant.O2 = Oxidant.O2./Cells;
Oxidant.T = T;
J = 4000*F*Oxidant.O2; % total current in A/cell
V = .85;
%%Initial guesses
r = .5;
h = enthalpy2(298);
hrxn1 = 2*h.H2O-2*h.H2-h.O2; %2H2 + O2 -->  2H2O
hrxn2 = h.CO2+h.H2-h.CO-h.H2O; %CO + H20 --> CO2 + H2 
hrxn3 = 3*h.H2+h.CO-h.CH4-h.H2O; %CH4+H2O --> CO + 3H2
Qgen = -hrxn1*Oxidant.O2 - J*V/1000; %total heat released in kW/cell
Qcool = hrxn3 + e2*hrxn2; %cooling in kJ/kmol of CH4
Fuel.CH4 = Qgen./Qcool; %flow rate in kmol/sec-cell
h = enthalpy2(T);
s = entropy2(T);
E0 = -((h.H2O-s.H2O.*T)-(h.H2-s.H2.*T)-.5*(h.O2-s.O2.*T))/(2*F); %reference voltage
%solve for distribution
i = (ones(n,1)).*(J./(L.*W)); %Initial current distribution per cell
J_int = zeros(n,1);
error = 1;
while abs(error)>1e-3
    r = S2C/((.5*S2C-1)*(1+e2)+2*Oxidant.O2/Fuel.CH4);
    for k=1:1:n
        J_int(k) = sum(i(1:k))*(W*L/n); %integral of current density as a function of length, total current thus far
    end

    X_H2 = 1+e2/3-2*Oxidant.O2*r/(3*Fuel.CH4)-(1-r)/(6000*F*Fuel.CH4)*J_int;
    X_H2O = 2*Oxidant.O2*r/(3*Fuel.CH4) - (1+e2)/3 + (1-r)/(6000*F*Fuel.CH4)*J_int;
    X_CO_L = (1-e2)/3;
    X_CO2_L = 1 - X_CO_L - X_H2(end) - X_H2O(end);
    SteamRatio = X_H2O(end)*3*Fuel.CH4*r/(1-r)/(Fuel.CH4+0.5*X_CO_L*3*Fuel.CH4*r/(1-r));
    
%     E = E0 + Ru*T/(2*F)*log((1/Pr).*(X_H2./X_H2O.^.5));%Nernst Potential
    E = E0 + Ru*T/(2*F)*log(X_H2./X_H2O.*Pr.^.5);%Nernst Potential
    error2 = 1;
    count=0;
    if error == 1
        V = mean(E - i*ASR);
    end
    while abs(mean(error2))>(J*1e-6) || count < 2
        i = max(0,(E-V)/ASR);%new current distribution
        error2 = (sum(i)/n*L*W) - J;%error in total current
        V = V + .5*(error2/(L*W)*ASR) ;
        count = count + 1;
    end
    Qgen = -J/(4000*F)*hrxn1 - V*J/1000;%heat release from electrochemistry
    NewFuel = Qgen/(e2*hrxn2+hrxn3); %energy balance heat generated = reformer cooling
    
    error = (Fuel.CH4-NewFuel)/Fuel.CH4; %change in fuel estimation on this iteration
    Fuel.CH4 = .8*Fuel.CH4+.2*NewFuel;
    
end
FuelFlow = Fuel.CH4*Cells;
FlowOut.T = T;
FlowOut.H2 = X_H2(end)*3*Fuel.CH4*Cells;
FlowOut.H2O = X_H2O(end)*3*Fuel.CH4*Cells;
FlowOut.CO = X_CO_L*3*Fuel.CH4*Cells;
FlowOut.CO2 = X_CO2_L*3*Fuel.CH4*Cells;
FlowOut.CH4 = 0;
Utilization = J/(2000*F)/(4*Fuel.CH4); %actual H2 use in kmol/s divided by ideal H2 production in kmol/s
<<<<<<< HEAD
P = (V*J)/1000*Cells;
=======
P = (V*J)/1000*Cells;

function r = solveRecirc(Oxidant,e2,Fuel,S2C,r)
error = 1;
while abs(error)>1e-2
    S2Cguess1 = (2*Oxidant.O2-(1+e2)*Fuel.CH4)*r/(1-r)/Fuel.CH4;
    S2Cguess=r.*((.5*S2Cguess1-1)*(1+e2)+2*Oxidant.O2/Fuel.CH4);
    error = S2C - S2Cguess;
    r = r + .05*error;
    
end
>>>>>>> origin/master
