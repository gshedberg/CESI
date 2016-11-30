function [i,r,FuelFlow,FlowOut,V,P] = FuelCell_struc(T,ASR,e2,S2C,Oxidant,L,W,n,Cells,Pr)
Ru = 8.314;
F = 96485; %Faraday constant C/mol
Oxidant.O2 = Oxidant.O2./Cells;
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
i = (ones(n,length(J))).*(J./(L.*W)); %Initial current distribution per cell
i_int = zeros(n,1);
error = 1;
while abs(error)>1e-3
     r = solveRecirc(Oxidant,e2,Fuel,S2C,r);
    for k=1:1:n
        i_int(k) = sum(i(1:k))*(L/n); %integral of current density as a function of length
    end
    X_H2 = 1+e2/3-2*Oxidant.O2*r/(3*Fuel.CH4)-W*(1-r)/(6000*F*Fuel.CH4)*i_int;
    X_H2O = 2*Oxidant.O2*r/(3*Fuel.CH4) - (1+e2)/3 + W*(1-r)/(6000*F*Fuel.CH4)*i_int;
    E = E0 - Ru*T/F*log((1./Pr).*X_H2O./(X_H2.*ones(n,1).^.5));%Nernst Potential
%     V =(mean(E)-mean(i)*ASR);
    error2 = 1;
    count=0;
    while abs(mean(error2))>(J*1e-4) || count<2
        i = (E-V)/ASR;%new current distribution
        error2 = (sum(i)*L*W/n) - J/n;
        V = V + error2/(L*W)*ASR ;
        count = count+1;
    end
    Qgen = -J/(4000*F)*hrxn1 - V*J/1000;%heat release from electrochemistry
    H2 = enthalpy2(T);
    H_298 = enthalpy2(298);
    X_CO2_L = e2*(1-r)/3;
    X_CO_L = 1 - X_CO2_L - X_H2(end) - X_H2O(end);
    Flow.T = T;
    Flow.H2 = X_H2(end)*3*Fuel.CH4/(1-r);
    Flow.H2O = X_H2O(end)*3*Fuel.CH4/(1-r);
    Flow.CO = X_CO_L*3*Fuel.CH4/(1-r);
    Flow.CO2 = X_CO2_L*3*Fuel.CH4/(1-r);
    Hout = enthalpy2(Flow);
    h_out = Hout/NetFlow(Flow);
%     NewFuel = (-Qgen-Oxidant.O2*h.O2)/(h_298.CH4 - 3*h_out - (e2*hrxn2+hrxn3)); %energy balance solved to determine fuel input
    NewFuel = (Oxidant.O2*H2.O2-Qgen)/(H_298.CH4 + 3*h_out - (e2*hrxn2+hrxn3)); %energy balance solved to determine fuel input

    error = (Fuel.CH4-NewFuel)/Fuel.CH4; %change in fuel estimation on this iteration
    Fuel.CH4 = .6*Fuel.CH4+.4*NewFuel;
end
FuelFlow = Fuel.CH4*Cells;
FlowOut.T = T;
FlowOut.H2 = Flow.H2*Cells;
FlowOut.H2O = Flow.H2O*Cells;
FlowOut.CO = Flow.CO*Cells;
FlowOut.CO2 = Flow.CO2*Cells;
FlowOut.CH4 = FuelFlow;
P = (V*J)/1000*Cells;

function r = solveRecirc(Oxidant,e2,Fuel,S2C,r)
error = 1;
while abs(error)>1e-2
    S2Cguess = ((2*Oxidant.O2*r-(1+e2)*Fuel.CH4)/(1-r) + 2*Oxidant.O2)*r/Fuel.CH4;
    error = S2C - S2Cguess;
    r = r + .05*error;
end