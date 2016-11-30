function [Wt,T9] = turbine_struc(Flow,EffTurb, Pr)
%% Inlet Conditions
Ru = 8.314;
Tin = Flow.T;
Tavg = (Tin+Pr.^(1-1/1.4).*Tin)/2;
Flow.T = Tavg;
Cp = SpecHeat(Flow);
gam = Cp./(Cp-Ru);
Hin = enthalpy2(Flow);
%% Isentropic Expansion
T9s = Pr.^((gam-1)./gam).*Tin;   %Isentropic Expansion Temperature
FlowS = Flow;
FlowS.T = T9s;
H9s = enthalpy2(FlowS);    %Isentropic enthaplpy out
%% Actual Enthalpy Change
H9 = Hin-EffTurb.*(Hin-H9s);           %Actual enthalpy out
T9 = Tin - EffTurb.*(T9s - Tin);
T_error = Tin*0 +100;
%% Solve for Temp Out
while min(abs(T_error)) > .1      %Reiteration to calculate temperature out
    Flow.T = T9;
    H_guess = enthalpy2(Flow);
    T_error = (H_guess-H9)./(Cp.*NetFlow(Flow));
    T9 = T9 - T_error; 
end
%% Power Out
Wt = Hin - H9;       %Isentropic expansion power generation
end

    
    




