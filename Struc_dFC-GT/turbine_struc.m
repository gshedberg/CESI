function [Wt,Tout,Flow] = turbine_struc(Flow,EffTurb, Pr)
%% Inlet Conditions
Ru = 8.314;
Tin = Flow.T;
[~,Hin] = enthalpy2(Flow);
Tavg = (Tin+Pr.^(1-1/1.4).*Tin)/2;
Flow.T = Tavg;
Cp = SpecHeat2(Flow);
gam = Cp./(Cp-Ru);
FlowS = Flow;
FlowS.T = Pr.^((gam-1)./gam).*Tin;   %Isentropic Expansion Temperature
[~,Hs] = enthalpy2(FlowS);    %Isentropic enthaplpy out
%% Actual Enthalpy Change
Hout = Hin-EffTurb.*(Hin-Hs);           %Actual enthalpy out
Tout = Tin - EffTurb.*(FlowS.T - Tin);
T_error = Tin*0 +100;
%% Solve for Temp Out
while min(abs(T_error)) > .1      %Reiteration to calculate temperature out
    Flow.T = Tout;
    [~,H_guess] = enthalpy2(Flow);
    T_error = (H_guess-Hout)./(Cp.*NetFlow(Flow));
    Tout = Tout - T_error; 
end
%% Power Out
Wt = Hin - Hout;       %Isentropic expansion power generation
Flow.T = Tout;