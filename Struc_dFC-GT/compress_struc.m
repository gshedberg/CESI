function [Wc, Tout, P2] = compress_struc(Flow, eff, Pr, Pin)
%% Inlet Conditions
Ru = 8.314;
P2 = Pr.*Pin;
[~,Hin] = enthalpy2(Flow);
Tin = Flow.T;
Tavg = (Tin+Pr.^(1-1/1.4).*Tin)/2;
Flow.T = Tavg;
Cp = SpecHeat2(Flow);
gam = Cp./(Cp-Ru);
FlowS = Flow;
FlowS.T = Pr.^(1-(1./gam)).*Tin; %Isnetropic Temperature Change
[~,Hs] = enthalpy2(FlowS); %Isentropic Enthalpy change
%% Actual Enthalpy Change
Hout = Hin + 1/eff*(Hs-Hin);
Tout = Tin + 1/eff*(FlowS.T-Tin);
T_error = Tin*0+ 100;
%% Solve for Temp out
while min(abs(T_error)) > .1
    Flow.T = Tout;
    [~,Hguess] = enthalpy2(Flow);
    T_error = (Hout - Hguess)./(Cp.*NetFlow(Flow));
    Tout  = Tout + T_error; %Reiteration to calculate temperature out
end
%% Power Out
Wc = Hout - Hin; %power taken for compression