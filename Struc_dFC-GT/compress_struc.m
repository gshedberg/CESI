function [Wc, T2, P2] = compress_struc(Flow, eff, Pr, Pin)
%% Inlet Conditions
Ru = 8.314;
P2 = Pr.*Pin;
H1 = enthalpy2(Flow);
Tin = Flow.T;
Tavg = (Tin+Pr.^(1-1/1.4).*Tin)/2;
Flow.T = Tavg;
Cp = SpecHeat2(Flow);
gam = Cp./(Cp-Ru);
%% Isentropic Compression
T2s = Pr.^(1-(1./gam)).*Tin; %Isnetropic Temperature Change
FlowS = Flow;
FlowS.T = T2s;
H2s = enthalpy2(FlowS); %Isentropic Enthalpy change
%% Actual Enthalpy Change
H2 = H1 + 1/eff*(H2s-H1);
T2 = Tin + 1/eff*(T2s-Tin);
T_error = Tin*0+ 100;
%% Solve for Temp out
while min(abs(T_error)) > .1
    Flow.T = T2;
    Hguess = enthalpy2(Flow);
    T_error = (H2 - Hguess)./(Cp.*NetFlow(Flow));
    T2  = T2 + T_error; %Reiteration to calculate temperature out
end
%% Power Out
Wc = H2 - H1; %power taken for compression

