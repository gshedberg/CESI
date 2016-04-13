function [Wt,T9, N_out] = turbine(TXNin, EffTurb, Pr,Q)

Ru = 8.314;
Tavg = (TXNin(:,1)+Pr.^(1-1/1.4).*TXNin(:,1))/2;
Cp = SpecHeat(Tavg, TXNin(:,2:8));
gam = Cp./(Cp-Ru);

Tin = TXNin(:,1);     %Temperature out of combustor
X_out = TXNin(:,2:8); %Composition coming out of combustor
N_out = TXNin(:,9);  %Flow rate coming out of combustor

T9s = Pr.^((gam-1)./gam).*Tin;   %Isentropic Expansion Temperature

[~,H2_s] = enthalpy(T9s, X_out, N_out);    %Isentropic enthaplpy out
[~, H1] = enthalpy(Tin,X_out,N_out);   %Enthalpy in to turbine from combustor

H2 = H1-EffTurb.*(H1-H2_s);           %Actual enthalpy out
if Pr < 1
    %Heatex
    h_hxout = N_out.*H2-Q;
    T9 = zeros(length(Tin),1)+500;
    T_error = 100;
    while abs(T_error) > .1      %Reiteration to calculate temperature out
        [~,H_guess] = enthalpy(T9,X_out,N_out);
        Cp = SpecHeat(T9,X_out);
        T_error = (H_guess-h_hxout)./(Cp.*N_out);
        T9 = T9 - T_error; 
    end
else
    T9 = zeros(length(Tin),1)+500;
    T_error = 100;
    while abs(T_error) > .1      %Reiteration to calculate temperature out
        [~,H_guess] = enthalpy(T9, X_out, N_out);
        T_error = (H_guess-H2)./(Cp.*N_out);
        T9 = T9 - T_error; 
    end
end
Wt = H1 - H2;       %Isentropic expansion power generation


    
    




