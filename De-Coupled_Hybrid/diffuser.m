function [P_out, v_out, T_out, X_out, N_out, M] = diffuser(P_in, v_in, TXNin, n_diffuser)
Ru = 8.314;

Cp = SpecHeat(TXNin(:,1), TXNin(:,2:8));
gam = Cp./(Cp-Ru);


[~,Hin] = enthalpy(TXNin(:,1),TXNin(:,2:8),TXNin(:,9));

X_out = TXNin(:,2:8);

c = (gam.*Ru.*TXNin(:,1)).^.5;

M = v_in./c;

N_out = TXNin(:,9);

T_out = zeros(length(v_in),1)+1000;
T_error = 100;
while abs(T_error) > 1
   [~,H_guess] = enthalpy(T_out, X_out, N_out);
   Cp = SpecHeat(T_out, Xout);

   T_error = (Hin-H_guess)./(Cp.*N_out);
   T_out = T_out+T_error;
end

[~,Hout] = enthalpy(T_out, X_out, N_out);

v_out = (v_in.^2 - 2*(Hout - Hin)).^.5;

P_out = P_in.*(T_out./TXNin(:,1).^(gam./gam-1));