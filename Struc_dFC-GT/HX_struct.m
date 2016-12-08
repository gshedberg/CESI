function [HotOut,ColdOut] = HX_struct(ColdIn,HotIn,effectiveness)
[~,H_coldIn] = enthalpy2(ColdIn);
ColdInIdeal = ColdIn;
ColdInIdeal.T = HotIn.T;
[~,H_coldIdeal] = enthalpy2(ColdInIdeal);
H_ColdOut = H_coldIn+ effectiveness.*(H_coldIdeal-H_coldIn);
ColdOut = ColdIn;
ColdOut.T = ColdIn.T + 100;
error = 1;
while abs(error)>1e-3
    [~,H_guess] = enthalpy2(ColdOut);
    Cp = SpecHeat2(ColdOut);
    error = (H_guess-H_ColdOut)./(Cp.*NetFlow(ColdIn));
    ColdOut.T = ColdOut.T - error;     
end
[~,H_Hotin] = enthalpy2(HotIn);
H_hotOut = H_Hotin - effectiveness.*(H_coldIdeal-H_coldIn);
HotOut = HotIn;
HotOut.T = HotIn.T + 100;
T_error = 1;
while abs(T_error)>1e-3
    [~,H_guess] = enthalpy2(HotOut);
    Cp = SpecHeat2(HotOut);
    T_error = (H_guess-H_hotOut)./(Cp.*NetFlow(HotIn));
    HotOut.T = HotOut.T - T_error;     
end