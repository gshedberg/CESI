function [ThotOut,TcoldOut] = HX(ColdIn,HotIn,effectiveness)
[~,h_coldIn] = enthalpy(ColdIn(:,1),ColdIn(:,2:8),ColdIn(:,9));
[~,h_coldIdeal] = enthalpy(HotIn(:,1),ColdIn(:,2:8),ColdIn(:,9));

h_coldOut = h_coldIn+ effectiveness.*(h_coldIdeal-h_coldIn);
TcoldOut = ColdIn(:,1)+100;
error = 1;
while abs(error)>1e-3
    [~,h_guess] = enthalpy(TcoldOut,ColdIn(:,2:8),ColdIn(:,9));
    Cp = SpecHeat(TcoldOut,ColdIn(:,2:8));
    error = (h_guess-h_coldOut)./(Cp.*ColdIn(:,9));
    TcoldOut = TcoldOut - error;     
end
[~,h_Hotin] = enthalpy(HotIn(:,1),HotIn(:,2:8),HotIn(:,9));
h_hotOut = h_Hotin - effectiveness.*(h_coldIdeal-h_coldIn);
ThotOut = HotIn(:,1)+100;
T_error = 1;
while abs(T_error)>1e-3
    [~,h_guess] = enthalpy(ThotOut,HotIn(:,2:8),HotIn(:,9));
    Cp = SpecHeat(ThotOut,HotIn(:,2:8));
    T_error = (h_guess-h_hotOut)./(Cp.*HotIn(:,9));
    ThotOut = ThotOut - T_error;     
end