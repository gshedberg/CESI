function [T_guess] = heatex(TXNin)

Tin = TXNin(:,1);       %Temp in from ITM
Xout = TXNin(:,2:8);       %Composition O2 from ITM
Nout = TXNin(:,9);      %FLow rate O2 from ITM
[~,H3] = enthalpy(Tin,Xout,Nout); %Enthalpy coming in to Heat Exchanger


T_guess = 300;
T_error = 100;
while T_error > .1
    H_guess = enthalpy(T_guess, Xout, Nout);
    Cp = SpecHeat(T_guess, Xout);
    
    T_error = (H_guess - H3)/(Cp*Nout);
    T_guess = T_guess+T_error;
end


