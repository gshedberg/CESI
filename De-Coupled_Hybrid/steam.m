function [Tout,Xout,Nout,W_pump] = steam(TXNinAir,TXNinWater,Eff_pump,Pr)

TinAir = TXNinAir(:,1);
XinAir = TXNinAir(:,2:8);
NinAir = TXNinAir(:,9);
TinWater = TXNinWater(:,1);
XinWater = TXNinWater(:,2:8);
NinWater = TXNinWater(:,9);


[~,h_water] = enthalpy(TinWater,XinWater,NinWater);
[~,h_air] = enthalpy(TinAir,XinAir,NinAir);

H1 = h_air;
H3 = h_water;
T2s = Pr.^(1-(1./gam)).*T1; %Isnetropic Temperature Change in Pump
[~,H2s] = enthalpy(T2s,XinAir,NinAir);
H2 = H1-Eff_pump.*(H1-H2s); 

W_pump = H2-H1;
H_needed = H3-H2; %Heat required to heat up water

T_guess = 300;
T_error = 100;
while T_error > .1
    H_guess = enthalpy(T_guess, Xout, Nout);
    Cp = SpecHeat(T_guess, Xout);
    
    T_error = (H_guess - H3)/(Cp*Nout);
    T_guess = T_guess+T_error;
end
Tout = T_guess;
Xout = Xinair;
Nout = Ninair;
