function [Tout,Q_transfer] = heatex(Pr,TXN)
Ru = 8.314;

T1 = TXN(:,1); %Inlet Temp'
Xout = TXN(:,2:8); %Inlet Composition
Nout = TXN(:,9); %Inlet Air Flow

Tavg = (T1+Pr.^(1-1/1.4).*T1)/2;
Cp = SpecHeat(Tavg, TXN(:,2:8));
gam = Cp./(Cp-Ru);
    %HeatEx
dT = 283; %K
T_hotin = (1./Pr).^((gam-1)./gam).*1200;   %Isentropic Expansion Temperature
T_out = T_hotin - dT;
[~,H_itmin] = enthalpy(T_out,Xout,Nout);
Q_transfer = Nout.*(H_itmin-H2);
end
