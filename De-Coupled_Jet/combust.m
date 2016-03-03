function [Tout, Xout, Nout] = combust(TXNin_air, TXNin_fuel)

erxn1 = 1;
erxn2 = 1;
erxn3 = 1;

h = enthalpy(TXNin_air(:,1)+200);

Tin_air = TXNin_air(:,1);
Xin_air = TXNin_air(:,2:8);
Nin_air = TXNin_air(:,9);

Tin_fuel = TXNin_fuel(:,1);
Xin_fuel = TXNin_fuel(:,2:8);
Nin_fuel = TXNin_fuel(:,9);
%Reactions occuring in combustor%
hrxn1 = 2*h(:,5)+h(:,2)-h(:,1)-1.5*h(:,7); %Ch4 + 1.5 O2 --> CO + 2 H2O
hrxn2 = h(:,3)-h(:,2)-.5*h(:,7); %CO + .5 O2 --> CO2 
hrxn3 = h(:,5)-h(:,4)-.5*h(:,7); %H2 + .5 O2 -->  H2O

%Total Enthalpy of combustion%
R1 =(Xin_fuel(:,1).*(Nin_fuel)*erxn1);
R2 = (Xin_fuel(:,2).*(Nin_fuel)+R1)*erxn2;
R3 =(Xin_fuel(:,4).*(Nin_fuel))*erxn3;

Nout = (TXNin_fuel(:,9) + TXNin_air(:,9) + .5*R1 - .5*R2 - R3);

Xout(:,1) = ((TXNin_fuel(:,2).*Nin_fuel + TXNin_air(:,2).*Nin_air)-R1)./Nout;
Xout(:,2) = ((TXNin_fuel(:,3).*Nin_fuel + TXNin_air(:,3).*Nin_air)+R1-R2)./Nout;
Xout(:,3) = ((TXNin_fuel(:,4).*Nin_fuel + TXNin_air(:,4).*Nin_air)+(R2))./Nout;
Xout(:,4) = ((TXNin_fuel(:,5).*Nin_fuel + TXNin_air(:,5).*Nin_air)-R3)./Nout;
Xout(:,5) = ((TXNin_fuel(:,6).*Nin_fuel + TXNin_air(:,6).*Nin_air)+R3+(2*R1))./Nout;
Xout(:,6) = ((TXNin_fuel(:,7).*Nin_fuel + TXNin_air(:,7).*Nin_air))./Nout;
Xout(:,7) = ((TXNin_fuel(:,8).*Nin_fuel + TXNin_air(:,8).*Nin_air)-((1.5*R1)-(.5*R2)-(.5*R3)))./Nout;

[~,H_air] = enthalpy(Tin_air, Xin_air, Nin_air);
[~,H_fuel] =  enthalpy(Tin_fuel, Xin_fuel, Nin_fuel);

H_out = H_air + H_fuel-(R1.*hrxn1)-(R2.*hrxn2)-(R3.*hrxn3);

Tout = zeros(length(Tin_air),1)+1000;
T_error = 100;
while abs(T_error) > 1
   [~,H_guess] = enthalpy(Tout, Xout, Nout);
   Cp = SpecHeat(Tout, Xout);

   T_error = (H_out-H_guess)./(Cp.*Nout);
   Tout = Tout+T_error;
end





