function [Tout, Xout, Nout, Q] = combust(TXNin_air, TXNin_fuel,Q,TIT)


h = enthalpy(TXNin_air(:,1)+200);
LHVH2 = zeros(length(TIT),1)+240420; %Lower HEating Value of H2
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
R1 =(Xin_fuel(:,1).*(Nin_fuel));
R2 = (Xin_fuel(:,2).*(Nin_fuel)+R1);
R3 =(Xin_fuel(:,4).*(Nin_fuel));

Nout = (Nin_fuel + Nin_air + .5*R1 - .5*R2 - .5*R3);

Xout(:,1) = ((Xin_fuel(:,1).*Nin_fuel + Xin_air(:,1).*Nin_air) - R1)./Nout;
Xout(:,2) = ((Xin_fuel(:,2).*Nin_fuel + Xin_air(:,2).*Nin_air) + R1-R2)./Nout;
Xout(:,3) = ((Xin_fuel(:,3).*Nin_fuel + Xin_air(:,3).*Nin_air) + R2)./Nout;
Xout(:,4) = ((Xin_fuel(:,4).*Nin_fuel + Xin_air(:,4).*Nin_air) - R3)./Nout;
Xout(:,5) = ((Xin_fuel(:,5).*Nin_fuel + Xin_air(:,5).*Nin_air) + R3 + 2*R1)./Nout;
Xout(:,6) = ((Xin_fuel(:,6).*Nin_fuel + Xin_air(:,6).*Nin_air))./Nout;
Xout(:,7) = ((Xin_fuel(:,7).*Nin_fuel + Xin_air(:,7).*Nin_air) - 1.5*R1 - .5*R2 - .5*R3)./Nout;

% Xout(:,1) = ((Xin_fuel(:,1).*Nin_fuel + Xin_air(:,1).*Nin_air));
% Xout(:,2) = ((Xin_fuel(:,2).*Nin_fuel + Xin_air(:,2).*Nin_air));
% Xout(:,3) = ((Xin_fuel(:,3).*Nin_fuel + Xin_air(:,3).*Nin_air));
% Xout(:,4) = ((Xin_fuel(:,4).*Nin_fuel + Xin_air(:,4).*Nin_air));
% Xout(:,5) = ((Xin_fuel(:,5).*Nin_fuel + Xin_air(:,5).*Nin_air));
% Xout(:,6) = ((Xin_fuel(:,6).*Nin_fuel + Xin_air(:,6).*Nin_air));
% Xout(:,7) = ((Xin_fuel(:,7).*Nin_fuel + Xin_air(:,7).*Nin_air));

% Nout2 = sum(Xout,2);
% Xout = Xout/Nout2;

[~,H_air] = enthalpy(Tin_air, Xin_air, Nin_air);
[~,H_fuel] =  enthalpy(Tin_fuel, Xin_fuel, Nin_fuel); %Normal H_fuel w/out stripping hydrogen
Hrxn = (R1.*hrxn1)-(R2.*hrxn2)-(R3.*hrxn3);
H_out = H_air + H_fuel+ Hrxn-Q; %Normal enthalpy out without taking energy back to preheater
% H_out = H_air + H_fuel+((R1.*hrxn1)-(R2.*hrxn2)-(R3.*hrxn3)-Q);
% [~,H_TIT] = enthalpy(TIT,Xout,Nout);
% A = find(H_out<H_TIT);
% B = find(H_out>H_TIT);
% Q_extra = zeros(length(TIT),1);
% n_h2 = zeros(length(TIT),1);
% Q_dump = zeros(length(TIT),1);
% if ~isempty(B)
%     Q_extra(B) = H_out(B) - H_TIT(B);
%     n_h2(B) = Q_extra(B)./LHVH2(B); %amount of H2 that can be produced from extra energy
%     Nin_fuel(B) = Nin_fuel(B) - n_h2(B); %Subtract the hydrogen not needed and recycled to preheater
%     [~,H_fuel] =  enthalpy(Tin_fuel, Xin_fuel, Nin_fuel); %Normal H_fuel w/out stripping hydrogen
%     Q_dump(B) = Q_extra(B) - Q(B);
%     Q(B) = 0;
%     H_out = H_air + H_fuel+ Hrxn ; %Normal enthalpy out without taking energy back to preheater   
% end
Tout = zeros(length(Tin_air),1)+1000;
T_error = 100;
while min(abs(T_error) > .001)
   [~,H_guess] = enthalpy(Tout, Xout, Nout);
   Cp = SpecHeat(Tout, Xout);

   T_error = (H_out-H_guess)./(Cp.*Nout);
   Tout = Tout+T_error;
end





