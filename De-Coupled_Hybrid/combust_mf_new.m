function [Xout, n_fuel, Nout, Tout] = combust_mf(TXNin_air, TXNin_fuel, TIT, Q_preheat)

erxn1 = 1;
erxn2 = 1;
erxn3 = 1;
Nin_air = TXNin_air(:,9);
Nin_fuel = TXNin_fuel(:,9);
h = enthalpy(TXNin_air(:,1)+200);

% Reactions occuring in combustor %
hrxn1 = 2*h(:,5)+h(:,2)-h(:,1)-1.5*h(:,7); %Ch4 + 1.5 O2 --> CO + 2 H2O
hrxn2 = h(:,3)-h(:,2)-.5*h(:,7); %CO + .5 O2 --> CO2 
hrxn3 = h(:,5)-h(:,4)-.5*h(:,7); %H2 + .5 O2 -->  H2O


[T_noFuel, x_noFuel, N_noFuel] = combust(TXNin_air, TXNin_fuel);
n_fuel = zeros(length(TXNin_air(:,1)),1);
Xout = x_noFuel;
Tout = T_noFuel;
Nout = N_noFuel;
AddFuel = (T_noFuel<TIT);
if min(T_noFuel)<TIT
    Q_hv = 8e5;
%     x_fuel = zeros(length(x_noFuel),7);
    x_fuel = zeros(length(TIT),7);
    x_fuel(:,1) = 1;
    T_fuel = zeros(length(Tout),1)+500;
    [~,H_TIT] = enthalpy(TIT, x_noFuel, N_noFuel);
    [~,H_noFuel] = enthalpy(T_noFuel, x_noFuel, N_noFuel);
    n_fuel = (H_TIT-H_noFuel)/(Q_hv).*AddFuel;
    error = 100;
    while max(abs(error))>.1
        R1 =(x_fuel(:,1).*n_fuel)*erxn1;
        R2 =((x_fuel(:,2).*n_fuel)+R1)*erxn2;
        R3 =(x_fuel(:,4).*n_fuel)*erxn3;

        Nout = N_noFuel + n_fuel + .5*R1 - .5*R2 - .5*R3;

        Xout(:,1) = ((x_noFuel(:,1).*N_noFuel + x_fuel(:,1).*n_fuel) - R1)./Nout;
        Xout(:,2) = ((x_noFuel(:,2).*N_noFuel ) + R1-R2)./Nout;
        Xout(:,3) = ((x_noFuel(:,3).*N_noFuel ) + R2)./Nout;
        Xout(:,4) = ((x_noFuel(:,4).*N_noFuel) - R3)./Nout;
        Xout(:,5) = ((x_noFuel(:,5).*N_noFuel ) + R3 + 2*R1)./Nout;
        Xout(:,6) = ((x_noFuel(:,6).*N_noFuel ))./Nout;
        Xout(:,7) = ((x_noFuel(:,7).*N_noFuel) - 1.5*R1 - .5*R2 - .5*R3)./Nout;
        
        [~,H_fuel] = enthalpy(T_fuel,x_fuel,n_fuel);
        H_out = H_noFuel+H_fuel-(R1.*hrxn1)-(R2.*hrxn2)-(R3.*hrxn3);
        
        Tout = max(TIT,T_noFuel);
        T_error = 100;
        while abs(T_error) > .1
           [~,H_guess] = enthalpy(Tout, Xout, Nout);
           Cp = SpecHeat(Tout, Xout);

           T_error = (H_out-H_guess)./(Cp.*Nout).*AddFuel;
           Tout = Tout+T_error.*AddFuel;
        end
        error = (TIT- Tout).*AddFuel;
        [~,H_TIT] = enthalpy(TIT,Xout,Nout);
        n_fuel = n_fuel + ((H_TIT-H_out).*AddFuel)./Q_hv;
    end
end
if T_noFuel > TIT && TIT == 1200;
    Q_extra = 0;
    Q_extra = H_TIT - H_noFuel;
    Q_preheat = Q_preheat - Q_extra;
    T_error = 100;
    while abs(min(T_error)) > 1
        [~,H_TIT] = enthalpy(TIT, x_noFuel, N_noFuel);
        [~,H_guess] = enthalpy(T_noFuel,Xout,Nout);
        Cp = SpecHeat(T_noFuel,X_out);
        T_error = (H_guess-H_TIT)./(Cp.*Nout);
        Tout = T_noFuel - T_error;
    end
end

if T_noFuel > TIT && TIT ~= 1200;
    

end
   



