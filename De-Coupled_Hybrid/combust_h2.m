function [Tout, Xout, Nout] = combust_h2(TXNin_air, TXNin_fuel)
erxn3 = 1;

h = enthalpy(TXNin_air(:,1)+200);
B = find(TXNin_fuel(:,9)== 0);
A = find(TXNin_fuel(:,9)~= 0);
if ~isempty(A) % h2 input
    Tin_air = TXNin_air(A,1);
    Xin_air = TXNin_air(A,2:8);
    Nin_air = TXNin_air(A,9);

    Tin_fuel = TXNin_fuel(A,1);
    Xin_fuel = TXNin_fuel(A,2:8);
    Nin_fuel = TXNin_fuel(A,9);
    %Reactions occuring in combustor%
    hrxn3 = h(:,5)-h(:,4)-.5*h(:,7); %H2 + .5 O2 -->  H2O

    %Total Enthalpy of combustion%
    R3 =(Xin_fuel(A,4).*(Nin_fuel))*erxn3;

    Nout = (Nin_fuel(A,1) + Nin_air(A,1) - .5*R3);

    Xout(A,1) = ((Xin_fuel(A,1).*Nin_fuel(A,1) + Xin_air(A,1).*Nin_air(A,1)))./Nout(A);
    Xout(A,2) = ((Xin_fuel(A,2).*Nin_fuel(A,1) + Xin_air(A,2).*Nin_air(A,1)))./Nout(A);
    Xout(A,3) = ((Xin_fuel(A,3).*Nin_fuel(A,1) + Xin_air(A,3).*Nin_air(A,1)))./Nout(A);
    Xout(A,4) = ((Xin_fuel(A,4).*Nin_fuel(A,1) + Xin_air(A,4).*Nin_air(A,1)) - R3)./Nout(A);
    Xout(A,5) = ((Xin_fuel(A,5).*Nin_fuel(A,1)+ Xin_air(A,5).*Nin_air(A,1)) + R3)./Nout(A);
    Xout(A,6) = ((Xin_fuel(A,6).*Nin_fuel(A,1)+ Xin_air(A,6).*Nin_air(A,1)))./Nout(A);
    Xout(A,7) = ((Xin_fuel(A,7).*Nin_fuel(A,1)+ Xin_air(A,7).*Nin_air(A,1)) - .5*R3)./Nout(A);

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
    [~,H_fuel] =  enthalpy(Tin_fuel, Xin_fuel, Nin_fuel);

    H_out = H_air(A,1) + H_fuel(A,1)-(R3.*hrxn3(A,1));

    Tout = zeros(length(Tin_air(A)),1)+800;
    T_error = 100;
    while max(abs(T_error) > .001)
    [~,H_guess] = enthalpy(Tout, Xout, Nout);
    Cp = SpecHeat(Tout, Xout);

    T_error = (H_out(A)-H_guess(A))./(Cp(A).*Nout(A));
    Tout(A) = Tout(A)+T_error(A);
    end
end
if ~isempty(B)
    Tout(B,1) = 0;
    Xout(B,1:7) = 0;
    Nout(B,1) = 0;
end