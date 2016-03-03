function [Xout, n_fuel, Nout, Tout,Q_recyc,n_h2,Q_extra] = combust_mf(TXNin_air, TXNin_fuel, TIT,Q_preheat)
erxn1 = 1;
erxn2 = 1;
erxn3 = 1;
Nin_air = TXNin_air(:,9);
Nin_fuel = TXNin_fuel(:,9);
h = enthalpy(TXNin_air(:,1)+200);
LHVH2 = zeros(length(TIT),1)+240420; %Lower HEating Value of H2
% Reactions occuring in combustor %
hrxn1 = 2*h(:,5)+h(:,2)-h(:,1)-1.5*h(:,7); %Ch4 + 1.5 O2 --> CO + 2 H2O
hrxn2 = h(:,3)-h(:,2)-.5*h(:,7); %CO + .5 O2 --> CO2 
hrxn3 = h(:,5)-h(:,4)-.5*h(:,7); %H2 + .5 O2 -->  H2O


[T_noFuel, x_noFuel, N_noFuel,Q_recyc] = combust(TXNin_air, TXNin_fuel,Q_preheat,TIT);
n_fuel = zeros(length(TXNin_air(:,1)),1);
n_h2 = zeros(length(TIT),1);
Xout = x_noFuel;
Tout = T_noFuel;
Nout = N_noFuel;
x_fuel = zeros(length(TIT),7);
A = find(T_noFuel<TIT);
B = find(T_noFuel>TIT);
H_out = zeros(length(TIT),1);
if ~isempty(A)%min(T_noFuel)<TIT
    Q_hv = 8e5;
    x_fuel(A,1) = 1;
    T_fuel = zeros(length(TIT),1)+500;
    [~,H_TIT] = enthalpy(TIT, x_noFuel, N_noFuel);
    [~,H_noFuel] = enthalpy(T_noFuel, x_noFuel, N_noFuel);
    n_fuel(A) = (H_TIT(A)-H_noFuel(A))/(Q_hv);
    error = 100;
    while max(abs(error))>.1
        R1 =(x_fuel(A,1).*n_fuel(A,1))*erxn1;
        R2 =((x_fuel(A,2).*n_fuel(A,1)+R1))*erxn2;
        R3 =(x_fuel(A,4).*n_fuel(A,1))*erxn3;

        Nout(A) = N_noFuel(A,1) + n_fuel(A,1) + .5*R1 - .5*R2 - .5*R3;

        Xout(A,1) = ((x_noFuel(A,1).*N_noFuel(A,1) + x_fuel(A,1).*n_fuel(A,1)) - R1)./Nout(A);
        Xout(A,2) = ((x_noFuel(A,2).*N_noFuel(A,1)) + R1-R2)./Nout(A);
        Xout(A,3) = ((x_noFuel(A,3).*N_noFuel(A,1)) + R2)./Nout(A);
        Xout(A,4) = ((x_noFuel(A,4).*N_noFuel(A,1)) - R3)./Nout(A);
        Xout(A,5) = ((x_noFuel(A,5).*N_noFuel(A,1)) + R3 + 2*R1)./Nout(A);
        Xout(A,6) = ((x_noFuel(A,6).*N_noFuel(A,1)))./Nout(A);
        Xout(A,7) = ((x_noFuel(A,7).*N_noFuel(A,1)) - 1.5*R1 - .5*R2 - .5*R3)./Nout(A);
        [~,H_fuel] = enthalpy(T_fuel,x_fuel,n_fuel);
        H_out(A) = H_noFuel(A,1)+H_fuel(A,1)-(R1.*hrxn1(A,1))-(R2.*hrxn2(A,1))-(R3.*hrxn3(A,1));
        
        Tout(A) = TIT(A);
        [~,H_TIT] = enthalpy(TIT,Xout,Nout);
        error = (H_TIT(A) - H_out(A))./(SpecHeat(TIT(A),Xout(A,:)).*Nout(A));
        
        n_fuel(A) = n_fuel(A) + (H_TIT(A) - H_out(A))./Q_hv;
        A = A(find(abs(error)>.1));
    end
end
Q_extra = zeros(length(TIT),1);
if ~isempty(B) %min(T_noFuel)>TIT
    [~,H_TIT] = enthalpy(TIT,Xout,Nout);
    [~,H_noFuel] = enthalpy(T_noFuel, x_noFuel, N_noFuel);
    Q_extra(B) = H_noFuel(B)-H_TIT(B);
%     Q_preheat(B) = Q_preheat(B) - Q_extra(B);
    Tout(B) = TIT(B);
    n_h2(B) = Q_extra(B)./LHVH2(B);
end
