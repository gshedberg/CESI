function [Tout, Xout, Nout, vout, n_noz] = nozzle(TXNnp,TXNanodein, Pin)
Ru = 8.314;
Pex = zeros(length(Pin),1)+101;
Tex = zeros(length(Pin),1)+330;

Cp = SpecHeat(TXNin(:,1), TXNin(:,2:8));
gam = Cp./(Cp-Ru);

T_np = TXNnp(:,1);
X_np = TXNnp(:,2:8);
N_np = TXNnp(:,9);

[~,Hnp] = enthaply(T_np,X_np,N_np);

T_anode = TXNanodein(:,1);
X_anode = TXNanodein(:,2:8);
N_anode = TXNanodein(:,9);

Nout = N_np+N_anode;
Xout = X_anode+X_np;

[~,Hanode] = enthalpy(T_anode,X_anode,N_anode);

Hin = Hnp + Hanode;
Tout = zeros(length(Pin),1)+1000;
Terror = 100;
while abs(Terror) > .001
   [~,H_guess] = enthalpy(Tout, Xout, Nout);
   Cp = SpecHeat(T_out, X_out);

   Terror = (Hin-H_guess)./(Cp.*N_out);
   Tout = Tout+Terror;
end

[~,Hout] = enthalpy(Tout,Xout,Nout);

NPR = Pin./Pex;

n_noz = Hin./Hex;

vout = (2.*Cp.*Tout.*n_noz.*(1-(1./NPR)^(gam-1/gam)));


