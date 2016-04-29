function [psi_destroyed]= exergy(TXNin,Win,Qin,TXNout,Wout,Qout)
T_0 = ones(length(TXNin(:,1)),1).*300;
[~,H_0] = enthalpy(T_0,TXNin(:,2:8),TXNin(:,9));
[~,S_0] = entropy(T_0,TXNin(:,2:8),TXNin(:,9));

Tin = TXNin(:,1);
Xin = TXNin(:,2:8);
Nin = TXNin(:,9);

[~,H_1] = enthalpy(Tin,Xin,Nin);
[~,S_1] = entropy(Tin,Xin,Nin);

Tout = TXNout(:,1);
Xout = TXNout(:,2:8);
Nout = TXNout(:,9);

psi_in = (H_1-H_0)-T_0.*(S_1-S_0);

[~,H_2] = enthalpy(Tout,Xout,Nout);
[~,S_2] = entropy(Tout,Xout,Nout);

psi_out = (H_2-H_0)-T_0.*(S_2-S_0);

psi_destroyed = (psi_in+Win+Qin) - (psi_out-Wout-Qout);