function [Thrust, Tout, Xout, Nout, Efficiency, Eff_FC] = hybridthrust(TXNin,Pr, P_ITMperm, n_h2, V_loss)

CompEff = .82;

T1 = TXNin(:,1);  %Initial conditions into hybrid
X1 = TXNin(:,2:8);
N1 = TXNin(:,9);
Pin = 101;

vectorLength = max([length(Pin), length(n_h2),length(V_loss)]);
LHVfuel = zeros(vectorLength,1)+2.42e5; %Lower heating value of H2

[P2, v2, T2, X2, N2] = diffuser(Pin, vin, [T1,X1,N1]);

[Wc1, T3 ,X3,N3, P3] = compress([T2,X2,N2], CompEff, Pr, P2);

[T4, X4, N4, N5, T5, Q_preheat] = ITM([T3,X3,N3], P3, P_ITMperm);

X5 = zeros(length(T1),7);
X5(:,7) = 1;
P4=P_ITMperm;
Pr2 = (Pr*Pin+25)./P4;
T4 = T3*0+300;

X6 = X5;
N6 = N5;
[~,H4] = enthalpy(T4,X4,N4);
[~,H5] = enthalpy(T5,X5,N5);
Q34 = H3-H4;

[Wc2, T5, X5, N5, P5] = compress([T5, X5, N5], CompEff, Pr2, P4);

T6 = zeros(length(T1),1)+1023;
TXfuel = zeros(length(T1),8);
TXfuel(:,1) = 300;
TXfuel(:,4) = 1;

[V_fc] = nernst(Pr,T6,V_loss);


[X7,N7, W_FC, n_fuel,Eff_FC] = FuelCellThrust(V_fc, T6,[T5,X5,N5],TXfuel);

[Tout, Xout, Nout, vout, n_noz] = nozzle([T4,X4,N4],[T6,X7,N7], P3);

W_net = W_FC - Wc1 - Wc2;
Efficiency = W_net./(n_fuel*LHVH2);

Thrust = ((Nout-N1)+n_fuel).*vout - N1.*vin;

