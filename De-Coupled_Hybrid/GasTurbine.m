function [Tout, Xout, Nout, Eff_GT, W_GT] = GasTurbine
Pr = linspace(15,15)'; % Compressor pressure ratio
TIT = linspace(1200,1200)';
EffComp = linspace(1,1)';
EffTurb = linspace(1,1)';

Pin = linspace(101,101)';

vectorLength = max([length(Pr), length(EffComp),length(EffTurb)]);

TXNin = zeros(vectorLength,9);
TXNin(:,1) = 300;
Xin = [0 0 0 0 0 .79 .21];
TXNin(:,2:8) = Xin;
TXNin(:,9) = .5;
LHVfuel = zeros(vectorLength,1)+8e5; %Lower heating value of CH4

[Wc,T2,X2,N2,P2] = compress([TXNin(:,1),Xin,TXNin(:,9)], EffComp, Pr, Pin);

T3 = 300;
X3 = [1 0 0 0 0 0 0];
N3 = .005;

[X4,nfuel,N4,T4] = combust_mf([T2,X2,N2], [T3,X3,N3], TIT);


[Wt,T5, N5] = turbine([T4,X4,N4], EffTurb, Pr);

Tout = T5;
Xout = X4;
Nout = N5;
W_GT = Wt-Wc;
Eff_GT = W_GT./(LHVfuel.*(n_fuel+N3));
