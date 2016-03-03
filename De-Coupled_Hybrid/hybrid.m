function [Efficiency,Eff_FC,Eff_GT,W_net,W_fc,W_gt,T_out,X8,N_out,V_fc, R_actual,T8] = hybrid(TXNin, Pr, P_ITMperm, TIT,V_loss)
TurbEff = .87;
CompEff = .82;

T1 = TXNin(:,1);  %Initial conditions into hybrid
X1 = TXNin(:,2:8);
N1 = TXNin(:,9);
Pin = 101;

vectorLength = max([length(Pr), length(P_ITMperm),length(V_loss)]);
LHVfuel = zeros(vectorLength,1)+8e5; %Lower heating value of CH4

[Wc1, T2 ,X2,N2, P2] = compress([T1,X1,N1], CompEff, Pr, Pin);

[T7, X7, N7, N3, T3, Q_preheat, R_actual, Rt] = ITM([T2,X2,N2], P2, P_ITMperm,1);

X3 = zeros(length(T1),7);
X3(:,7) = 1;
P3=P_ITMperm;
Pr2 = (Pr*Pin+25)./P3;
T4 = T3*0+300;

X4 = X3;
N4 = N3;
[~,H3] = enthalpy(T3,X3,N3);
[~,H4] = enthalpy(T4,X4,N4);
Q34 = H3-H4;

% if P_ITMperm < Pin
[Wc2, T5, X5, N5, P5] = compress([T4, X4, N4], CompEff, Pr2, P3);

T6 = zeros(length(T1),1)+1023;
TXfuel = zeros(length(T1),8);
TXfuel(:,1) = 300;
TXfuel(:,2) = 1;

[V_fc] = nernst(Pr,T6,V_loss);
[X6,N6, W_fc, FC_fuel, Eff_FC] = FuelCell(V_fc,T6,[T5,X5,N5],TXfuel,2);


LHVCO = 303000;
LHVH2 = 240420;

    % [~,HfuelIn] = enthalpy(300,TXfuel(2:8),FC_fuel);
    % [~,HfuelOut] = enthalpy(T6,TXfuel(2:8),FC_fuel);
    % Error = FC_fuel.*LHVfuel-W_fc-N6.*(X6(:,1).*LHVfuel+X6(:,2).*LHVCO+X6(:,4).*LHVH2)+HfuelIn-HfuelOut;

LHVanode = zeros(length(T1),1)+LHVfuel.*X6(:,1)+LHVCO.*X6(:,2)+LHVH2.*X6(:,4);

TIT = zeros(length(T1),1)+TIT;

[X8, CombustFuel, N8, T8] = combust_mf([T7,X7,N7], [T6,X6,N6], TIT);

[Wt,T_out, N_out] = turbine([T8,X8,N8], TurbEff, 1./Pr);
W_gt = Wt-Wc1;
Eff_GT = (Wt-Wc1)./((N6.*LHVanode)+(CombustFuel.*LHVfuel)+Q_preheat);
W_net = W_fc + W_gt - Wc2;

Efficiency = W_net./(((FC_fuel+CombustFuel).*LHVfuel)+Q_preheat);


% else
%     X6 = [0 0 0 0 0 0 0];
%     
%     TIT = zeros(length(T1),1)+TIT;
% 
%     [X8, CombustFuel, N8, T8] = combust_mf([T7,X7,N7], [0,X6,0], TIT);
% 
%     [Wt,T_out, N_out] = turbine([T8,X8,N8], TurbEff, 1./Pr);
%     W_gt = Wt-Wc1;
%     Eff_GT = (Wt-Wc1)./((CombustFuel.*LHVfuel)+Q_preheat);
%     W_net = W_gt;
% 
%     Efficiency = W_net./((CombustFuel).*LHVfuel+Q_preheat);
%     
% end
