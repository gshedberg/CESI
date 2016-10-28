function [Efficiency,Eff_FC,Eff_GT,W_net,W_fc,W_gt,T_out,X8,N_out,V_fc,R_actual,recovery,FC_util,Qimbalance,E0,n_h2,Wc2,Rt] = hybrid_final(varargin)
TXNin = varargin{1}; %Temp, Composition, Flow in to system
Pr = varargin{2};       %Pressure Ratio across turbomachinary
P_ITMperm = varargin{3}; %Back pressure of OTM
V_loss =  varargin{4};      %Current Density across Fuel cell
TIT = varargin{5};          %Ideal Turbine Inlet Temp

if length(varargin)>5
    recovery = varargin{6};
else recovery = linspace(.4,.4)'; %Intial value of recovery
end

S2C = linspace(2,2,length(Pr))';
TurbEff = linspace(.90,.90,length(Pr))';      %Turbine Efficiency
CompEff = linspace(.85,.85,length(Pr))';      %Compressor Efficiency
LHVCO = 303000; %Lower Heating Value of CO
LHVH2 = 240420; %Lower HEating Value of H2

T1 = TXNin(:,1);  %Initial conditions into hybrid
X1 = TXNin(:,2:8);
N1 = TXNin(:,9);
Pin = 101;
TIT = zeros(length(T1),1)+TIT;

vectorLength = max([length(Pr), length(P_ITMperm),length(V_loss)]);
LHVfuel = zeros(vectorLength,1)+8e5; %Lower heating value of CH4
n_h2 = zeros(vectorLength,1);
[Wc1, T2 ,X2,N2, P2] = compress([T1,X1,N1], CompEff, Pr, Pin); %Compressor Model

if length(varargin)<6   %Condition to decide whether to run at constant recovery or adjust to meet TIT%
    n_h2 = 0;
    error = 1;
    while max(abs(error))> .99% loop to solve for TIT by varying the recovery percentage of o2
        B = find(recovery>.99);
        if ~isempty(B)
           recovery(B) = 1;
        end

        [T7, X7, N7, N3, T3, Q_preheat,R_actual,Rt] = ITM([T2,X2,N2], P2, P_ITMperm,recovery); %ITM Model
    
        X3 = zeros(length(T1),7);
        X3(:,7) = 1;
        P3=P_ITMperm;
        Pr2 = (Pr*Pin+25)./P3;      %Initialization for parasitic compressor to FC
        T4 = T3*0+300;
         X4 = X3;
        N4 = N3;
    
        [Wc2, T5, X5, N5, P5] = compress([T4, X4, N4], CompEff, Pr2, P3);   %Compressor2 Model

        T6 = zeros(length(T1),1)+1023;
        TXfuel = zeros(length(T1),8);
        TXfuel(:,1) = 300;              %Initialization of Fuel (CH4) inlet to FC
        TXfuel(:,2) = 1;

        [V_fc,E0] = nernst(Pr,T6,V_loss);       %Calculation of Open Circuit voltage based on Pressure
        [X6,N6, W_fc, FC_fuel, Eff_FC,FC_util,Qimbalance,I] = FuelCell(V_fc,T6,[T5,X5,N5],TXfuel,S2C);%FC model

        LHVanode = zeros(length(T1),1)+LHVfuel.*X6(:,1)+LHVCO.*X6(:,2)+LHVH2.*X6(:,4);%LHV of anode off products

        [T8, X8, N8,~] = combust([T7,X7,N7], [T6,X6,N6],Q_preheat,TIT); %No fuel combustor
        A = find(recovery>.99);
        if ~isempty(A)
            recovery(A) = 1;
            error = T8-TIT;
            error(A) = 0;
            newrecov = -error./TIT;
            recovery = (recovery + newrecov);%Adjust recovery percentage to get TIT to 1200%
        else
            error = T8-TIT;
            newrecov = -error./TIT;
            recovery = (recovery + newrecov);%Adjust recovery percentage to get TIT to 1200%
        end
    end
    CombustFuel = zeros(length(T1),1);
else %Run at constant recovery
    [T7, X7, N7, N3, T3, Q_preheat, R_actual, Rt] = ITM([T2,X2,N2], P2, P_ITMperm,recovery);%inputs recovery

    X3 = zeros(length(T1),7);
    X3(:,7) = 1;
    P3=P_ITMperm;
    Pr2 = (Pr*Pin+25)./P3;
    T4 = T3*0+300;

    X4 = X3;
    N4 = N3;

    [Wc2, T5, X5, N5, P5] = compress([T4, X4, N4], CompEff, Pr2, P3);

    T6 = zeros(length(T1),1)+1023;
    TXfuel = zeros(length(T1),8);
    TXfuel(:,1) = 300;
    TXfuel(:,2) = 1;

    [V_fc,E0] = nernst(Pr,T6,V_loss);
    [X6,N6, W_fc, FC_fuel, Eff_FC,FC_util,Qimbalance,I] = FuelCell(V_fc,T6,[T5,X5,N5],TXfuel,S2C);

    LHVanode = zeros(length(T1),1)+LHVfuel.*X6(:,1)+LHVCO.*X6(:,2)+LHVH2.*X6(:,4);
    [X8, CombustFuel, N8, T8,~,n_h2,Q_extra] = combust_mf([T7,X7,N7], [T6,X6,N6], TIT,Q_preheat);%Fuel injected (if needed) combustor% And H2 generation (extra)
%     if max(n_h2)~=0
%         A = find(n_h2 ~= 0);
%         x_h2 = zeros(length(X2),7);
%         x_h2(:,4) = 1;
%         [Tpre, X2pre, N2pre] = combust_h2([T2,X2,N2], [T6,x_h2,n_h2]);
%         [~,Hpre] = enthalpy(Tpre,X2pre,N2pre);
%         Q_extra = Hpre - Q_preheat;
%         Q_preheat(A) = 0;
%     end
end

[Wt,T_out, N_out] = turbine([T8,X8,N8], TurbEff, 1./Pr);
W_gt = Wt-Wc1;  %Gas turbine power minus parasitic compression%
Eff_GT = (Wt-Wc1)./((N6.*LHVanode)); %Efficiency of Gas Turbine
W_net = W_fc + W_gt - Wc2;  %Net power output of hybrid
Cogen = 0; %Decide whether cogeneration of h2 is wanted
if Cogen == 0
    Efficiency = (W_net)./((FC_fuel+CombustFuel).*LHVfuel-Qimbalance); %Hybrid Efficiency without cogeneration
else
    Efficiency = (W_net+(n_h2.*LHVH2))./((FC_fuel+CombustFuel).*LHVfuel-Qimbalance);%hybrid efficiency with cogeneration
end


   
    

