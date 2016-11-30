function [Efficiency,Eff_FC,Eff_GT,W_net,W_fc,W_gt,Wc2,T_out,FlowOut,ReactMix,V,R_actual,Rt,recovery,Qextra,i,recirc] = hybrid_final_struc(varargin)
Tin = varargin{1}; %Temp, Composition, Flow in to system
Pr = varargin{2};       %Pressure Ratio across turbomachinary
P_ITMperm = varargin{3}; %Back pressure of OTM
% V_loss =  varargin{4};      %Current Density across Fuel cell
TIT = varargin{4};%Ideal Turbine Inlet Temp
M_air = varargin{5}; %Mass Flow Rate of GT In
Air.O2 = .21;
Air.N2 = .79;

if length(varargin)>5
    recovery = varargin{6};
else recovery = linspace(.4,.4)'; %Intial value of recovery
end

S2C = linspace(2,2,length(Pr))';    %design steam to carbon ratio
TurbEff = linspace(.88,.88,length(Pr))';      %Turbine Efficiency
CompEff = linspace(.8,.8,length(Pr))';      %Compressor Efficiency
LHVCO = 303000; %Lower Heating Value of CO
LHVH2 = 240420; %Lower HEating Value of H2
vectorLength = max([length(Pr), length(P_ITMperm),length(TIT)]);
LHVfuel = zeros(vectorLength,1)+8e5; %Lower heating value of CH4
TIT = zeros(length(TIT),1)+TIT;
%% Compressor
MMass = MassFlow(Air)/NetFlow(Air);
spec = fieldnames(Air);
for i = 1:1:length(spec)
    CompFlow.(spec{i}) = Air.(spec{i})*M_air/MMass;
end
CompFlow.T = Tin;
T1 = Tin;  %Initial conditions into hybrid
Pin = ones(length(T1),1).*101;
[Wc1,T2,P2] = compress_struc(CompFlow,CompEff,Pr,Pin); %Compressor Model
CompFlow.T = T2;
%% Decide whether to run at constant recovery or adjust to meet TIT
if length(varargin)<6   
    error = 1;
    while max(abs(error))> .99% loop to solve for TIT by varying the recovery percentage of o2
        B = find(recovery>.99);
        if ~isempty(B)
           recovery(B) = 1;
        end
        %% OTM for Oxygen Separation
        [NonPerm,O2Flow,Q_preheat,R_actual,Rt] = ITM_struc(CompFlow,P2,P_ITMperm,recovery); %ITM Model
 %% Parasitic Compressor
        P3=P_ITMperm;
        Pr2 = (Pr.*Pin+25)./P3;      %Initialization for parasitic compressor to FC
        CompFlow2 = O2Flow;
        [Wc2,T5,P5] = compress_struc(CompFlow2,CompEff,Pr2,P3);   %Compressor2 Model
 %% Initialization to Fuel Cell
     T6 = 1023;
     ASR = 2;        %Design Area Specific Resistance
     eff_2 = .6;      %Effectiveness of WGS
     Oxidant.O2 = O2Flow.O2;
     L = 10;     %Length of cells
     W = 10;     %Width of cells
     n = 10;     %Nodes to analyze
     Cells = 1e5;
     i_array = zeros(n,length(Oxidant.O2));
     recirc_vec = zeros(length(Oxidant.O2),1);
     FC_Fuel_vec = zeros(length(Oxidant.O2),1);
     V_vec = zeros(length(Oxidant.O2),1);
     W_vec = zeros(length(Oxidant.O2),1);
     TFlowOut=[];
     TotalFlow = zeros(length(T1),1);
     TFlowOut.T = zeros(length(T1),1);
     TFlowOut.H2 = zeros(length(T1),1);
     TFlowOut.H2O = zeros(length(T1),1);
     TFlowOut.CO = zeros(length(T1),1);
     TFlowOut.CO2 = zeros(length(T1),1);
     TFlowOut.CH4 = zeros(length(T1),1);
     %% Fuel Cell Model
     for k = 1:1:length(Oxidant.O2)
         Supply.O2 = Oxidant.O2(k);
         S2C2 = S2C(k);
         Pr_FC = Pr(k);
        [i,recirc,FC_Fuel,FlowOut,V,W_fc] = FuelCell_struc(T6,ASR,eff_2,S2C2,Supply,L,W,n,Cells,Pr_FC);
        i_array(:,k) = i;
        recirc_vec(k) = recirc;
        FC_Fuel_vec(k) = FC_Fuel;
        V_vec(k) = V;
        W_vec(k) = W_fc;
        TFlowOut.T(k) = FlowOut.T;
        TFlowOut.H2(k) = FlowOut.H2;
        TFlowOut.H2O(k) = FlowOut.H2O;
        TFlowOut.CO(k) = FlowOut.CO;
        TFlowOut.CO2(k) = FlowOut.CO2;
        TFlowOut.CH4(k) = FlowOut.CH4;
        TotalFlow(k) = NetFlow(FlowOut);
     end
     %% Initialize Combustor
     LHVanode = zeros(length(T1),1)+(TFlowOut.H2./TotalFlow).*LHVH2+...
         (TFlowOut.CO./TotalFlow).*LHVCO+...
         (TFlowOut.CH4./TotalFlow).*LHVfuel;%LHV of anode off products
     %% Combustor Model
     [ReactMix, Qextra] = combust_struc(NonPerm,TFlowOut,Q_preheat,TIT);
%         [T8,X8,N8,~] = combust(T7,X7,N7, Tfuel,X6,N6,Q_preheat,TIT); %No fuel combustor
        A = find(recovery>.99);
        if ~isempty(A)
            recovery(A) = 1;
            error = ReactMix.T-TIT;
            error(A) = 0;
            newrecov = -error./TIT;
            recovery = (recovery + newrecov);%Adjust recovery percentage to get TIT to 1200%
        else
            error = ReactMix.T-TIT;
            newrecov = -error./TIT;
            recovery = (recovery + newrecov);%Adjust recovery percentage to get TIT to 1200%
        end
    end
else
    %% Run at constant recovery
 %% OTM for Oxygen Separation
        [NonPerm,O2Flow,Q_preheat,R_actual,Rt] = ITM_struc(CompFlow,P2,P_ITMperm,recovery); %ITM Model
 %% Parasitic Compressor
        P3=P_ITMperm;
        Pr2 = (Pr.*Pin+25)./P3;      %Initialization for parasitic compressor to FC
        CompFlow2 = O2Flow;
        [Wc2,T5,P5] = compress_struc(CompFlow2,CompEff,Pr2,P3);   %Compressor2 Model
 %% Initialization to Fuel Cell
     T6 = 1023;
     ASR = 2;        %Design Area Specific Resistance
     eff_2 = .6;      %Effectiveness of WGS
     Oxidant.O2 = O2Flow.O2;
     L = 10;     %Length of cells
     W = 10;     %Width of cells
     n = 10;     %Nodes to analyze
     Cells = 1e5;
     i_array = zeros(n,length(Oxidant.O2));
     recirc_vec = zeros(length(Oxidant.O2),1);
     FC_Fuel_vec = zeros(length(Oxidant.O2),1);
     V_vec = zeros(length(Oxidant.O2),1);
     W_vec = zeros(length(Oxidant.O2),1);
     TFlowOut=[];
     TotalFlow = zeros(length(T1),1);
     TFlowOut.T = zeros(length(T1),1);
     TFlowOut.H2 = zeros(length(T1),1);
     TFlowOut.H2O = zeros(length(T1),1);
     TFlowOut.CO = zeros(length(T1),1);
     TFlowOut.CO2 = zeros(length(T1),1);
     TFlowOut.CH4 = zeros(length(T1),1);
     %% Fuel Cell Model
     for k = 1:1:length(Oxidant.O2)
         Supply.O2 = Oxidant.O2(k);
         S2C2 = S2C(k);
         Pr_FC = Pr(k);
        [i,recirc,FC_Fuel,FlowOut,V,W_fc] = FuelCell_struc(T6,ASR,eff_2,S2C2,Supply,L,W,n,Cells,Pr_FC);
        i_array(:,k) = i;
        recirc_vec(k) = recirc;
        FC_Fuel_vec(k) = FC_Fuel;
        V_vec(k) = V;
        W_vec(k) = W_fc;
        TFlowOut.T(k) = FlowOut.T;
        TFlowOut.H2(k) = FlowOut.H2;
        TFlowOut.H2O(k) = FlowOut.H2O;
        TFlowOut.CO(k) = FlowOut.CO;
        TFlowOut.CO2(k) = FlowOut.CO2;
        TFlowOut.CH4(k) = FlowOut.CH4;
        TotalFlow(k) = NetFlow(FlowOut);
     end
     %% Initialize Combustor
     LHVanode = zeros(length(T1),1)+(TFlowOut.H2./TotalFlow).*LHVH2+...
         (TFlowOut.CO./TotalFlow).*LHVCO+...
         (TFlowOut.CH4./TotalFlow).*LHVfuel;%LHV of anode off products
     %% Combustor Model
     [ReactMix, Qextra] = combust_struc(NonPerm,TFlowOut,Q_preheat,TIT);
end
%% Turbine Model
[Wt,T_out] = turbine_struc(ReactMix,TurbEff, 1/Pr);
%% Results
W_gt = Wt-Wc1;  %Gas turbine power minus compression%
Eff_FC = W_vec./(FC_Fuel_vec.*LHVfuel);
Eff_GT = (Wt-Wc1)./((TotalFlow.*LHVanode)); %Efficiency of Gas Turbine
W_net = W_vec + W_gt - Wc2;  %Net power output of hybrid
%% Decide whether cogeneration of h2 is wanted
Cogen = 0; 
if Cogen == 0
    Efficiency = (W_net)./((FC_Fuel_vec.*LHVfuel)-Qextra); %Hybrid Efficiency without cogeneration
else
    Efficiency = (W_net+(FlowOut.H2.*LHVH2))./((FC_Fuel_vec.*LHVfuel)-Qextra);%hybrid efficiency with cogeneration
end


   
    

