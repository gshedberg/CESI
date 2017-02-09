function [Efficiency,Eff_FC,Eff_GT,W_net,Wfc_vec,W_gt,Wc2,T_out,FCFlowOut,FC_Fuel_vec,combustorCH4,ReactMix,V_vec,Utilization,R_actual,Rt,recovery,Qextra,i_array,recirc_vec,nO2] = dFC_GT(varargin)
Tin = varargin{1}; %Temp, Composition, Flow in to system
Pr = varargin{2};       %Pressure Ratio across turbomachinary
P_ITMperm = varargin{3}; %Back pressure of OTM
TIT = varargin{4};%Ideal Turbine Inlet Temp
M_air = varargin{5}; %Mass Flow Rate of GT In
iDen = varargin{6};
Air.O2 = .21;
Air.N2 = .79;

if length(varargin)>6
    recovery = varargin{7};         %run at constant value of recovery
%     const = FALSE;
else recovery = ones(length(Tin),1).*.51; %Intial value of recovery
%     const = TRUE;  %will run loop to solve for recovery to match TIT
end
Cogen = 1; %Decide whether hydrogen will be generated
S2C = linspace(2,2,length(Pr))';    %design steam to carbon ratio
TurbEff = linspace(.88,.88,length(Pr))';      %Turbine Efficiency
CompEff = linspace(.8,.8,length(Pr))';      %Compressor Efficiency
LHVH2 = 240420; %Lower HEating Value of H2
vectorLength = max([length(Pr), length(P_ITMperm),length(TIT)]);
LHVfuel = zeros(vectorLength,1)+802301; %Lower heating value of CH4
TIT = zeros(length(TIT),1)+TIT;
%% Compressor
MMass = MassFlow(Air)/NetFlow(Air);
spec = fieldnames(Air);
for i = 1:1:length(spec)
    CompFlow.(spec{i}) = Air.(spec{i})*M_air/MMass;
end
CompFlow.T = Tin;
T1 = Tin;  %Initial conditions into hybrid
n = length(T1);
Pin = ones(n,1).*101;
[Wc1,T2,P2] = compress_struc(CompFlow,CompEff,Pr,Pin); %Compressor Model
CompFlow.T = T2;
error = 100;
while max(abs(error))> 1e-4% loop to solve for TIT by varying the recovery percentage of o2
    B = find(recovery>.99);
    if ~isempty(B)
       recovery(B) = 1;
    end
    %% OTM for Oxygen Separation
    [NonPerm,O2Flow,Q_preheat,R_actual,Rt] = ITM_struc(CompFlow,P2,P_ITMperm,recovery); %ITM Model
 %% Parasitic Compressor
    Pr2 = (Pr.*Pin+25)./P_ITMperm;      %Initialization for parasitic compressor to FC
    O2Flow.T = 320; %cooling back to near ambient conditions
    [Wc2,T5,P5] = compress_struc(O2Flow,CompEff,Pr2,P_ITMperm);   %Compressor2 Model
 %% Initialization to Fuel Cell
     L = 10;     %Length of cells
     W = 10;     %Width of cells
     nodes = 10;     %Nodes to analyze
     F = 96485;
     i_array = zeros(nodes,n);
     recirc_vec = zeros(n,1);
     FC_Fuel_vec = zeros(n,1);
     V_vec = zeros(n,1);
     Wfc_vec = zeros(n,1);
     Utilization = zeros(n,1);
     FCFlowOut=[];
     Q_HVanodeOut = zeros(n,1);
     FCFlowOut.T = zeros(n,1);
     FCFlowOut.H2 = zeros(n,1);
     FCFlowOut.H2O = zeros(n,1);
     FCFlowOut.CO = zeros(n,1);
     FCFlowOut.CO2 = zeros(n,1);
     FCFlowOut.CH4 = zeros(n,1);
     %% Fuel Cell Model
     for k = 1:1:n
        Supply.O2 = O2Flow.O2(k);
        %%# of cells will determine current density and thus voltage !!
%         Cells = .1231./(iDen(k)*L*W/(4000*F));
%         Cells = O2Flow.O2(k)./(iDen(k)*L*W/(4000*F));
        Cells = [510000];
        Cells_vec(k,:) = Cells;
        [i,recirc,FC_Fuel,FlowOut,V,W_fc,U] = FuelCell_struc(1023,.25,.6,S2C(k),Supply,L,W,nodes,Cells,Pr(k));
        i_array(:,k) = i;
        recirc_vec(k) = recirc;
        FC_Fuel_vec(k) = FC_Fuel;
        V_vec(k) = V;
        Wfc_vec(k) = W_fc;
        Utilization(k) = U;
        FCFlowOut.T(k) = FlowOut.T;
        FCFlowOut.H2(k) = FlowOut.H2;
        FCFlowOut.H2O(k) = FlowOut.H2O;
        FCFlowOut.CO(k) = FlowOut.CO;
        FCFlowOut.CO2(k) = FlowOut.CO2;
        FCFlowOut.CH4(k) = FlowOut.CH4;
        Q_HVanodeOut(k) = FlowOut.H2*241826.4 + FlowOut.CO*283004.7 + FlowOut.CH4*802301;%LHV of anode off products
     end
     %% Combustor Model
     [ReactMix, Qextra] = combust_struc(NonPerm,FCFlowOut,Q_preheat,TIT);
     combustorCH4 = Qextra./LHVfuel;
     combustorCH4(combustorCH4<0) = 0;
     H2co_produced = -Qextra./LHVH2;
     H2co_produced(H2co_produced<0) = 0;
     if length(varargin)<7   %% Adjust recovery to meet TIT
        a = .25;
        error = Qextra./(FC_Fuel*LHVfuel);
        newrecov = (1-a)*recovery + a*(recovery + error);
        recovery = min(newrecov,1);
        error(newrecov>=1) = 0;

     else %% Run at constant recovery
         error = 0;
     end
end
%% Turbine Model
[Wt,T_out,TurbFlow] = turbine_struc(ReactMix,TurbEff, 1./Pr);
%% Results
nO2 = O2Flow.O2;
W_gt = Wt-Wc1;  %Gas turbine power minus compression%
Eff_FC = Wfc_vec./(FC_Fuel_vec.*LHVfuel);
Eff_GT = (Wt-Wc1)./(Q_HVanodeOut + combustorCH4.*LHVfuel - H2co_produced.*LHVH2); %Efficiency of Gas Turbine
W_net = Wfc_vec + W_gt - Wc2;  %Net power output of hybrid
%% Decide whether cogeneration of h2 is wanted
if Cogen == 0
    Efficiency = (W_net)./((FC_Fuel_vec + combustorCH4).*LHVfuel); %Hybrid Efficiency without cogeneration
else
    Efficiency = (W_net+H2co_produced*LHVH2)./((FC_Fuel_vec + combustorCH4).*LHVfuel);%hybrid efficiency with cogeneration
end


   
    

