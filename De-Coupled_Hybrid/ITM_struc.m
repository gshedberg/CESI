function [FlowOut,FlowO2,Q_preheat,R_actual,Rt] = ITM_struc(Flow,Pin,P_ITMperm,recovery)
%% Inputs to OTM
T = Flow.T;
X_feed.O2 = Flow.O2; %Feed composition from compressor
X_feed.N2 = Flow.N2;
% Nin = TXNin(:,9);
Hin = enthalpy2(Flow);       %Enthalpy in from compressor

%% Required Preheating
T_preheated = max(T,1073);
Flow.T = T_preheated;
Hheated = enthalpy2(Flow);
Q_preheat = (Hheated)-(Hin.O2);

%% Output Conditions
Tout = T_preheated;
TO2 = T_preheated;
FlowO2.T = TO2;
%% Recovery of Oxygen
Rt = 1-(1-X_feed.O2).*P_ITMperm./(X_feed.O2.*(Pin-P_ITMperm));%Theoretcial O2 recovery (Air Prod)
R_actual = max(0,(recovery) .* Rt); %(Actual O2 Recovery)
NO2 = (R_actual).*X_feed.O2.*Nin; %Molar Flow O2
FlowO2.O2 = NO2;
%% Non-permeate to Combustor
FlowOut = Flow-NO2; %Molar flow of Non-Permeate
FlowOut.T = Tout;
% for i = 1:1:7    
%     if i ==7
%         X_np(:,i) = (X_feed(:,i).*Nin - NO2)./Nout;%Composition of Non-Permeate
%     else
%         X_np(:,i) = (X_feed(:,i).*Nin )./Nout;%Composition of Non-Permeate
%     end
% end
end