
function [Tout, X_np, Nout, NO2, TO2, Q_preheat,R_actual,Rt] = ITM_heatex(TXNin, Pin,P_ITMperm, recovery,Q)
T = TXNin(:,1);
X_feed =TXNin(:,2:8); %Feed composition from compressor
Nin = TXNin(:,9);


[~,Hin] = enthalpy(T,X_feed,Nin);       %Enthalpy in from compressor
T_preheated = max(T,1073);
[~,Hheated] = enthalpy(T_preheated,X_feed,Nin);
Q_preheat = Hheated-Hin;
Tout = T_preheated;
TO2 = T_preheated;


Rt = 1-(1-X_feed(:,7)).*P_ITMperm./(X_feed(:,7).*(Pin-P_ITMperm));%Theoretcial O2 recovery (Air Prod)


R_actual = max(0,(recovery) .* Rt); %(Actual O2 Recovery)

NO2 = (R_actual).*X_feed(:,7).*Nin; %Molar Flow O2

Nout = Nin-NO2; %Molar flow of Non-Permeate

for i = 1:1:7    
    if i ==7
        X_np(:,i) = (X_feed(:,i).*Nin - NO2)./Nout;%Composition of Non-Permeate
    else
        X_np(:,i) = (X_feed(:,i).*Nin )./Nout;%Composition of Non-Permeate
    end
end
end