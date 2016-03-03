
function [Tout, X_np, Nout, NO2, TO2, Q_preheat] = ITM(TXNin, Pin,P_ITMperm)
T = TXNin(:,1);
X_feed =TXNin(:,2:8); %Feed composition from compressor
Nin = TXNin(:,9);


[~,Hin] = enthalpy(T,X_feed,Nin);       %Enthalpy in from compressor

% P = X_feed(:,7).*(P-P_perm);

% X = (1-(1-X_feed(:,7))).*P_perm;

Rt = 1-(1-X_feed(:,7)).*P_ITMperm./(X_feed(:,7).*(Pin-P_ITMperm));%Theoretcial O2 recovery (Air Prod)

% Rt = (X./P)*10; 

R_actual = max(0,.5 * Rt); %(Actual O2 Recovery)

NO2 = R_actual.*X_feed(:,7).*Nin; %Molar Flow O2

Nout = Nin-NO2; %Molar flow of Non-Permeate

for i = 1:1:7    
    if i ==7
        X_np(:,i) = (X_feed(:,i).*Nin - NO2)./Nout;%Composition of Non-Permeate
    else
        X_np(:,i) = (X_feed(:,i).*Nin )./Nout;%Composition of Non-Permeate
    end
end

Tout = T*0 + 800+273;
TO2 = Tout;

[~,H_preheated] = enthalpy(Tout,X_feed,Nin);

Q_preheat = H_preheated-Hin;

end