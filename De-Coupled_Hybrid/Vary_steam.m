P_ITMperm = linspace(50,50)'; %ITM back pressure in kPa
V_loss = linspace(.21,.21)'; %Fuel cell voltage
recovery = linspace(.4,.4)'; %Fixed value of recovery
TIT = linspace(1200,1200)';
Fuel = 0; % 0 for no suplemental fuel into combustor, 1 for fixed % recovery
Pr = linspace(5,25)'; % Compressor pressure ratio
vectorLength = max([length(Pr), length(P_ITMperm),length(V_loss)]); %set length of vectors to correspond to given inputs

TXNin = zeros(vectorLength,9);
TXNin(:,1) = 300;
Xin = [0 0 0 0 0 .79 .21];      %initial conditions of temp, x, and flow
for i =1:1:7
    TXNin(:,i+1) = Xin(i);
end
TXNin(:,9) = .5;


[Efficiency,Eff_FC,Eff_GT,W_net,W_fc,W_gt,T_out,X8,N9,V_fc,R_actual,recovery,FC_util,Qimbalance,E0,n_h2] = hybrid_final_steam(varargin);


