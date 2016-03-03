Pr = linspace(1,15);
P_ratio = linspace(15,30);
V_fc = linspace(.7,1.15);

Efficiency = zeros(100,1);
W_net = Efficiency;
% W_fc = W_net;
% W_gt = W_net;
% T_out = W_net;
% N_out = W_net



for j = 1:3
    for i = 1:3
        [Efficiency(i,j), W_net(i,j),W_fc(i,j),W_gt(i), T_out, X8,N_out] = hybrid([300, 0 0 0 0 0 .79 .21, .5], Pr(i), P_ratio(i), 1200,V_fc(j));
    end
end


% 
% [Efficiency, W_net,W_fc,W_gt, T_out, X8,N_out] = hybrid(TXNin, Pr, P_ratio, TIT,V_fc);
% 
figure(1)
plot(Pr,Efficiency)
ylabel('PressureRatio')
xlabel('Net Efficiency of de-coupled FC-GT')

% figure(2)
% plot(Pr,W_fc)
% ylabel('PressureRatio')
% xlabel('Net Power of Fuel Cell')
% 
% 
% figure(3)
% plot(Pr,W_gt)
% ylabel('PressureRatio')
% xlabel('Net Power of GT')
% 
% figure(4)
% plot(Pr,N_net)
% ylabel('PressureRatio')
% xlabel('Net Molar Flow')
% 
% figure(5)
% plot(Pr,H_out)
% ylabel('PressureRatio')
% xlabel('Enthalpy out')
% figure(2)
% 
% figure (6)
% surf(Pr,n_air,W_net)
% xlabel('PressureRatio')
% ylabel('Mol/sec')
% zlabel('Net Power Fc-GT')