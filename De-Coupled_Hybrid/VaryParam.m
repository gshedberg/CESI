Pr = linspace(3, 30);
Nt = 1;
Nc = 1;
T3 = linspace(800,1600);

% Wnet = zeros(100,1);
% Wt = Wnet;
% Wc = Wt;
% Eta = Wt;
% for j = 1:1:100
%     for i = 1:1:100
%         [Wnet(i,j), Wt(i,j), Wc(i,j), Eta(i,j)] = GasTurb(Pr(i), Nt, Nc, T3(j));
%     end
% end
figure(1)
hold off
plot(Pr,Wnet,'m--')
hold on
% plot(Pr,Wt,'g:','Linewidth',3)
% plot(Pr,Wc,'ro-')
ylabel('Pressure Ratio')
xlabel('Net Work')
% figure(2)
% hold off
% [AX] = plotyy(Pr,Eta,Pr,Wnet);
% ylabel(AX(1),'Efficiency')
% ylabel(AX(2),'Net Work')


figure(3)
surf(Pr,T3,Eta)
xlabel('Pressure Ratio')
ylabel('Inlet Turbine Temp')
zlabel('Net Efficiency')
