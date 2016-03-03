P_feed =  GasTurb(Pr);
P_ratio = linspace(1,10);
[Rt, XnpMat, Re, P_perm] = ITM(P_feed, P_ratio);

% Rt = zeros(100,1);
% XnpMat = Rt;
% Re = XnpMat;
% P_Perm = XnpMat;
% for j = 1:1:100
%     for i = 1:1:100
%         [Rt(i,j), XnpMat(i,j), Re(i,j), P_Perm(i,j)] = ITM(P_feed(i),P_ratio(j));
%     end
% end

figure(1)
hold on
plot(P_ratio,Re*100,'m--')
xlabel('Pressure Ratio of O2 across ITM')
ylabel('% Recovery of 02')
hold off
figure(2)
surf(P_ratio,Re*100,XnpMat)
xlabel('Pressure Ratio of O2 across ITM')
ylabel('% Recovery of 02')
zlabel('% Non-Permeate Air')