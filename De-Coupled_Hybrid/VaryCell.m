ph = 1;
po = .21;
dVdT = linspace(.001, .010);
T1 =  linspace(400, 1500);
[Pnet, E, Real] = FuelCell(ph, po, dVdT, T1);
% figure(1)
% plot(dVdT,Pnet,'m--')
% figure(2)
% hold off
% [AX] = plotyy(dVdT,Real,dVdT,Pnet);
% ylabel(AX(1),'Efficiency')
% ylabel(AX(2),'Net Power')


figure(3)
surf(dVdT , T1, Real)