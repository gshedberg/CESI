load Econ
load Time
load DemandE
load Battery
% N = length(DemandE);
% Ts = Time.Seconds(2)-Time.Seconds(1);
% Days = (Time.Seconds(end)+Ts)/(3600*24);
% fs = 1/Ts;
% f = (0:1/length(DemandE):(1/2-1/length(DemandE)))*fs;
% outlier = zeros(N,1);
% outScale =[];
% for t = 1:1:N
%     outlier(t) = DemandE(t)>(1.5*mean(DemandE(max(1,t-5):min(N,t+5))));
%     if outlier(t)
%         outScale(end+1) = DemandE(t)/mean(DemandE(max(1,t-5):min(N,t+5)));
%         DemandE(t) = mean(DemandE(max(1,t-5):min(N,t+5)));
%     end
% end
% d = round(rand(1)*N/(3600*24/Ts));
% t = max(1,d)*(3600*24/Ts); %shift by a random # of days
% altDemand = [DemandE(t:end);DemandE(1:t-1)];
% s(1) = rand(1)*.25+1;
% for d = 1:1:Days %stretch days by random amount(up to 25%)
%     s(d+1) = rand(1)*.25+1;
%     base = min(altDemand((d-1)*3600*24/Ts+1:d*3600*24/Ts));
%     altDemand((d-1)*3600*24/Ts+1:d*3600*24/Ts) = linspace(s(d),s(d+1),3600*24/Ts)'.*(altDemand((d-1)*3600*24/Ts+1:d*3600*24/Ts)-base) + base;
% end
% figure(1)
% hold off
% plot(Time.Days,DemandE)
% % xlim([0,4])
% hold on
% plot(Time.Days,altDemand,'r')
% 
% fdom = fft(altDemand,N);
% fdom(1) = 0;
% scaleFreq = ones(N,1);
% 
% r = 2*rand(N,1);
% r(r>1) = r(r>1).^2;
% scaleFreq(2:30*24) = r(2:30*24); %dont mess with frequencies <1hr
% scaleFreq(Days+1:Days:end) = 1; %same scale for things periodic with # of days
% for t = 2:1:(N/2)
%     scaleFreq(end-t+1) = scaleFreq(t);
% end
% 
% fdom2 = real(fdom).*scaleFreq + i*imag(fdom).*scaleFreq;
% figure(2)
% hold off
% plot(f*3600*24,abs(fdom(1:length(DemandE)/2)))
% xlim([0,5])
% hold on
% plot(f*3600*24,abs(fdom2(1:length(DemandE)/2)),'r')
% 
% newDemand = real(ifft(fdom2))+mean(DemandE);
% newDemand = newDemand - min(newDemand) + min(DemandE);
% outliers = find(outlier);
% for j = 1:1:2*length(outliers)
%     t = round(rand(1)*N);
%     b=max(1,round(rand(1)*length(outScale)));
%     newDemand(t) = (rand(1)*.4+.8)*newDemand(t)*outScale(b);
% end
% 
% 
% figure(3)
% hold off
% plot(Time.Days,DemandE)
% xlim([0,7])
% hold on
% plot(Time.Days,newDemand,'r')
% title('Comparison of microsoft trace and new trace')
% ylabel('% of rated demand');
% xlabel('Day of Month')
% 
% 
% newDemand = newDemand*100/max(newDemand);
% figure(4)
% plot(Time.Days,DemandE)
% xlim([0,7])
% title('Microsoft Trace')
% ylabel('% of rated demand');
% xlabel('Day of Month')
% 
% figure(5)
% plot(Time.Days,newDemand)
% xlim([0,7])
% title('New Trace')
% ylabel('% of rated demand');
% xlabel('Day of Month')

load ModifiedDemand
DemandE = newDemand;
tic
Ts = Time.Seconds(2) - Time.Seconds(1);
TsNew = 30;
newDemand = zeros(30*24*3600/TsNew,1);

for t = 1:1:length(DemandE)
    newDemand(Ts/TsNew*(t-1)+1:Ts/TsNew*t) = DemandE(t);
end
DemandE = newDemand;
Time.Seconds = linspace(0,30*24*3600-TsNew,30*24*3600/TsNew);
Time.Minutes = linspace(0,30*24*60-TsNew/60,30*24*3600/TsNew);
Time.Hours = linspace(0,30*24-TsNew/3600,30*24*3600/TsNew);
Time.Days = linspace(0,30-TsNew/24/3600,30*24*3600/TsNew);
toc
tic
%% baseline
DemandE = DemandE/100*1e4;
Dispatch.Elec = DemandE;
Dispatch.Fuel =  0*DemandE;
Dispatch.SysSize = 0;
Dispatch.BatterySize = 0;
Econ.reserveCharge = Econ.reserveChargeRate*1e4; %  $ for reserving capacity on the grid
[~, BarChart(:,1)] = FinancialCalcs(Dispatch,Econ);
figure(1)
plot(Time.Days,DemandE,'k');

toc
tic
%% Case A
SysSize = 2.5; %rated SOFC size
RampRate = 100; % change in power in %/hr
TurnDown = 6; % ratio of peak power to minimum power
SOFC_Eff = 0.65;
BatSize = .25;
Battery.Size = 100/Battery.MaxDOD*BatSize; %size in kWh;
Battery.PeakDisch =  15; %discharge rate
Battery.PeakCharge =  2; %charge rate
[SOFC,Bat,Grid] = ControllerLoadFollow(DemandE/1000,Battery,SysSize,RampRate,TurnDown);
Dispatch.Elec = Grid*1000;
Dispatch.Fuel = SOFC/SOFC_Eff*1000; %kW of fuel
Dispatch.SysSize = SysSize*1000;
Dispatch.BatterySize = BatSize*1000;
Econ.reserveCharge = 0; %  $ for reserving capacity on the grid
Econ.InstallCost = 3753;
Econ.StackReplaceCost = 354;
Econ.StackLife = 10;
[~, BarChart(:,2)] = FinancialCalcs(Dispatch,Econ);
figure(2)
AX = plotyy(Time.Days,SOFC,Time.Days,Bat);
hold on
plot(AX(1),Time.Days,DemandE/1000,'k');
ylim(AX(2),[0,1])

toc
tic
%% Case B
SysSize = 2.5; %rated SOFC size
RampRate = 100; % change in power in %/hr
TurnDown = 6; % ratio of peak power to minimum power
SOFC_Eff = 0.67;
BatSize = .25;
Battery.Size = 100/Battery.MaxDOD*BatSize; %size in kWh;
Battery.PeakDisch =  15; %discharge rate
Battery.PeakCharge =  2; %charge rate
[SOFC,Bat,Grid] = ControllerLoadFollow(DemandE/1000,Battery,SysSize,RampRate,TurnDown);
Dispatch.Elec = Grid*1000;
Dispatch.Fuel = SOFC/SOFC_Eff*1000; %kW of fuel
Dispatch.SysSize = SysSize*1000;
Dispatch.BatterySize = BatSize*1000;
Econ.reserveCharge = 0; %  $ for reserving capacity on the grid
Econ.InstallCost = 3500;
Econ.StackReplaceCost = 354;
Econ.StackLife = 10;
[~, BarChart(:,3)] = FinancialCalcs(Dispatch,Econ);
figure(3)
AX = plotyy(Time.Days,SOFC,Time.Days,Bat);
hold on
plot(AX(1),Time.Days,DemandE/1000,'k');
ylim(AX(2),[0,1])

toc
tic
%% Case C
SysSize = 2.5; %rated SOFC size
RampRate = 100; % change in power in %/hr
TurnDown = 6; % ratio of peak power to minimum power
SOFC_Eff = 0.66;
BatSize = .25;
Battery.Size = 100/Battery.MaxDOD*BatSize; %size in kWh;
Battery.PeakDisch =  15; %discharge rate
Battery.PeakCharge =  2; %charge rate
[SOFC,Bat,Grid] = ControllerLoadFollow(DemandE/1000,Battery,SysSize,RampRate,TurnDown);
Dispatch.Elec = Grid*1000;
Dispatch.Fuel = SOFC/SOFC_Eff*1000; %kW of fuel
Dispatch.SysSize = SysSize*1000;
Dispatch.BatterySize = BatSize*1000;
Econ.reserveCharge = 0; %  $ for reserving capacity on the grid
Econ.InstallCost = 2500;
Econ.StackReplaceCost = 354;
Econ.StackLife = 7;
[~, BarChart(:,4)] = FinancialCalcs(Dispatch,Econ);
figure(4)
AX = plotyy(Time.Days,SOFC,Time.Days,Bat);
hold on
plot(AX(1),Time.Days,DemandE/1000,'k');
ylim(AX(2),[0,1])

toc
tic
%% Case D
SysSize = 2000; %rated SOFC size
RampRate = 100; % change in power in %/hr
TurnDown = 5; % ratio of peak power to minimum power
SOFC_Eff = 0.6;
BatSize = 250;
Battery.Size = 100/Battery.MaxDOD*BatSize; %size in kWh;
Battery.PeakDisch =  15; %discharge rate
Battery.PeakCharge =  2; %charge rate
[SOFC,Bat,Grid] = ControllerLoadFollow(DemandE,Battery,SysSize,RampRate,TurnDown);
Dispatch.Elec = Grid;
Dispatch.Fuel = SOFC/SOFC_Eff; %kW of fuel
Dispatch.SysSize = SysSize;
Dispatch.BatterySize = BatSize;
Econ.reserveCharge = 0; %  $ for reserving capacity on the grid
Econ.InstallCost = 705;
Econ.StackReplaceCost = 196;
Econ.StackLife = 7;
[~, BarChart(:,5)] = FinancialCalcs(Dispatch,Econ);
figure(5)
AX = plotyy(Time.Days,SOFC,Time.Days,Bat);
hold on
plot(AX(1),Time.Days,DemandE,'k');
ylim(AX(2),[0,1])

toc
tic
%% Case Hybrid
SysSize = 1000; %rated SOFC size
RampRate = 100; % change in power in %/hr
TurnDown = 3; % ratio of peak power to minimum power
SOFC_Eff = 0.6;
BatSize = 0;
% Battery.Size = 100/Battery.MaxDOD*BatSize; %size in kWh;
% Battery.PeakDisch =  15; %discharge rate
% Battery.PeakCharge =  2; %charge rate
[SOFC,Bat,Grid] = ControllerLoadFollow(DemandE,[],SysSize,RampRate,TurnDown);
Dispatch.Elec = Grid;
Dispatch.Fuel = SOFC/SOFC_Eff; %kW of fuel
Dispatch.SysSize = SysSize;
Dispatch.BatterySize = BatSize;
Econ.reserveCharge = Econ.reserveChargeRate*1e4; %  $ for reserving capacity on the grid
Econ.InstallCost = 705;
Econ.StackReplaceCost = 196;
Econ.StackLife = 6;
[~, BarChart(:,6)] = FinancialCalcs(Dispatch,Econ);
figure(6)
plot(Time.Days,SOFC);
hold on
plot(Time.Days,DemandE,'k');
toc
%plot together
figure(7)
bar(BarChart','stacked')
legend('Reserve capacity charge','Demand charges', 'Use charges', 'Fuel cost', 'O & M', 'Finance')
ylabel('$100,000 per year')
