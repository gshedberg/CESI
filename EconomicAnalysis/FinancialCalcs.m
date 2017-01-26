function [out, BarChart] = FinancialCalcs(Dispatch,Econ)
% First do the grid analysis for the baseline demand. Scale to the average 
% state electric rate of each month. Then use those scaled electric rates on 
% the dispatched generator grid demand. Then apply the financial calculations.

%% Financial Calculations
startYear = 2017;
Inflation = (1+Econ.Inflation/100);
Years = startYear:1:startYear+19;
inflationPriceFactor=(Inflation.^(Years-startYear));

steps = length(Dispatch.Elec);
% monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 365];
monthDays = [0 30];
hours = monthDays(end)*24; 
Ts = hours/steps;
if hours == 30*24;
    DE = Dispatch.Elec;
    DF = Dispatch.Fuel;
    Dispatch.Fuel = zeros(12,1);
    monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 365];
    for m= 1:1:length(monthDays)-1
        days = monthDays(m+1) - monthDays(m);
        if days<30
            Dispatch.Elec = [Dispatch.Elec; DE(1:days*24/Ts)];
            Dispatch.Fuel(m) = sum(DF(1:days*24/Ts))*Ts;
        elseif days == 30
            Dispatch.Elec = [Dispatch.Elec; DE];
            Dispatch.Fuel(m) = sum(DF)*Ts;
        else
            Dispatch.Elec = [Dispatch.Elec; DE; DE(1:(days-30)*24/Ts)];
            Dispatch.Fuel(m) =  (sum(DF)+sum(DF(1:(days-30)*24/Ts)))*Ts;
        end
    end
end

steps = length(Dispatch.Elec);
hours = monthDays(end)*24; 
tmx = [ones(steps,1)*[startYear 1 1]  linspace(0,hours-Ts,steps)' ones(steps,1)*[0 0]];
date = datenum(tmx);
month = datevec(date);
month = month(:,2);

gridEscalator(1) = 1;
for i = 2:1:Years(end)-startYear+1
    gridEscalator(i) = gridEscalator(i-1)*1.02;
end
gasEscalator = (Econ.AnnualNatGas/Econ.AnnualNatGas(1));
for i = 1:1:length(Years)
    FuelPrice(1+12*(i-1):12*i,1) = Econ.AnnualNatGas(i);
end

DemCharge = zeros(12,1);
MonthUseCost = zeros(12,1);
ReserveCharge = zeros(12,1);
%% Dispatch Grid Cost
for i = 1:12
    DemCharge(i) = max(Dispatch.Elec.*(month==i))*Econ.demandCharge;
    MonthUseCost(i) = sum(Econ.CentskWh/100.*Dispatch.Elec.*(month==i))*Ts;
    if max(Dispatch.Elec)>0
        ReserveCharge(i) = Econ.reserveCharge;
    end
end
out.Dispatch.TotalDemandCharges = sum(DemCharge);
out.Dispatch.TotalUseCharges = sum(MonthUseCost);
out.Dispatch.TotalGridBill = out.Dispatch.TotalDemandCharges+out.Dispatch.TotalUseCharges + sum(ReserveCharge);

%% Annual Costs & revenue
FuelCosts = sum(Dispatch.Fuel.*(FuelPrice(1:12)/293.01))*gasEscalator;%converts million BTU to kWh
GridBill = out.Dispatch.TotalGridBill*gridEscalator;
ReserveCharge = sum(ReserveCharge)*gridEscalator;
UseCharges = out.Dispatch.TotalUseCharges*gridEscalator;
DemCharges = out.Dispatch.TotalDemandCharges*gridEscalator;

%%%stack replacement
StackReplace = zeros(length(Years),1);
for i = 1:1:length(Years)
    lifespan = Econ.StackLife;
    if mod(i,lifespan)<1
        StackReplace(i) = Econ.StackReplaceCost*Dispatch.SysSize;
    else StackReplace(i) = 0;
    end
end
SystemCost = (Econ.InstallCost-Econ.Incentive)*Dispatch.SysSize +sum(inflationPriceFactor.*StackReplace') + Dispatch.BatterySize*Econ.Battery; 

%%Operation and Maintenance Costs
Total_OM = Econ.CHP_OM*Dispatch.SysSize + Econ.BatteryOM*Dispatch.BatterySize;
YearlyOandM = inflationPriceFactor*Total_OM;

InterestMonth = Econ.Interest/12/100;
months = Econ.FinanceYrs*12;
debtPaymentYearly=(12*(InterestMonth*SystemCost/(1-(1+InterestMonth)^(-months))))*(Years<startYear+Econ.FinanceYrs);

% Reserve capacity charge; demand charges, use charges, fuel cost, O & M, Finance
BarChart = [ReserveCharge(1); DemCharges(1); UseCharges(1); FuelCosts(1); YearlyOandM(1); debtPaymentYearly(1);];

out.InstalledCHPcost = SystemCost;
out.CostPerYear = debtPaymentYearly+YearlyOandM+FuelCosts+GridBill;
out.CostPerYearkW = out.CostPerYear/Dispatch.SysSize;
out.OandMFinance = debtPaymentYearly+YearlyOandM;

out.NPVnewGridCost=NPV(Inflation,GridBill);
out.NPVnewReserveCharge=NPV(Inflation,ReserveCharge);
out.NPVnewUseCharges=NPV(Inflation,UseCharges);
out.NPVnewDemCharges=NPV(Inflation,DemCharges);
out.NPVnewFuelCost=NPV(Inflation,FuelCosts);
out.NPVnewOandM = NPV(Inflation,YearlyOandM);
out.NPVnewFinance = NPV(Inflation,debtPaymentYearly);
out.NPVnewOandMandFinance = out.NPVnewOandM+out.NPVnewFinance;

out.actualYears = Years;
out.FuelCosts = FuelCosts;
out.NewGridBill = GridBill;

function value=NPV(irr,cashflows)
value=0;
for i=1:length(cashflows)
    value=value+cashflows(i)/(irr)^i;
end