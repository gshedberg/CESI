function [SOFC,SOC,Grid] =ControllerLoadFollow(DemandE,Battery,SysSize,RampRate,TurnDown)

steps = length(DemandE);
% monthDays = [0 31 59 90 120 151 181 212 243 273 304 334 365];
monthDays = [0 30];
hours = monthDays(end)*24; 
Ts = hours/steps;

%% Load SOFC parameters
CHPsize=sum(SysSize);
MinPower = sum(SysSize./TurnDown);
SOFC = zeros(steps,1);
SpareCap = zeros(steps,1);
Grid = zeros(steps,1);
SOC = zeros(steps,length(Battery));
Recharge = zeros(steps,1);
SOFC(1) = min(CHPsize,DemandE(1));
RampTs = sum(RampRate/100.*SysSize)*Ts;
%% Load Battery Parameters
%% Battery Characteristics
if ~isempty(Battery)
    for i = 1:1:length(Battery)
        Size(i)         = Battery(i).Size;
        PeakDisch(i)    = Battery(i).PeakDisch;
        PeakCharge(i)   = Battery(i).PeakCharge;
        DischResist(i)  = Battery(i).DischResist;
        ChargeResist(i) = Battery(i).ChargeResist;
        MaxDOD(i)       = Battery(i).MaxDOD;
        Voltage(i)      = Battery(i).Voltage;
        VoltCurveX(:,i) = Battery(i).VoltCurve(:,1);
        VoltCurveY(:,i) = Battery(i).VoltCurve(:,2);
    end
    DischCurrent = PeakDisch.*Size./Voltage*1000;
    DischResistScaled = (100./DischCurrent).*DischResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
    DischVoltLoss = DischCurrent.*DischResistScaled;
    PeakDischargePower = sum(DischCurrent.*(Voltage-DischVoltLoss))/1000;
    UseableSize  = Size.*(MaxDOD/100).*(Voltage-DischVoltLoss)./Voltage;
    totalUseSize = sum(UseableSize);
    ChargeCurrent = PeakCharge.*Size./Voltage*1000;
    ChargeResistScaled = (100./ChargeCurrent).*ChargeResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
    ChargeVoltLoss = ChargeCurrent.*ChargeResistScaled;
    ChargeEff = (Voltage-DischVoltLoss)./(Voltage+ChargeVoltLoss);
    PeakChargePower = sum(ChargeCurrent.*(Voltage+ChargeVoltLoss))/1000;

    SOC(1,:)                        = 1;
    Volts(1:length(Battery))        = 0;
    CurrentEst(1:length(Battery))   = 0;
    TerminalVolt(1:length(Battery)) = 0;
    Current(1:length(Battery))      = 0;
else
    PeakDischargePower = 0;
end
    
for t = 2:1:steps
    if DemandE(t) > SOFC(t-1)
        SOFC(t) = min([CHPsize SOFC(t-1)+RampTs DemandE(t)]);
    else
        SOFC(t) = max([MinPower SOFC(t-1)-RampTs DemandE(t)]);
    end 
    SpareCap(t) = max(0,min(CHPsize,SOFC(t-1)+RampTs)- DemandE(t)); %spare capacity to charge battery
    BatPower = max(0,min(PeakDischargePower,DemandE(t) - SOFC(t)));
    Grid(t) = max(0,DemandE(t) - SOFC(t) - BatPower);
    if ~isempty(Battery)
        if BatPower>0%use battery
            for i = 1:1:length(Battery)
                Volts(i) = interp1(VoltCurveX(:,i),VoltCurveY(:,i),SOC(t-1,i)*100);
                CurrentEst(i) = BatPower*(UseableSize(i)/totalUseSize)/Volts(i)*1000;
                TerminalVolt(i) = Volts(i)-CurrentEst(i)*DischResistScaled(i);
                Current(i) = BatPower*(UseableSize(i)/totalUseSize)/TerminalVolt(i)*1000;
                SOC(t,i) = SOC(t-1,i) -(Current(i)*Volts(i)/1000)*Ts/Size(i);
            end
        elseif min(SOC(t-1,:))<1 && SpareCap(t)>0 %charge battery if capacity exists
            ChargePow = min([SpareCap(t),sum((1-SOC(t-1,:))'.*Size),sum(PeakChargePower)]);
            for i = 1:1:length(Battery)
                Volts(i) = interp1(VoltCurveX(:,i),VoltCurveY(:,i),SOC(t-1,i)*100);
                if strcmp(Battery(i).ChargeMeth,'smooth') %Smoothest demand profile
                    CurrentEst(i) = ChargePow*(UseableSize(i)/totalUseSize)/Volts(i)*1000;
                    TerminalVolt(i) = Volts(i)+CurrentEst(i)*ChargeResistScaled(i);
                    Current(i) = ChargePow*(UseableSize(i)/totalUseSize)/TerminalVolt(i)*1000;
                elseif strcmp(Battery(i).ChargeMeth,'CC') || strcmp(Battery(i).ChargeMeth,'CV')
                    if strcmp(Battery(i).ChargeMeth,'CC') %Constant current
                        if sum(ChargeCurrent(1:i).*Volts(1:i)/1000)<Threshold-DemandE(t)
                            Current(i) = ChargeCurrent(i);
                        else Current(i) = 0;
                        end
                        TerminalVolt(i) = Volts(i)+Current(i)*ChargeResistScaled(i);
                    elseif strcmp(Battery(i).ChargeMeth,'CV') %Constant voltage
                        TerminalVolt(i) = interp1(VoltCurveX(:,i),VoltCurveY(:,i),100)+.2;% a voltage slightly higher than peak voltage to ensure 100% charging
                        Current(i) = (TerminalVolt(i)-Volts(i))/ChargeResistScaled(i);
                        if sum(Current(1:i).*TerminalVolt(1:i)/1000)>Threshold-DemandE(t)
                            Current(i) = 0;
                        end
                    end
                end
                ChargeFrac(i) = min((1-SOC(t-1,i))/((Current(i)*Volts(i)/1000)*Ts/Size(i)),1); %prevents overcharging battery
                SOC(t,i) = SOC(t-1,i) +(Current(i)*Volts(i)/1000)*Ts/Size(i)*ChargeFrac(i);
            end
            Recharge(t) = sum(TerminalVolt.*Current.*ChargeFrac)/1000;
            SOFC(t) = DemandE(t)+Recharge(t);
        else SOC(t,:) = 1;
        end
    end
    if isnan(SOFC(t))
        disp('WTF')
    end
end