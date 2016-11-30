function Out = Oxidizer(t,Y, Inlet,block,string1)
% an oxidizer with n inlet flows 
% Complete combustion is assumed for all fuels (if enough O2 present), I want to add equilibrium calculation for for CO, CO2
% two or more (2+) inlets: Flows and outletl pressure
% Two (2) outlets: Flow and pressure at the inlet
% Two (2) states: Temperature, Inlet Pressure
global Ru Tags
Y = Y.*block.Scale;
Pin = Y(end);
n = block.inlets;
inlets = fieldnames(Inlet);
%merge flows
H_in = 0;
NetIn =[];

for j = 1:1:n
    H_in = H_in + enthalpy(Inlet.(inlets{j}));
    Spec = fieldnames(Inlet.(inlets{j}));
    for i = 1:1:length(block.spec)
        if ~isfield(NetIn,block.spec{i})
            NetIn.(block.spec{i}) = 0;
        end
        if ismember(block.spec{i},Spec)
            NetIn.(block.spec{i}) = NetIn.(block.spec{i}) + Inlet.(inlets{j}).(block.spec{i});
        end
    end
end
%% 3 reaction:
% CH4 + 1.5O2 --> CO + 2H2O
% CO + .5O2 --> CO2
% H2 + .5O2 --> H2O
R.CH4 = NetIn.CH4*block.Rxn1;
R.CO = (NetIn.CO+R.CH4)*block.Rxn2;
R.H2 = NetIn.H2*block.Rxn3;
sumR = 1.5*R.CH4 + 0.5*R.CO + 0.5*R.H2;
r = fieldnames(R);
phi = sumR/NetIn.O2;
if phi>1 %rich combustion
    for i = 1:1:length(r)
        R.(r{i}) = R.(r{i})/phi; %rich combustion limited by O2
    end
end

for i = 1:1:length(block.spec)
    if strcmp(block.spec{i},'CH4')
        ReactMix.CH4 = NetIn.CH4 - R.CH4;
    elseif strcmp(block.spec{i},'CO')
        ReactMix.CO = NetIn.CO + R.CH4- R.CO;
    elseif strcmp(block.spec{i},'CO2')
        ReactMix.CO2 = NetIn.CO2 + R.CO;
    elseif strcmp(block.spec{i},'H2')
        ReactMix.H2 = NetIn.H2  - R.H2;
    elseif strcmp(block.spec{i},'H2O')
        ReactMix.H2O = NetIn.H2O + 2*R.CH4 + R.H2;
    elseif strcmp(block.spec{i},'O2')
        ReactMix.O2 = NetIn.O2 - 1.5*R.CH4 - .5*R.CO - .5*R.H2;
    else
        ReactMix.(block.spec{i}) = NetIn.(block.spec{i});
    end
end
ReactMix.T = Y(1);

NetOut = block.Pfactor*(Pin-Inlet.Pout);%total cold flow out
Out.Flow.T = Y(1);
for i = 1:length(block.spec)
%     Out.Flow.(block.spec{i}) = ReactMix.(block.spec{i})*FlowFraction;%scale outflow by p-factor
    Out.Flow.(block.spec{i}) = Y(1+i);
end

if strcmp(string1,'Outlet')
    Out.Pin = Pin;
    Tags.(block.name).EquivelanceRatio = phi;
    Tags.(block.name).Temperatures = Y(1);
    Tags.(block.name).MassFlow = MassFlow(Out.Flow);
elseif strcmp(string1,'dY')
    dY = 0*Y;
    Cp = SpecHeat(Out.Flow);
    dY(1) = (H_in-enthalpy(ReactMix))*Ru*Y(1)/(Cp*block.Vol*Inlet.Pout);
    for i = 1:1:length(block.spec)
        dY(1+i) = (ReactMix.(block.spec{i}) - Y(1+i))*Ru*Y(1)/(block.Vol*Pin);
    end
    dY(2+i) = (NetFlow(ReactMix) - NetOut)*Ru*Y(1)/(block.Vol); %change in Pressure
    dY = dY./block.Scale;
    Out = dY;
end