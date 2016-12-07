function [ReactMix, Qextra] = combust_struc(Air,Fuel,Q,TIT)

NetIn =[];
Hair = enthalpy2(Air);
Hanode = enthalpy2(Fuel);
spec = fieldnames(Air);
spec = spec(~strcmp(spec,'T'));
for i = 1:1:length(spec)
    NetIn.(spec{i}) = Air.(spec{i});
end

spec2 = fieldnames(Fuel);
spec2 = spec2(~strcmp(spec2,'T'));
spec = unique([spec;spec2]);
for i = 1:1:length(spec2)
    if ~isfield(NetIn,spec2{i})
        NetIn.(spec2{i}) = Fuel.(spec2{i});
    else NetIn.(spec2{i}) = NetIn.(spec2{i}) + Fuel.(spec2{i});
    end
end


%% 3 reaction:
% CH4 + 1.5O2 --> CO + 2H2O
% CO + .5O2 --> CO2
% H2 + .5O2 --> H2O
R.CH4 = NetIn.CH4;
R.CO = (NetIn.CO+R.CH4);
R.H2 = NetIn.H2;
sumR = 1.5*R.CH4 + 0.5*R.CO + 0.5*R.H2;
r = fieldnames(R);
phi = sumR/NetIn.O2;
if phi>1 %rich combustion
    for i = 1:1:length(r)
        R.(r{i}) = R.(r{i})/phi; %rich combustion limited by O2
    end
end

for i = 1:1:length(spec)
    if strcmp(spec{i},'CH4')
        ReactMix.CH4 = NetIn.CH4 - R.CH4;
    elseif strcmp(spec{i},'CO')
        ReactMix.CO = NetIn.CO + R.CH4- R.CO;
    elseif strcmp(spec{i},'CO2')
        ReactMix.CO2 = NetIn.CO2 + R.CO;
    elseif strcmp(spec{i},'H2')
        ReactMix.H2 = NetIn.H2  - R.H2;
    elseif strcmp(spec{i},'H2O')
        ReactMix.H2O = NetIn.H2O + 2*R.CH4 + R.H2;
    elseif strcmp(spec{i},'O2')
        ReactMix.O2 = NetIn.O2 - 1.5*R.CH4 - .5*R.CO - .5*R.H2;
    else
        ReactMix.(spec{i}) = NetIn.(spec{i});
    end
end
Qextra = zeros(length(TIT),1);
if ~isempty(TIT)
    ReactMix.T=TIT;
    Hout = enthalpy2(ReactMix);
    Qextra = Hout - Hair  - (Hanode - Q); 
else
    ReactMix.T =  zeros(length(Air.T),1)+1000;
    T_error = 100;
    while min(abs(T_error) > .001)
       Hout = enthalpy2(ReactMix);
       Cp = SpecHeat2(ReactMix);
       T_error = (Hair + (Hanode - Q) - Hout)./(Cp.*NetFlow(ReactMix));
       ReactMix.T = ReactMix.T+T_error;
    end
end